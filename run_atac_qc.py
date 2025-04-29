#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import argparse
import subprocess
from pathlib import Path
import sys
import shutil # 用于检查命令是否存在
import glob # 用于查找文件
import pandas as pd # For creating FRiP summary table

def find_executable(executable):
    """检查可执行文件是否存在于系统的 PATH 中。(Checks if an executable exists in the system's PATH.)"""
    return shutil.which(executable) is not None

def parse_args():
    """
    解析命令行参数 (Parses command-line arguments)
    """
    parser = argparse.ArgumentParser(
        description="ATAC-seq 质量控制脚本：使用主流程的输出生成 QC 图和指标。"
                    "(ATAC-seq Quality Control Script: Generates QC plots and metrics using outputs from the main pipeline.)"
    )
    # --- Input ---
    parser.add_argument(
        "-b", "--base_dir",
        required=True,
        help="主 ATAC-seq 流程的根输出目录路径。(Path to the base output directory of the main ATAC-seq pipeline.)"
    )
    parser.add_argument(
        "--tss_bed",
        required=True,
        help="[TSS富集需要/TSS Enrichment requires] 提供一个包含TSS坐标的BED文件路径。(Path to a BED file containing TSS coordinates.)"
    )
    # --- Output ---
    parser.add_argument(
        "-o", "--output_dir",
        required=True,
        help="保存 QC 结果的输出目录路径。(Path to the output directory for saving QC results.)"
    )
    # --- Tool Settings ---
    parser.add_argument(
        "--cores",
        type=int,
        default=4,
        help="执行工具 (如 deepTools, samtools) 时使用的线程数。(Number of threads for tools like deepTools, samtools) (默认/Default: 4)"
    )
    # --- Parallelization ---
    parser.add_argument(
        "--parafly_jobs",
        type=int,
        default=4,
        help="使用 ParaFly 并行运行每个样本的任务数 (例如片段大小计算)。(Number of per-sample jobs to run in parallel using ParaFly, e.g., for fragment size calculation) (默认/Default: 4)"
    )
    # --- Execution Control ---
    parser.add_argument(
        "--no_run",
        action="store_true",
        help="如果指定，则只打印命令而不实际执行，用于调试。(If specified, only print commands without actual execution, for debugging.)"
    )
    parser.add_argument(
        "--skip_fragsize",
        action="store_true",
        help="如果指定，则跳过片段大小分布计算。(If specified, skip fragment size distribution calculation.)"
    )
    parser.add_argument(
        "--skip_tss",
        action="store_true",
        help="如果指定，则跳过 TSS 富集图生成。(If specified, skip TSS enrichment plot generation.)"
    )
    parser.add_argument(
        "--skip_frip",
        action="store_true",
        help="如果指定，则跳过 FRiP 分数计算。(If specified, skip FRiP score calculation.)"
    )

    return parser.parse_args()

def find_files(base_dir, pattern):
    """在指定的基础目录下查找匹配模式的文件，返回绝对路径列表。
    (Finds files matching a pattern within a base directory, returns list of absolute paths.)"""
    found_files = sorted([str(p.resolve()) for p in Path(base_dir).glob(pattern)])
    return found_files

def run_single_command(step_name, cmd, output_dir, no_run=False, check_error=True, executable_name=None, capture_stdout=False):
    """
    执行单个命令。如果 capture_stdout=True，则返回 stdout 内容。
    (Executes a single command. Returns stdout content if capture_stdout=True.)
    """
    print(f"\n--- 运行步骤 (Running Step): {step_name} ---")
    print(f"[命令/CMD] {cmd}")

    # Check executable existence
    if not no_run:
        tool_to_check = executable_name or cmd.split()[0]
        if not find_executable(tool_to_check):
            print(f"[错误/ERROR] 命令 '{tool_to_check}' 未找到。")
            if check_error:
                print(f"[错误/ERROR] 因缺少命令而无法执行步骤 {step_name}。终止。")
                sys.exit(1)
            else:
                print(f"[警告/WARN] 因缺少命令而跳过步骤 {step_name}。")
                return False, None if capture_stdout else False # Return failure indication

    # Execute command
    if not no_run:
        try:
            os.makedirs(output_dir, exist_ok=True)
            process = subprocess.Popen(cmd, shell=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=output_dir)
            stdout, stderr = process.communicate()
            exit_code = process.returncode

            # Always print stderr for context
            if stderr: print(f"[{step_name} STDERR]\n{stderr.strip()}")
            # Print stdout only if not capturing or if it's short
            if stdout and (not capture_stdout or len(stdout) < 500):
                 print(f"[{step_name} STDOUT]\n{stdout.strip()}")
            elif stdout and capture_stdout:
                 print(f"[{step_name} STDOUT] (Captured, {len(stdout)} bytes)")


            if exit_code != 0:
                raise subprocess.CalledProcessError(exit_code, cmd, output=stdout, stderr=stderr)

            print(f"[信息/INFO] 步骤 {step_name} 成功执行。")
            return True, stdout.strip() if capture_stdout else True # Return success and captured output if requested

        except subprocess.CalledProcessError as e:
            print(f"[错误/ERROR] 步骤 {step_name} 失败，退出代码 {e.returncode}")
            print(f"    命令 (Command): {e.cmd}")
            if check_error:
                print(f"[错误/ERROR] 因步骤 {step_name} 失败而退出。")
                sys.exit(1)
            else:
                print(f"[警告/WARN] 步骤 {step_name} 失败，但由于 check_error=False，执行将继续。")
                return False, None if capture_stdout else False # Return failure
        except Exception as e:
            print(f"[错误/ERROR] 步骤 {step_name} 执行期间发生意外错误: {e}")
            if check_error:
                sys.exit(1)
            else:
                return False, None if capture_stdout else False # Return failure
    else:
        print(f"[信息/INFO] --no_run 已指定。步骤 {step_name} 未执行。")
        # Return success indication for --no_run, but None for captured output
        return True, None if capture_stdout else True

def run_parafly_step(step_name, command_list, command_file_path, log_file_path, cpu, output_dir, no_run=False):
    """
    写入命令文件并执行ParaFly。
    (Writes command file and executes ParaFly.)
    """
    if not command_list:
        print(f"[信息/INFO] 步骤 {step_name} 没有生成任何命令，跳过 ParaFly 执行。")
        return True

    print(f"\n--- 准备并行步骤 (Preparing Parallel Step): {step_name} ---")
    print(f"    写入 {len(command_list)} 条命令到: {command_file_path}")
    try:
        os.makedirs(command_file_path.parent, exist_ok=True)
        with open(command_file_path, 'w') as f:
            f.write("\n".join(command_list) + "\n")
    except IOError as e:
        print(f"[错误/ERROR] 无法写入命令文件 {command_file_path}: {e}")
        return False

    if not no_run and not find_executable("ParaFly"):
         print(f"[错误/ERROR] 命令 'ParaFly' 未找到。")
         print(f"[错误/ERROR] 无法执行并行步骤 {step_name}。")
         # Don't exit for QC steps, just return failure
         return False

    cmd_parafly = (
        f"ParaFly -c {command_file_path} -CPU {cpu} "
        f"-failed_cmds {log_file_path} -v"
    )

    print(f"--- 运行 ParaFly 步骤 (Running ParaFly for Step): {step_name} ---")
    print(f"[命令/CMD] {cmd_parafly}")
    if not no_run:
        try:
            os.makedirs(output_dir, exist_ok=True)
            process = subprocess.Popen(cmd_parafly, shell=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=output_dir)
            stdout, stderr = process.communicate()
            exit_code = process.returncode

            if stdout: print(f"[ParaFly STDOUT]\n{stdout.strip()}")
            if stderr: print(f"[ParaFly STDERR]\n{stderr.strip()}")

            failed_log_exists = log_file_path.exists() and log_file_path.stat().st_size > 0

            if exit_code != 0:
                if failed_log_exists:
                    print(f"[错误/ERROR] ParaFly 在步骤 {step_name} 中失败 (exit code {exit_code}) 并记录了失败的命令。")
                    print(f"    请检查日志文件: {log_file_path}")
                else:
                    print(f"[错误/ERROR] ParaFly 在步骤 {step_name} 中失败 (exit code {exit_code}) 但未记录失败的命令。")
                return False
            else:
                if failed_log_exists:
                    print(f"[警告/WARN] ParaFly 步骤 {step_name} 完成，但报告了失败的命令。")
                    print(f"    请检查日志文件: {log_file_path}")
                else:
                    print(f"[信息/INFO] ParaFly 步骤 {step_name} 成功完成。")
                return True

        except Exception as e:
            print(f"[错误/ERROR] 运行 ParaFly 步骤 {step_name} 时发生意外错误: {e}")
            return False
    else:
        print(f"[信息/INFO] --no_run 已指定。ParaFly 命令步骤 {step_name} 未执行。")
        return True

def main():
    args = parse_args()
    base_dir = Path(args.base_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    qc_log_dir = output_dir / "QC_Scripts_Logs" # Specific log dir for this script

    # --- Create Output and Log Directories ---
    try:
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(qc_log_dir, exist_ok=True)
        print(f"[信息/INFO] QC 输出目录 (QC Output directory): {output_dir}")
        print(f"[信息/INFO] QC 日志目录 (QC Log directory): {qc_log_dir}")
    except OSError as e:
        print(f"[错误/ERROR] 无法创建输出或日志目录: {e}")
        sys.exit(1)

    # --- Define expected input subdirectories and file patterns ---
    bam_dir_pattern = "Run4_Filtered_BAM" # Use final filtered BAMs
    peaks_dir_pattern = "Run7_Peaks"
    bw_dir_pattern = "Run8_BigWig" # Needed for TSS plot consistency

    # --- Find Input Files ---
    print("\n--- 查找输入文件 (Finding Input Files) ---")
    bam_files = find_files(base_dir, f"{bam_dir_pattern}/*.last.bam")
    peak_files = find_files(base_dir, f"{peaks_dir_pattern}/*_peaks.narrowPeak")
    bw_files = find_files(base_dir, f"{bw_dir_pattern}/*.last.bw") # Find bigwigs for TSS plot

    # Create a mapping from sample name to files (assuming consistent naming)
    sample_files = {}
    for bam_path in bam_files:
        sample_name = Path(bam_path).name.replace("_Run4.last.bam", "")
        sample_files[sample_name] = {"bam": bam_path}

    for peak_path in peak_files:
        sample_name = Path(peak_path).name.replace("_Run7_peaks.narrowPeak", "")
        if sample_name in sample_files:
            sample_files[sample_name]["peak"] = peak_path
        else:
            print(f"[警告/WARN] 找到 Peak 文件但没有对应的 BAM 文件: {peak_path}")

    for bw_path in bw_files:
        sample_name = Path(bw_path).name.replace("_Run8.last.bw", "")
        if sample_name in sample_files:
             sample_files[sample_name]["bw"] = bw_path
        # else: # Don't warn if bigwig is missing, TSS plot handles it
        #     print(f"[警告/WARN] Found BigWig file without corresponding BAM: {bw_path}")


    if not sample_files:
        print(f"[错误/ERROR] 在 '{base_dir / bam_dir_pattern}' 中未找到 BAM 文件 (*.last.bam)。无法继续 QC。")
        sys.exit(1)
    else:
         print(f"  找到 {len(sample_files)} 个样本的 BAM 文件。")

    # --- Resolve TSS BED Path ---
    tss_bed_path_abs = None
    if not args.skip_tss:
        try:
            tss_bed_path_abs = Path(args.tss_bed).resolve(strict=True) # strict=True checks existence
            print(f"[信息/INFO] 使用 TSS BED 文件进行富集分析: {tss_bed_path_abs}")
        except FileNotFoundError:
            print(f"[错误/ERROR] 提供的 TSS BED 文件未找到: {args.tss_bed}。跳过 TSS 富集分析。")
            args.skip_tss = True # Force skip if file not found
        except Exception as e:
            print(f"[错误/ERROR] 处理 TSS BED 文件路径时出错 '{args.tss_bed}': {e}。跳过 TSS 富集分析。")
            args.skip_tss = True # Force skip on other errors


    # =========================================================================
    #                          QC STEPS START
    # =========================================================================
    all_qc_successful = True # Track overall success

    # ==================== 1. Fragment Size Distribution ====================
    if not args.skip_fragsize:
        print("\n--- 开始计算片段大小分布 (Starting Fragment Size Distribution Calculation) ---")
        step_key = "frag_size"
        fragsize_output_dir = output_dir / "Fragment_Sizes"
        dt_frag_exe = "bamPEFragmentSize"

        if not find_executable(dt_frag_exe):
            print(f"[警告/WARN] 跳过片段大小计算，因为 '{dt_frag_exe}' 未在 PATH 中找到。")
            all_qc_successful = False
        else:
            os.makedirs(fragsize_output_dir, exist_ok=True)
            commands = []
            expected_outputs = []
            for sample_name, files in sample_files.items():
                bam_file = files.get("bam")
                if not bam_file: continue # Skip if BAM is missing

                # Define output files relative to fragsize_output_dir
                hist_out = f"{sample_name}.fragmentSize.txt"
                plot_out = f"{sample_name}.fragmentSize.png" # bamPEFragmentSize can also plot
                expected_outputs.append(fragsize_output_dir / hist_out)

                # Use absolute path for input BAM
                # Output paths are relative because ParaFly will run with cwd=fragsize_output_dir
                cmd = (
                    f"{dt_frag_exe} --bamfiles {bam_file} "
                    f"--histogram {hist_out} "
                    f"--plotFile {plot_out} " # Generate plot directly
                    f"--numberOfProcessors {args.cores} " # Use main cores arg here
                    f"--maxFragmentLength 1000 " # Common setting for ATAC
                    f"--samplesLabel '{sample_name}'" # Ensure label is quoted if it contains spaces
                )
                commands.append(cmd)

            if commands:
                cmd_file = qc_log_dir / f"parafly_commands_{step_key}.txt"
                log_file = qc_log_dir / f"parafly_failed_{step_key}.log"
                frag_success = run_parafly_step(step_key, commands, cmd_file, log_file, args.parafly_jobs, fragsize_output_dir, args.no_run)
                if not frag_success: all_qc_successful = False

                # Optional: Check if output files were created
                if not args.no_run and frag_success:
                    missing_count = 0
                    for outfile in expected_outputs:
                        if not outfile.is_file():
                            missing_count += 1
                            print(f"[警告/WARN] 预期输出文件未生成: {outfile}")
                    if missing_count > 0:
                         print(f"[警告/WARN] {missing_count} 个片段大小文件未生成，请检查 ParaFly 日志。")
                         all_qc_successful = False
                    else:
                         print(f"[信息/INFO] 片段大小数据和图表已生成在: {fragsize_output_dir}")
                         print(f"       你可以使用这些 .txt 文件进行自定义绘图。")

            else:
                 print("[信息/INFO] 没有为片段大小计算生成命令。")

    # ==================== 2. TSS Enrichment Plot ====================
    if not args.skip_tss:
        print("\n--- 开始 TSS 富集分析 (Starting TSS Enrichment Analysis) ---")
        step_key = "tss_enrichment"
        dt_compute_exe = "computeMatrix"
        dt_heatmap_exe = "plotHeatmap"
        dt_profile_exe = "plotProfile"
        tss_output_dir = output_dir / "TSS_Enrichment"

        # Check tools
        if not find_executable(dt_compute_exe) or not find_executable(dt_heatmap_exe) or not find_executable(dt_profile_exe):
            missing_tools = [tool for tool in [dt_compute_exe, dt_heatmap_exe, dt_profile_exe] if not find_executable(tool)]
            print(f"[警告/WARN] 跳过 TSS 富集分析，因为以下命令未找到: {', '.join(missing_tools)}")
            all_qc_successful = False
        # Check if we have BigWig files to plot
        elif not any('bw' in files for files in sample_files.values()):
             print(f"[警告/WARN] 跳过 TSS 富集分析，因为未找到 BigWig 输入文件。")
             # Don't mark as failure, just can't perform this QC
        else:
            # Collect valid BigWig files and labels
            valid_bw_files = []
            valid_samples = []
            for sample, files in sample_files.items():
                 if 'bw' in files:
                      valid_bw_files.append(files['bw']) # Absolute path
                      valid_samples.append(sample)

            if not valid_bw_files:
                 print(f"[警告/WARN] 跳过 TSS 富集分析，因为未找到有效的 BigWig 文件。")
            else:
                bw_files_str = " ".join(valid_bw_files)
                sample_labels = " ".join(valid_samples)

                matrix_tss_out = tss_output_dir / "matrix_tss.gz"
                heatmap_tss_out = tss_output_dir / "heatmap_tss.png"
                profile_tss_out = tss_output_dir / "profile_tss.png"

                # --- Compute Matrix ---
                cmd_compute = (
                    f"{dt_compute_exe} scale-regions -S {bw_files_str} -R {tss_bed_path_abs} "
                    f"--beforeRegionStartLength 3000 --regionBodyLength 10 --afterRegionStartLength 3000 "
                    f"--skipZeros --missingDataAsZero -o {matrix_tss_out} -p {args.cores} "
                    f"--samplesLabel {sample_labels}"
                )
                compute_success, _ = run_single_command(f"{step_key}_computeMatrix", cmd_compute, tss_output_dir, args.no_run, check_error=False, executable_name=dt_compute_exe)

                if compute_success or args.no_run:
                    if args.no_run or matrix_tss_out.is_file():
                        # --- Plot Heatmap ---
                        cmd_heatmap = (
                            f"{dt_heatmap_exe} -m {matrix_tss_out} -out {heatmap_tss_out} "
                            f"--colorMap Blues --sortRegions keep --whatToShow 'heatmap and colorbar' "
                            f"--heatmapHeight 15 --heatmapWidth 5 --dpi 300"
                        )
                        heatmap_success, _ = run_single_command(f"{step_key}_plotHeatmap", cmd_heatmap, tss_output_dir, args.no_run, check_error=False, executable_name=dt_heatmap_exe)
                        if not heatmap_success: all_qc_successful = False

                        # --- Plot Profile ---
                        cmd_profile = (
                            f"{dt_profile_exe} -m {matrix_tss_out} -out {profile_tss_out} "
                            f"--plotTitle 'ATAC-seq Signal Profile around TSS' --perGroup "
                            f"--plotHeight 8 --plotWidth 10 --dpi 300"
                        )
                        profile_success, _ = run_single_command(f"{step_key}_plotProfile", cmd_profile, tss_output_dir, args.no_run, check_error=False, executable_name=dt_profile_exe)
                        if not profile_success: all_qc_successful = False

                        if heatmap_success and profile_success:
                             print(f"[信息/INFO] TSS 富集图已生成在: {tss_output_dir}")

                    else:
                        print(f"[警告/WARN] 跳过 TSS 绘图：矩阵文件 '{matrix_tss_out}' 未找到或 computeMatrix 失败。")
                        all_qc_successful = False
                else:
                    print(f"[警告/WARN] 跳过 TSS 绘图，因为 computeMatrix 步骤失败。")
                    all_qc_successful = False

    # ==================== 3. FRiP Score Calculation ====================
    if not args.skip_frip:
        print("\n--- 开始计算 FRiP 分数 (Starting FRiP Score Calculation) ---")
        step_key = "frip_score"
        frip_output_dir = output_dir / "FRiP_Scores"
        samtools_exe = "samtools"
        bedtools_exe = "bedtools"

        if not find_executable(samtools_exe) or not find_executable(bedtools_exe):
            missing_tools = [tool for tool in [samtools_exe, bedtools_exe] if not find_executable(tool)]
            print(f"[警告/WARN] 跳过 FRiP 计算，因为以下命令未找到: {', '.join(missing_tools)}")
            all_qc_successful = False
        elif not any('peak' in files for files in sample_files.values()):
             print(f"[警告/WARN] 跳过 FRiP 计算，因为未找到 Peak 输入文件。")
             # Don't mark as failure, just can't perform this QC
        else:
            os.makedirs(frip_output_dir, exist_ok=True)
            frip_results = []
            calculation_successful = True

            for sample_name, files in sample_files.items():
                print(f"  计算样本 (Calculating for sample): {sample_name}")
                bam_file = files.get("bam")
                peak_file = files.get("peak")

                if not bam_file or not peak_file:
                    print(f"[警告/WARN] 样本 {sample_name} 缺少 BAM 或 Peak 文件，跳过 FRiP 计算。")
                    frip_results.append({"Sample": sample_name, "Total Reads": "N/A", "Reads in Peaks": "N/A", "FRiP Score": "N/A", "Error": "Missing Input"})
                    calculation_successful = False
                    continue

                # 1. Get total mapped reads (primary alignments only, excluding supplementary/secondary/unmapped)
                # Using -F 2308 = exclude unmapped, secondary, supplementary, QC fail
                total_reads_cmd = f"{samtools_exe} view -c -F 2308 {bam_file}"
                success_total, total_reads_str = run_single_command(f"{step_key}_total_{sample_name}", total_reads_cmd, frip_output_dir, args.no_run, check_error=False, executable_name=samtools_exe, capture_stdout=True)

                # 2. Get reads overlapping peaks
                # bedtools intersect -u counts each read only once if it overlaps any peak
                # Use absolute paths for safety
                # Pipe to samtools view to count
                peak_reads_cmd = f"{bedtools_exe} intersect -u -abam {bam_file} -b {peak_file} | {samtools_exe} view -c -F 2308"
                success_peak, peak_reads_str = run_single_command(f"{step_key}_peak_{sample_name}", peak_reads_cmd, frip_output_dir, args.no_run, check_error=False, executable_name=bedtools_exe, capture_stdout=True) # Check bedtools

                # 3. Calculate FRiP
                frip_score = "N/A"
                error_msg = None
                total_reads = -1
                peak_reads = -1

                if not success_total or not success_peak:
                    error_msg = "Command Execution Failed"
                    calculation_successful = False
                elif total_reads_str is None or peak_reads_str is None:
                     error_msg = "Failed to capture command output"
                     calculation_successful = False
                else:
                    try:
                        total_reads = int(total_reads_str)
                        peak_reads = int(peak_reads_str)
                        if total_reads > 0:
                            frip_score = f"{(peak_reads / total_reads):.4f}" # Format to 4 decimal places
                        elif total_reads == 0:
                             frip_score = 0.0
                             print(f"[警告/WARN] 样本 {sample_name} 的总读数为 0。")
                        else:
                             error_msg = "Invalid Read Counts"
                             calculation_successful = False

                    except ValueError:
                        error_msg = "Non-integer Read Counts"
                        calculation_successful = False
                    except Exception as e:
                         error_msg = f"Calculation Error: {e}"
                         calculation_successful = False

                frip_results.append({
                    "Sample": sample_name,
                    "Total Reads": total_reads if total_reads >= 0 else "Error",
                    "Reads in Peaks": peak_reads if peak_reads >= 0 else "Error",
                    "FRiP Score": frip_score,
                    "Error": error_msg if error_msg else ""
                })
                print(f"    结果 (Result): Total={total_reads_str}, Peak={peak_reads_str}, FRiP={frip_score}{' (' + error_msg + ')' if error_msg else ''}")

            # Save FRiP results to a file
            if frip_results:
                frip_df = pd.DataFrame(frip_results)
                frip_summary_file = frip_output_dir / "frip_scores_summary.tsv"
                try:
                    frip_df.to_csv(frip_summary_file, sep="\t", index=False)
                    print(f"\n[信息/INFO] FRiP 分数摘要已保存到: {frip_summary_file}")
                except Exception as e:
                    print(f"[错误/ERROR] 无法保存 FRiP 摘要文件: {e}")
                    all_qc_successful = False

            if not calculation_successful:
                 all_qc_successful = False


    # --- Final Summary ---
    print("\n-------------------------------------------------------------")
    if all_qc_successful:
        print("ATAC-seq QC 脚本执行完毕。")
        print("(ATAC-seq QC Script finished execution.)")
    else:
        print("ATAC-seq QC 脚本在执行过程中遇到警告或错误。请检查上面的日志。")
        print("(ATAC-seq QC Script encountered warnings or errors during execution. Please review logs above.)")
    print("-------------------------------------------------------------")
    print(f"QC 结果保存在 (QC results saved in): {output_dir}")
    if (output_dir / "Fragment_Sizes").exists():
        print(f"  - 片段大小数据/图 (Fragment size data/plots): {output_dir / 'Fragment_Sizes'}")
    if (output_dir / "TSS_Enrichment").exists():
        print(f"  - TSS 富集图 (TSS enrichment plots): {output_dir / 'TSS_Enrichment'}")
    if (output_dir / "FRiP_Scores").exists():
        print(f"  - FRiP 分数摘要 (FRiP score summary): {output_dir / 'FRiP_Scores' / 'frip_scores_summary.tsv'}")
    print(f"  - 脚本日志 (Script logs): {qc_log_dir}")
    print("-------------------------------------------------------------")


if __name__ == "__main__":
    # Check for pandas dependency
    try:
        import pandas as pd
    except ImportError:
        print("[错误/ERROR] 此脚本需要 'pandas' 库来生成 FRiP 摘要。")
        print("请安装它，例如使用: conda install pandas  或  pip install pandas")
        sys.exit(1)
    main()

