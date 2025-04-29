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

def find_executable(executable):
    """检查可执行文件是否存在于系统的 PATH 中。(Checks if an executable exists in the system's PATH.)"""
    return shutil.which(executable) is not None

def parse_args():
    """
    解析命令行参数 (Parses command-line arguments)
    """
    parser = argparse.ArgumentParser(
        description="ATAC-seq 可视化脚本：使用主流程的输出生成 DeepTools 图和 HOMER 注释。"
                    "(ATAC-seq Visualization Script: Generates DeepTools plots and HOMER annotations using outputs from the main pipeline.)"
    )
    # --- Input ---
    parser.add_argument(
        "-b", "--base_dir",
        required=True,
        help="主 ATAC-seq 流程的根输出目录路径。(Path to the base output directory of the main ATAC-seq pipeline.)"
    )
    parser.add_argument(
        "--tss_bed",
        required=True, # Make TSS BED required if running DeepTools part
        help="[DeepTools 需要/DeepTools requires] 提供一个包含TSS坐标的BED文件路径。(Path to a BED file containing TSS coordinates.)"
    )
    parser.add_argument(
        "--homer_genome",
        default=None,
        help="[HOMER 需要/HOMER requires] 提供HOMER兼容的基因组标识符 (例如: hg19, mm10)。(HOMER-compatible genome identifier (e.g., hg19, mm10).)"
    )
    # --- Output ---
    parser.add_argument(
        "-o", "--output_dir",
        required=True,
        help="保存可视化结果的输出目录路径。(Path to the output directory for saving visualization results.)"
    )
    # --- Tool Settings ---
    parser.add_argument(
        "--dt_cores",
        type=int,
        default=4,
        help="执行 DeepTools computeMatrix 时使用的线程数。(Number of threads for DeepTools computeMatrix) (默认/Default: 4)"
    )
    parser.add_argument(
        "--homer_cores",
        type=int,
        default=2,
        help="执行 HOMER annotatePeaks.pl 时每个作业使用的线程数。(Number of threads per job for HOMER annotatePeaks.pl) (默认/Default: 2)"
    )
    # --- Parallelization (for HOMER) ---
    parser.add_argument(
        "--parafly_jobs",
        type=int,
        default=4,
        help="如果需要并行运行多个 HOMER 作业，指定 ParaFly 的并行作业数。(Number of parallel HOMER jobs if using ParaFly) (默认/Default: 4)"
    )
    parser.add_argument(
        "--run_homer_sequentially",
        action="store_true",
        help="如果指定，则串行运行 HOMER 注释，而不是使用 ParaFly。(If specified, run HOMER annotations sequentially instead of using ParaFly.)"
    )
    # --- Execution Control ---
    parser.add_argument(
        "--no_run",
        action="store_true",
        help="如果指定，则只打印命令而不实际执行，用于调试。(If specified, only print commands without actual execution, for debugging.)"
    )
    return parser.parse_args()

def find_files(base_dir, pattern):
    """在指定的基础目录下查找匹配模式的文件。(Finds files matching a pattern within a base directory.)"""
    # Use glob to find files matching the pattern, including wildcards
    # sorted() ensures a consistent order
    found_files = sorted([str(p.resolve()) for p in Path(base_dir).glob(pattern)])
    return found_files

def run_single_command(step_name, cmd, output_dir, no_run=False, check_error=True, executable_name=None):
    """
    执行单个命令，如 DeepTools 绘图。
    (Executes a single final command like DeepTools plot)
    """
    print(f"\n--- 运行步骤 (Running Step): {step_name} ---")
    print(f"[命令/CMD] {cmd}")

    if not no_run:
        tool_to_check = executable_name or cmd.split()[0]
        if not find_executable(tool_to_check):
            print(f"[错误/ERROR] 命令 '{tool_to_check}' 未找到。(Command '{tool_to_check}' not found.)")
            print(f"    请确保它已安装并在您的 PATH 环境变量中。(Please ensure it is installed and in your PATH environment variable.)")
            if check_error:
                print(f"[错误/ERROR] 因缺少命令而无法执行步骤 {step_name}。终止。(Cannot execute step {step_name} due to missing command. Terminating.)")
                sys.exit(1)
            else:
                print(f"[警告/WARN] 因缺少命令而跳过步骤 {step_name}。(Skipping step {step_name} due to missing command.)")
                return False

    if not no_run:
        try:
            os.makedirs(output_dir, exist_ok=True)
            # Note: Running the command with shell=True means paths should be handled correctly by the shell,
            # but explicitly using absolute paths for inputs is safer.
            # We modify the calling code to pass absolute paths for input files.
            # The 'cwd' argument changes the working directory *for the command*, which caused the FileNotFoundError.
            # We keep cwd=output_dir so that output files are written there directly.
            process = subprocess.Popen(cmd, shell=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=output_dir)
            stdout, stderr = process.communicate()
            exit_code = process.returncode

            if stdout: print(f"[{step_name} STDOUT]\n{stdout.strip()}")
            if stderr: print(f"[{step_name} STDERR]\n{stderr.strip()}")

            if exit_code != 0:
                raise subprocess.CalledProcessError(exit_code, cmd, output=stdout, stderr=stderr)

            print(f"[信息/INFO] 步骤 {step_name} 成功执行。(Step {step_name} executed successfully.)")
            return True

        except subprocess.CalledProcessError as e:
            print(f"[错误/ERROR] 步骤 {step_name} 失败，退出代码 (failed with exit code) {e.returncode}")
            print(f"    命令 (Command): {e.cmd}")
            if check_error:
                print(f"[错误/ERROR] 因步骤 {step_name} 失败而退出。(Exiting due to failure in step {step_name}.)")
                sys.exit(1)
            else:
                print(f"[警告/WARN] 步骤 {step_name} 失败，但由于 check_error=False，执行将继续。(Step {step_name} failed, but execution continues as check_error=False.)")
                return False
        except Exception as e:
            print(f"[错误/ERROR] 步骤 {step_name} 执行期间发生意外错误 (An unexpected error occurred during step {step_name}): {e}")
            if check_error:
                sys.exit(1)
            else:
                return False
    else:
        print(f"[信息/INFO] --no_run 已指定。(specified.) 步骤 {step_name} 未执行。(Step {step_name} not executed.)")
        return True

def run_parafly_step(step_name, command_list, command_file_path, log_file_path, cpu, output_dir, no_run=False):
    """
    (复用自 atac_pipeline_visualize_v2) 写入命令文件并执行ParaFly。
    (Reused from atac_pipeline_visualize_v2) Writes command file and executes ParaFly.
    """
    if not command_list:
        print(f"[信息/INFO] 步骤 {step_name} 没有生成任何命令，跳过 ParaFly 执行。(No commands generated for step {step_name}. Skipping ParaFly execution.)")
        return True

    print(f"\n--- 准备并行步骤 (Preparing Parallel Step): {step_name} ---")
    print(f"    写入 {len(command_list)} 条命令到 (Writing {len(command_list)} commands to): {command_file_path}")
    try:
        os.makedirs(command_file_path.parent, exist_ok=True)
        with open(command_file_path, 'w') as f:
            f.write("\n".join(command_list) + "\n")
    except IOError as e:
        print(f"[错误/ERROR] 无法写入命令文件 (Failed to write command file) {command_file_path}: {e}")
        return False

    if not no_run and not find_executable("ParaFly"):
         print(f"[错误/ERROR] 命令 'ParaFly' 未找到。(Command 'ParaFly' not found.)")
         print(f"[错误/ERROR] 无法执行并行步骤 (Cannot execute parallel step) {step_name}. 终止。(Terminating.)")
         sys.exit(1)

    cmd_parafly = (
        f"ParaFly -c {command_file_path} -CPU {cpu} "
        f"-failed_cmds {log_file_path} -v"
    )

    print(f"--- 运行 ParaFly 步骤 (Running ParaFly for Step): {step_name} ---")
    print(f"[命令/CMD] {cmd_parafly}")
    if not no_run:
        try:
            os.makedirs(output_dir, exist_ok=True)
            # ParaFly executes commands listed in the command file.
            # Ensure paths within the command file are absolute or resolvable from the 'cwd'.
            # We modify the calling code (HOMER part) to use absolute paths in commands.
            process = subprocess.Popen(cmd_parafly, shell=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=output_dir)
            stdout, stderr = process.communicate()
            exit_code = process.returncode

            if stdout: print(f"[ParaFly STDOUT]\n{stdout.strip()}")
            if stderr: print(f"[ParaFly STDERR]\n{stderr.strip()}")

            failed_log_exists = log_file_path.exists() and log_file_path.stat().st_size > 0

            if exit_code != 0:
                if failed_log_exists:
                    print(f"[错误/ERROR] ParaFly 在步骤 {step_name} 中失败 (exit code {exit_code}) 并记录了失败的命令。(ParaFly failed in step {step_name} (exit code {exit_code}) AND logged failed commands.)")
                    print(f"    请检查日志文件 (Please check the log file): {log_file_path}")
                else:
                    print(f"[错误/ERROR] ParaFly 在步骤 {step_name} 中失败 (exit code {exit_code}) 但未记录失败的命令。(ParaFly failed in step {step_name} (exit code {exit_code}) but did not log failed commands.)")
                return False
            else:
                if failed_log_exists:
                    print(f"[警告/WARN] ParaFly 步骤 {step_name} 完成，但报告了失败的命令。(ParaFly step {step_name} completed BUT reported failed commands.)")
                    print(f"    请检查日志文件 (Please check the log file): {log_file_path}")
                else:
                    print(f"[信息/INFO] ParaFly 步骤 {step_name} 成功完成。(ParaFly for step {step_name} completed successfully.)")
                return True

        except Exception as e:
            print(f"[错误/ERROR] 运行 ParaFly 步骤 {step_name} 时发生意外错误 (An unexpected error occurred while running ParaFly for step {step_name}): {e}")
            return False
    else:
        print(f"[信息/INFO] --no_run 已指定。(specified.) ParaFly 命令步骤 {step_name} 未执行。(ParaFly command for step {step_name} not executed.)")
        return True

def main():
    args = parse_args()
    base_dir = Path(args.base_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    script_log_dir = output_dir / "Visualization_Scripts_Logs"

    # --- Create Output and Log Directories ---
    try:
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(script_log_dir, exist_ok=True)
        print(f"[信息/INFO] 输出目录 (Output directory): {output_dir}")
        print(f"[信息/INFO] 日志目录 (Log directory): {script_log_dir}")
    except OSError as e:
        print(f"[错误/ERROR] 无法创建输出或日志目录 (Failed to create output/log directory): {e}")
        sys.exit(1)

    # --- Define expected subdirectories from the main pipeline ---
    peaks_dir_pattern = "Run7_Peaks"
    bw_dir_pattern = "Run8_BigWig"

    # --- Find Input Files ---
    print("\n--- 查找输入文件 (Finding Input Files) ---")
    # find_files already returns absolute paths
    peak_files = find_files(base_dir, f"{peaks_dir_pattern}/*_peaks.narrowPeak")
    bw_files = find_files(base_dir, f"{bw_dir_pattern}/*.last.bw")

    if not peak_files:
        print(f"[警告/WARN] 在 '{base_dir / peaks_dir_pattern}' 中未找到 Peak 文件 (*_peaks.narrowPeak)。(No peak files found in '{base_dir / peaks_dir_pattern}'.) HOMER 注释将不可用。(HOMER annotation will be unavailable.)")
    else:
        print(f"  找到 {len(peak_files)} 个 Peak 文件。(Found {len(peak_files)} peak files.)")

    if not bw_files:
        print(f"[警告/WARN] 在 '{base_dir / bw_dir_pattern}' 中未找到 BigWig 文件 (*.last.bw)。(No BigWig files found in '{base_dir / bw_dir_pattern}'.) DeepTools 绘图将不可用。(DeepTools plotting will be unavailable.)")
    else:
        print(f"  找到 {len(bw_files)} 个 BigWig 文件。(Found {len(bw_files)} BigWig files.)")

    # --- Resolve TSS BED Path ---
    # Crucial Fix: Get the absolute path for the TSS BED file
    try:
        tss_bed_path_abs = Path(args.tss_bed).resolve(strict=True) # strict=True checks existence
        print(f"[信息/INFO] 使用 TSS BED 文件 (Using TSS BED file): {tss_bed_path_abs}")
    except FileNotFoundError:
        print(f"[错误/ERROR] 提供的 TSS BED 文件未找到 (Provided TSS BED file not found): {args.tss_bed}")
        sys.exit(1)
    except Exception as e:
         print(f"[错误/ERROR] 处理 TSS BED 文件路径时出错 (Error processing TSS BED file path) '{args.tss_bed}': {e}")
         sys.exit(1)


    # =========================================================================
    #                      VISUALIZATION STEPS START
    # =========================================================================
    visualization_successful = True

    # ==================== DeepTools Plots ====================
    # No need for 'if args.tss_bed:' check here as it's now a required argument
    print("\n--- 开始 DeepTools 绘图 (Starting DeepTools Plotting) ---")
    step_key = "deeptools_plots"
    dt_compute_exe = "computeMatrix"
    dt_heatmap_exe = "plotHeatmap"
    dt_profile_exe = "plotProfile"
    dt_output_dir = output_dir / "DeepTools_Plots"

    if not find_executable(dt_compute_exe) or not find_executable(dt_heatmap_exe) or not find_executable(dt_profile_exe):
        missing_tools = [tool for tool in [dt_compute_exe, dt_heatmap_exe, dt_profile_exe] if not find_executable(tool)]
        print(f"[警告/WARN] 跳过 DeepTools 绘图步骤，因为以下命令未在 PATH 中找到 (Skipping DeepTools plotting step because the following commands were not found in PATH): {', '.join(missing_tools)}")
        visualization_successful = False
    elif not bw_files:
         print(f"[警告/WARN] 跳过 DeepTools 绘图，因为未找到 BigWig 输入文件。(Skipping DeepTools plotting as no input BigWig files were found.)")
         visualization_successful = False
    else:
        # Files are already absolute paths from find_files
        bw_files_str = " ".join(bw_files)
        sample_labels = " ".join([Path(bw).name.replace("_Run8.last.bw", "") for bw in bw_files])

        matrix_tss_out = dt_output_dir / "matrix_tss.gz"
        heatmap_tss_out = dt_output_dir / "heatmap_tss.png"
        profile_tss_out = dt_output_dir / "profile_tss.png"

        # --- Compute Matrix ---
        # Use the absolute path for the TSS BED file
        cmd_compute = (
            f"{dt_compute_exe} scale-regions -S {bw_files_str} -R {tss_bed_path_abs} " # Use absolute path here
            f"--beforeRegionStartLength 3000 --regionBodyLength 10 --afterRegionStartLength 3000 "
            f"--skipZeros --missingDataAsZero -o {matrix_tss_out} -p {args.dt_cores} "
            f"--samplesLabel {sample_labels}"
        )
        # Pass dt_output_dir as the directory where the command should *write* outputs
        compute_success = run_single_command(f"{step_key}_computeMatrix", cmd_compute, dt_output_dir, args.no_run, check_error=False, executable_name=dt_compute_exe)

        if compute_success or args.no_run:
            if args.no_run or matrix_tss_out.is_file():
                # --- Plot Heatmap ---
                cmd_heatmap = (
                    f"{dt_heatmap_exe} -m {matrix_tss_out} -out {heatmap_tss_out} "
                    f"--colorMap Blues --sortRegions keep --whatToShow 'heatmap and colorbar' "
                    f"--heatmapHeight 15 --heatmapWidth 5 --dpi 300"
                )
                heatmap_success = run_single_command(f"{step_key}_plotHeatmap", cmd_heatmap, dt_output_dir, args.no_run, check_error=False, executable_name=dt_heatmap_exe)
                if not heatmap_success: visualization_successful = False

                # --- Plot Profile ---
                cmd_profile = (
                    f"{dt_profile_exe} -m {matrix_tss_out} -out {profile_tss_out} "
                    f"--plotTitle 'ATAC-seq Signal Profile around TSS' --perGroup "
                    f"--plotHeight 8 --plotWidth 10 --dpi 300"
                )
                profile_success = run_single_command(f"{step_key}_plotProfile", cmd_profile, dt_output_dir, args.no_run, check_error=False, executable_name=dt_profile_exe)
                if not profile_success: visualization_successful = False
            else:
                print(f"[警告/WARN] 跳过 DeepTools 绘图：矩阵文件 '{matrix_tss_out}' 未找到或 computeMatrix 失败。(Skipping DeepTools plotting: Matrix file '{matrix_tss_out}' not found or computeMatrix failed.)")
                visualization_successful = False
        else:
            print(f"[警告/WARN] 跳过 DeepTools 绘图，因为 computeMatrix 步骤失败。(Skipping DeepTools plotting because the computeMatrix step failed.)")
            visualization_successful = False


    # ==================== HOMER Annotation ====================
    if args.homer_genome:
        print("\n--- 开始 HOMER Peak 注释 (Starting HOMER Peak Annotation) ---")
        step_key = "homer_annotation"
        homer_executable = "annotatePeaks.pl"
        homer_output_dir = output_dir / "HOMER_Annotation"

        if not find_executable(homer_executable):
            print(f"[警告/WARN] 跳过 HOMER 注释步骤，因为 '{homer_executable}' 未在 PATH 中找到。(Skipping HOMER annotation step because '{homer_executable}' was not found in PATH.)")
            visualization_successful = False # Mark as less successful if optional step is skipped due to missing tool
        elif not peak_files:
            print(f"[警告/WARN] 跳过 HOMER 注释，因为未找到 Peak 输入文件。(Skipping HOMER annotation as no input Peak files were found.)")
            # Don't mark as unsuccessful if input is missing, just skip
        else:
            commands = []
            for peak_in_path_str in peak_files:
                # peak_files contains absolute paths already
                peak_in_path = Path(peak_in_path_str)
                sample_name = peak_in_path.name.replace("_Run7_peaks.narrowPeak", "")
                # Define output path relative to homer_output_dir
                anno_out_rel = f"{sample_name}.peak_annotation.txt" # Relative path for output file
                anno_out_abs = homer_output_dir / anno_out_rel # Absolute path for checking later

                # Create command using absolute path for input peak file
                # Output path can be relative as the command runs with cwd=homer_output_dir
                cmd = (
                    f"({homer_executable} {peak_in_path} {args.homer_genome} -gid -cpu {args.homer_cores} "
                    f"> {anno_out_rel} " # Output relative path
                    f"|| echo '[警告/WARN] HOMER annotation failed for {sample_name}' >&2)"
                )
                commands.append(cmd)

            if commands:
                if args.run_homer_sequentially or not find_executable("ParaFly"):
                    if args.run_homer_sequentially:
                         print("[信息/INFO] 串行运行 HOMER 注释。(Running HOMER annotations sequentially.)")
                    else:
                         print("[警告/WARN] 未找到 ParaFly，将串行运行 HOMER 注释。(ParaFly not found, running HOMER annotations sequentially.)")

                    homer_success_flag = True
                    for i, cmd in enumerate(commands):
                         sample_name = cmd.split("echo '[警告/WARN] HOMER annotation failed for ")[1].split("'")[0]
                         # Pass homer_output_dir as the directory where the command should run and write outputs
                         success = run_single_command(f"{step_key}_{sample_name}", cmd, homer_output_dir, args.no_run, check_error=False, executable_name=homer_executable)
                         if not success: homer_success_flag = False
                    if not homer_success_flag: visualization_successful = False

                else:
                    print(f"[信息/INFO] 使用 ParaFly 并行运行 {len(commands)} 个 HOMER 注释作业。(Running {len(commands)} HOMER annotation jobs in parallel using ParaFly.)")
                    cmd_file = script_log_dir / f"parafly_commands_{step_key}.txt"
                    log_file = script_log_dir / f"parafly_failed_{step_key}.log"
                    # Pass homer_output_dir as the directory where ParaFly should run the commands
                    homer_success_flag = run_parafly_step(step_key, commands, cmd_file, log_file, args.parafly_jobs, homer_output_dir, args.no_run)
                    if not homer_success_flag: visualization_successful = False
            else:
                 print("[信息/INFO] 没有为 HOMER 注释生成命令。(No commands generated for HOMER annotation.)")

    # --- Final Summary ---
    print("\n-------------------------------------------------------------")
    if visualization_successful:
        print("ATAC-seq 可视化脚本执行完毕。")
        print("(ATAC-seq Visualization Script finished execution.)")
    else:
        print("ATAC-seq 可视化脚本在执行过程中遇到警告或错误。请检查上面的日志。")
        print("(ATAC-seq Visualization Script encountered warnings or errors during execution. Please review logs above.)")
    print("-------------------------------------------------------------")
    print(f"结果保存在 (Results saved in): {output_dir}")
    # Check if directories were created, implying the step was attempted
    if (output_dir / "DeepTools_Plots").exists():
        print(f"  - DeepTools 图 (DeepTools plots): {output_dir / 'DeepTools_Plots'}")
    if (output_dir / "HOMER_Annotation").exists():
        print(f"  - HOMER 注释 (HOMER annotations): {output_dir / 'HOMER_Annotation'}")
    print(f"  - 脚本日志 (Script logs): {script_log_dir}")
    print("-------------------------------------------------------------")


if __name__ == "__main__":
    main()
