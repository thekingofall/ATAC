#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import argparse
import subprocess
from pathlib import Path
import sys

def parse_args():
    """
    解析命令行参数
    """
    parser = argparse.ArgumentParser(
        description="ATAC-seq分析流程 (分步并行版): 使用ParaFly分步并行处理样本，确保步骤间同步。"
    )
    # --- Input/Output ---
    parser.add_argument(
        "-i", "--input_dir",
        required=True,
        help="输入fastq文件所在文件夹路径。"
    )
    parser.add_argument(
        "-g", "--genome_index",
        default="/home/maolp/mao/Ref/AllnewstarRef/Homo/HG19/HG19BT/hg19",
        help="Bowtie2的参考基因组索引前缀(默认: /home/maolp/mao/Ref/AllnewstarRef/Homo/HG19/HG19BT/hg19)"
    )
    # --- Tool Settings (Per Job) ---
    parser.add_argument(
        "-f", "--fastqc_cores",
        type=int,
        default=4,
        help="执行trim_galore时FastQC使用的线程数 (每个并行作业)。(默认: 4)"
    )
    parser.add_argument(
        "-b", "--bowtie2_cores",
        type=int,
        default=6,
        help="执行bowtie2时使用的线程数 (每个并行作业)。(默认: 6)"
    )
    parser.add_argument(
        "-s", "--samtools_cores",
        type=int,
        default=3,
        help="执行samtools/sambamba时使用的线程数 (每个并行作业)。(默认: 3)"
    )
    parser.add_argument(
        "-p", "--peak_genome_size",
        default="hs",
        help="MACS2物种基因组大小参数(默认: hs)"
    )
    # --- Parallelization ---
    parser.add_argument(
        "--parafly_jobs",
        type=int,
        default=4,
        help="使用ParaFly并行运行的样本作业数。(默认: 4)"
    )
    # --- Visualization Options ---
    parser.add_argument(
        "--visualize",
        action="store_true",
        help="如果指定，则在流程结束后尝试生成额外的可视化图表 (需要DeepTools和可选的HOMER)。"
    )
    parser.add_argument(
        "--tss_bed",
        default=None,
        help="[可视化需要] 提供一个包含TSS坐标的BED文件路径，用于生成DeepTools的热图和谱图。"
    )
    parser.add_argument(
        "--homer_genome",
        default=None,
        help="[可视化需要] 提供HOMER兼容的基因组标识符 (例如: hg19, mm10)，用于annotatePeaks.pl。"
    )
    # --- Execution Control ---
    parser.add_argument(
        "--no_run",
        action="store_true",
        help="如果指定，则只生成命令文件和打印ParaFly命令，而不实际执行，用于调试。"
    )
    return parser.parse_args()

def find_paired_fastq_files(input_dir):
    """
    在 input_dir 中查找所有支持的fastq文件，自动配对R1和R2。
    """
    # (代码同前一个版本，保持不变)
    pattern = re.compile(
        r"^(?P<sample>.+?)"
        r"(?:[._](?:R|_)?(?P<read>[12]))"
        r"\.(?:fastq|fq)\.gz$"
    )
    paired_files = {}
    input_path = Path(input_dir)
    if not input_path.is_dir():
        print(f"[ERROR] 输入路径不是一个有效的文件夹: {input_dir}")
        sys.exit(1)

    for file_path in input_path.iterdir():
        if not file_path.is_file():
            continue
        match = pattern.match(file_path.name)
        if match:
            sample = match.group("sample")
            read = match.group("read")
            if sample not in paired_files:
                paired_files[sample] = ["", ""]

            current_path_str = str(file_path.resolve())
            if read == "1":
                if paired_files[sample][0]:
                    print(f"警告：样本 '{sample}' 发现多个 R1 文件。将使用: {file_path.name} (覆盖之前的: {os.path.basename(paired_files[sample][0])})")
                paired_files[sample][0] = current_path_str
            elif read == "2":
                if paired_files[sample][1]:
                    print(f"警告：样本 '{sample}' 发现多个 R2 文件。将使用: {file_path.name} (覆盖之前的: {os.path.basename(paired_files[sample][1])})")
                paired_files[sample][1] = current_path_str

    final_paired_files = {}
    for k, v in paired_files.items():
        if v[0] and v[1]:
            final_paired_files[k] = tuple(v)
        else:
            missing = "R1" if not v[0] else "R2"
            found = v[1] if not v[0] else v[0]
            print(f"警告：样本 '{k}' 未能成功配对 (缺少 {missing} 文件). 找到的文件: '{os.path.basename(found)}'. 将跳过此样本。")

    return final_paired_files

def run_parafly_step(step_name, command_list, command_file_path, log_file_path, cpu, no_run=False):
    """
    写入命令文件并执行ParaFly，等待完成。
    """
    if not command_list:
        print(f"[INFO] No commands generated for step {step_name}. Skipping ParaFly execution.")
        return

    print(f"\n--- Preparing Step: {step_name} ---")
    print(f"  Writing {len(command_list)} commands to: {command_file_path}")
    try:
        with open(command_file_path, 'w') as f:
            f.write("\n".join(command_list) + "\n")
    except IOError as e:
        print(f"[ERROR] Failed to write command file {command_file_path}: {e}")
        sys.exit(1)

    cmd_parafly = (
        f"ParaFly -c {command_file_path} -CPU {cpu} "
        f"-failed_cmds {log_file_path} -v"
    )

    print(f"--- Running ParaFly for Step: {step_name} ---")
    print(f"[CMD] {cmd_parafly}")
    if not no_run:
        try:
            # Run ParaFly and wait for it to complete
            process = subprocess.Popen(cmd_parafly, shell=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()
            exit_code = process.returncode

            if stdout: print(f"[STDOUT]\n{stdout.strip()}")
            if stderr: print(f"[STDERR]\n{stderr.strip()}")

            if exit_code != 0:
                 # Check the failed commands log even if ParaFly itself exits non-zero
                 if log_file_path.exists() and log_file_path.stat().st_size > 0:
                     print(f"[WARN] ParaFly exited non-zero ({exit_code}) AND reported failed commands. Check log: {log_file_path}")
                 else:
                     print(f"[ERROR] ParaFly command failed with exit code {exit_code} and no failed commands logged.")
                 print(f"    Command: {cmd_parafly}")
                 print(f"[ERROR] Exiting due to ParaFly failure in step {step_name}.")
                 sys.exit(1) # Exit if ParaFly fails for a step
            else:
                 # Check log even on success, ParaFly might exit 0 but still log failures
                 if log_file_path.exists() and log_file_path.stat().st_size > 0:
                     print(f"[WARN] ParaFly completed BUT reported failed commands for step {step_name}. Check log: {log_file_path}")
                     # Decide if this is critical - maybe allow continuing? For now, we warn.
                     # sys.exit(1) # Uncomment to make any logged failure critical
                 else:
                     print(f"[INFO] ParaFly for step {step_name} completed successfully.")

        except FileNotFoundError:
            print(f"[ERROR] Command not found: ParaFly. Is it installed and in PATH?")
            sys.exit(1)
        except Exception as e:
            print(f"[ERROR] An unexpected error occurred while running ParaFly for step {step_name}: {e}")
            sys.exit(1)
    else:
        print(f"[INFO] --no_run specified. ParaFly command for step {step_name} not executed.")


def run_single_command(step_name, cmd, no_run=False, check_error=True):
    """ Executes a single final command like MultiQC or DeepTools plot """
    print(f"\n--- Running Step: {step_name} ---")
    print(f"[CMD] {cmd}")
    if not no_run:
        try:
            process = subprocess.Popen(cmd, shell=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()
            exit_code = process.returncode

            if stdout: print(f"[STDOUT]\n{stdout.strip()}")
            if stderr: print(f"[STDERR]\n{stderr.strip()}")

            if exit_code != 0:
                raise subprocess.CalledProcessError(exit_code, cmd, output=stdout, stderr=stderr)
            print(f"[INFO] Step {step_name} executed successfully.")

        except subprocess.CalledProcessError as e:
            print(f"[ERROR] Step {step_name} failed with exit code {e.returncode}")
            print(f"    Command: {e.cmd}")
            if check_error:
                print(f"[ERROR] Exiting due to failure in step {step_name}.")
                sys.exit(1)
            else:
                print(f"[WARN] Step {step_name} failed, but execution continues as check_error=False.")
        except FileNotFoundError:
            print(f"[ERROR] Command not found for step {step_name}. Is the tool installed and in PATH? Command: {cmd.split()[0]}")
            if check_error:
                 print(f"[ERROR] Exiting due to missing command in step {step_name}.")
                 sys.exit(1)
            else:
                print(f"[WARN] Command not found for step {step_name}, but execution continues as check_error=False.")
        except Exception as e:
            print(f"[ERROR] An unexpected error occurred during step {step_name}: {e}")
            if check_error:
                 sys.exit(1)


def main():
    args = parse_args()

    # 1. Find paired files
    sample_dict = find_paired_fastq_files(args.input_dir)
    if not sample_dict:
        print("未找到任何可配对的FASTQ文件。流程终止。")
        return
    print(f"找到 {len(sample_dict)} 个配对样本。")

    # 2. Define and create directory structure
    base_out_dir = Path(".").resolve()
    dir_structure = {
        "run0_scripts": base_out_dir / "Run0_Scripts_Logs",
        "run1_trim": base_out_dir / "Run1_Trimmed_QC",
        "run2_align": base_out_dir / "Run2_Alignment",
        "run3_rmdup": base_out_dir / "Run3_RmDup_BAM",
        "run4_filter": base_out_dir / "Run4_Filtered_BAM",
        "run5_stats": base_out_dir / "Run5_Stats",
        "run6_bed": base_out_dir / "Run6_BED",
        "run7_peaks": base_out_dir / "Run7_Peaks",
        "run8_bw": base_out_dir / "Run8_BigWig",
        "run9_multiqc": base_out_dir / "Run9_MultiQC",
        "run10a_homer": base_out_dir / "Run10a_PeakAnnotation",
        "run10b_plots": base_out_dir / "Run10b_AggregatePlots",
    }
    for step_key, dir_path in dir_structure.items():
        if step_key.startswith("run10") and not args.visualize:
            continue
        os.makedirs(dir_path, exist_ok=True)

    # --- Helper to get paths ---
    def get_path(step_key, sample, filename_pattern):
        return dir_structure[step_key] / filename_pattern.format(sample=sample)

    # --- Define command file and log paths ---
    cmd_log_dir = dir_structure["run0_scripts"]
    def get_cmd_log_paths(step_key):
        cmd_file = cmd_log_dir / f"parafly_commands_{step_key}.txt"
        log_file = cmd_log_dir / f"parafly_failed_{step_key}.log"
        return cmd_file, log_file

    # --- Store expected outputs for later steps ---
    bw_files_expected = []
    peak_files_expected = {}


    # ==================== STEP 1: Trim Galore ====================
    step_key = "run1_trim"
    commands = []
    for sample, (fastq1, fastq2) in sample_dict.items():
        fastq1_abs = Path(fastq1).resolve()
        fastq2_abs = Path(fastq2).resolve()
        cmd = (
            f"trim_galore --paired --fastqc --cores {args.fastqc_cores} "
            f"--output_dir {dir_structure[step_key]} {fastq1_abs} {fastq2_abs}"
        )
        commands.append(cmd)
    cmd_file, log_file = get_cmd_log_paths(step_key)
    run_parafly_step(step_key, commands, cmd_file, log_file, args.parafly_jobs, args.no_run)

    # ==================== STEP 2: Alignment & Sort ====================
    step_key = "run2_align"
    commands = []
    for sample, (fastq1, fastq2) in sample_dict.items():
        # Determine expected trimmed filenames from Run1
        base1_name = Path(fastq1).name
        base2_name = Path(fastq2).name
        trimmed_base1, trimmed_base2 = None, None
        for suffix in [".fastq.gz", ".fq.gz"]:
            if base1_name.endswith(suffix): trimmed_base1 = base1_name[:-len(suffix)] + "_val_1.fq.gz"
            if base2_name.endswith(suffix): trimmed_base2 = base2_name[:-len(suffix)] + "_val_2.fq.gz"

        if not trimmed_base1 or not trimmed_base2:
             print(f"[WARN] 无法确定 {sample} 的 Trim Galore 输出名，跳过 {step_key}。")
             continue

        fastq1_trimmed = dir_structure["run1_trim"] / trimmed_base1
        fastq2_trimmed = dir_structure["run1_trim"] / trimmed_base2
        bam_out = get_path(step_key, sample, "{sample}_Run2.bam")

        # Check if input files exist before adding command (optional but good practice)
        if not args.no_run and (not fastq1_trimmed.is_file() or not fastq2_trimmed.is_file()):
            print(f"[WARN] Trimmed input for {sample} not found ({fastq1_trimmed}, {fastq2_trimmed}), skipping {step_key}.")
            continue

        cmd = (
            f"(bowtie2 -p {args.bowtie2_cores} -X 1000 --no-mixed --no-discordant --no-unal "
            f"-x {args.genome_index} -1 {fastq1_trimmed} -2 {fastq2_trimmed} | "
            f"samtools sort -O bam -@ {args.samtools_cores} -o {bam_out} -)"
        )
        commands.append(cmd)
    cmd_file, log_file = get_cmd_log_paths(step_key)
    run_parafly_step(step_key, commands, cmd_file, log_file, args.parafly_jobs, args.no_run)

    # ==================== STEP 3: Remove Duplicates ====================
    step_key = "run3_rmdup"
    commands = []
    for sample in sample_dict.keys():
        bam_in = get_path("run2_align", sample, "{sample}_Run2.bam")
        bam_out = get_path(step_key, sample, "{sample}_Run3.sambamba.rmdup.bam")
        if not args.no_run and not bam_in.is_file():
            print(f"[WARN] Input BAM for {sample} ({bam_in}) not found, skipping {step_key}.")
            continue
        cmd = (
            f"sambamba markdup -r --nthreads={args.samtools_cores} "
            f"{bam_in} {bam_out}"
        )
        commands.append(cmd)
    cmd_file, log_file = get_cmd_log_paths(step_key)
    run_parafly_step(step_key, commands, cmd_file, log_file, args.parafly_jobs, args.no_run)

    # ==================== STEP 4: Filter BAM ====================
    step_key = "run4_filter"
    commands = []
    for sample in sample_dict.keys():
        bam_in = get_path("run3_rmdup", sample, "{sample}_Run3.sambamba.rmdup.bam")
        bam_out = get_path(step_key, sample, "{sample}_Run4.last.bam")
        if not args.no_run and not bam_in.is_file():
             print(f"[WARN] Input BAM for {sample} ({bam_in}) not found, skipping {step_key}.")
             continue
        cmd = (
             f"(samtools view -h -f 2 -q 30 -@ {args.samtools_cores} {bam_in} | "
             f"grep -v 'chrM' | "
             f"samtools sort -O bam -@ {args.samtools_cores} -o {bam_out} -)"
        )
        commands.append(cmd)
    cmd_file, log_file = get_cmd_log_paths(step_key)
    run_parafly_step(step_key, commands, cmd_file, log_file, args.parafly_jobs, args.no_run)

    # ==================== STEP 4b: Index Final BAM ====================
    step_key = "run4b_index" # Treat indexing as a separate parallel step
    commands = []
    for sample in sample_dict.keys():
        bam_to_index = get_path("run4_filter", sample, "{sample}_Run4.last.bam")
        if not args.no_run and not bam_to_index.is_file():
            print(f"[WARN] BAM file for {sample} ({bam_to_index}) not found, skipping indexing.")
            continue
        cmd = f"samtools index -@ {args.samtools_cores} {bam_to_index}"
        commands.append(cmd)
    cmd_file, log_file = get_cmd_log_paths(step_key)
    run_parafly_step(step_key, commands, cmd_file, log_file, args.parafly_jobs, args.no_run)

    # ==================== STEP 5a: Index Run3 BAM (for Stats) ====================
    step_key = "run5a_index_run3"
    commands = []
    for sample in sample_dict.keys():
        bam_to_index = get_path("run3_rmdup", sample, "{sample}_Run3.sambamba.rmdup.bam")
        if not args.no_run and not bam_to_index.is_file():
             print(f"[WARN] Run3 BAM file for {sample} ({bam_to_index}) not found, skipping indexing for stats.")
             continue
        # Allow index command to fail if file already exists using || true
        cmd = f"(samtools index -@ {args.samtools_cores} {bam_to_index} || true)"
        commands.append(cmd)
    cmd_file, log_file = get_cmd_log_paths(step_key)
    run_parafly_step(step_key, commands, cmd_file, log_file, args.parafly_jobs, args.no_run)


    # ==================== STEP 5b: Calculate Stats ====================
    step_key = "run5b_stats"
    commands = []
    for sample in sample_dict.keys():
        bam_in = get_path("run3_rmdup", sample, "{sample}_Run3.sambamba.rmdup.bam")
        stat_out = get_path("run5_stats", sample, "{sample}_Run5.rmdup.stat") # Output to Run5 dir
        if not args.no_run and not bam_in.is_file():
            print(f"[WARN] Run3 BAM for {sample} ({bam_in}) not found, skipping {step_key}.")
            continue
        # Requires Run3 BAM index from step 5a
        cmd = f"samtools flagstat -@ {args.samtools_cores} {bam_in} > {stat_out}"
        commands.append(cmd)
    cmd_file, log_file = get_cmd_log_paths(step_key)
    run_parafly_step(step_key, commands, cmd_file, log_file, args.parafly_jobs, args.no_run)


    # ==================== STEP 6: BAM to BED ====================
    step_key = "run6_bed"
    commands = []
    for sample in sample_dict.keys():
        bam_in = get_path("run4_filter", sample, "{sample}_Run4.last.bam")
        bed_out = get_path(step_key, sample, "{sample}_Run6.last.bed")
        if not args.no_run and not bam_in.is_file():
             print(f"[WARN] Run4 BAM for {sample} ({bam_in}) not found, skipping {step_key}.")
             continue
        cmd = f"bedtools bamtobed -i {bam_in} > {bed_out}"
        commands.append(cmd)
    cmd_file, log_file = get_cmd_log_paths(step_key)
    run_parafly_step(step_key, commands, cmd_file, log_file, args.parafly_jobs, args.no_run)

    # ==================== STEP 7: Peak Calling (MACS2) ====================
    step_key = "run7_peaks"
    commands = []
    for sample in sample_dict.keys():
        bed_in = get_path("run6_bed", sample, "{sample}_Run6.last.bed")
        macs2_out_prefix = f"{sample}_Run7"
        macs2_peak_file = dir_structure[step_key] / f"{macs2_out_prefix}_peaks.narrowPeak"
        peak_files_expected[sample] = macs2_peak_file # Store expected path

        if not args.no_run and not bed_in.is_file():
            print(f"[WARN] Input BED for {sample} ({bed_in}) not found, skipping {step_key}.")
            continue

        cmd = (
            f"macs2 callpeak -t {bed_in} -f BED -g {args.peak_genome_size} "
            f"--nomodel --shift -100 --extsize 200 "
            f"-n {macs2_out_prefix} --outdir {dir_structure[step_key]}"
        )
        commands.append(cmd)
    cmd_file, log_file = get_cmd_log_paths(step_key)
    run_parafly_step(step_key, commands, cmd_file, log_file, args.parafly_jobs, args.no_run)

    # ==================== STEP 8: BigWig Generation ====================
    step_key = "run8_bw"
    commands = []
    for sample in sample_dict.keys():
        bam_in = get_path("run4_filter", sample, "{sample}_Run4.last.bam")
        bw_out = get_path(step_key, sample, "{sample}_Run8.last.bw")
        bw_files_expected.append(str(bw_out)) # Store expected path

        if not args.no_run and not bam_in.is_file():
            print(f"[WARN] Run4 BAM for {sample} ({bam_in}) not found, skipping {step_key}.")
            continue
        # Requires Run4 BAM index from step 4b
        cmd = (
            f"bamCoverage --normalizeUsing CPM --numberOfProcessors {args.samtools_cores} "
            f"-b {bam_in} -o {bw_out}"
        )
        commands.append(cmd)
    cmd_file, log_file = get_cmd_log_paths(step_key)
    run_parafly_step(step_key, commands, cmd_file, log_file, args.parafly_jobs, args.no_run)

    # ==================== STEP 10a: HOMER Annotation (Optional) ====================
    if args.visualize and args.homer_genome:
        step_key = "run10a_homer"
        commands = []
        for sample in sample_dict.keys():
            # Use the expected peak file path generated in step 7
            peak_in = peak_files_expected.get(sample)
            anno_out = get_path(step_key, sample, "{sample}_Run10a.peak_annotation.txt")

            if not peak_in:
                 print(f"[WARN] Expected peak file path for {sample} not found, skipping {step_key}.")
                 continue
            if not args.no_run and not peak_in.is_file():
                 print(f"[WARN] Input Peak file for {sample} ({peak_in}) not found, skipping {step_key}.")
                 continue

            # Allow homer to fail without stopping the script
            cmd = (
                f"(annotatePeaks.pl {peak_in} {args.homer_genome} -gid -cpu {args.samtools_cores} "
                f"> {anno_out} || echo '[WARN] HOMER annotation failed for {sample}')"
            )
            commands.append(cmd)
        cmd_file, log_file = get_cmd_log_paths(step_key)
        run_parafly_step(step_key, commands, cmd_file, log_file, args.parafly_jobs, args.no_run)

    # --- Final Single-Threaded Steps ---

    # ==================== STEP 9: MultiQC ====================
    step_key = "run9_multiqc"
    cmd_multiqc = f"multiqc {base_out_dir} -o {dir_structure[step_key]}"
    run_single_command(step_key, cmd_multiqc, args.no_run, check_error=False) # Don't exit if MultiQC fails

    # ==================== STEP 10b: DeepTools Plots (Optional) ====================
    if args.visualize and args.tss_bed:
        step_key = "run10b_plots"
        print(f"\n--- Preparing Step: {step_key} ---")
        if not Path(args.tss_bed).is_file():
            print(f"[WARN] Skipping DeepTools plots: Provided TSS BED file not found: {args.tss_bed}")
        else:
            # Check which expected BigWig files actually exist
            valid_bw_files = [bw for bw in bw_files_expected if Path(bw).exists()]
            if not valid_bw_files:
                 print(f"[WARN] Skipping DeepTools plots: No BigWig files found after Run8.")
            else:
                matrix_tss_out = dir_structure[step_key] / "Run10b_matrix_tss.gz"
                heatmap_tss_out = dir_structure[step_key] / "Run10b_heatmap_tss.png"
                profile_tss_out = dir_structure[step_key] / "Run10b_profile_tss.png"
                bw_files_str = " ".join(valid_bw_files)
                valid_samples = [Path(bw).name.replace("_Run8.last.bw", "") for bw in valid_bw_files]
                sample_labels = " ".join(valid_samples)

                # Compute Matrix (can use more cores here)
                total_cores_available = args.parafly_jobs * args.samtools_cores # Rough estimate
                compute_cores = max(1, total_cores_available // 2) # Use half, at least 1
                cmd_compute = (
                    f"computeMatrix scale-regions -S {bw_files_str} -R {args.tss_bed} "
                    f"--beforeRegionStartLength 3000 --regionBodyLength 10 --afterRegionStartLength 3000 "
                    f"--skipZeros --missingDataAsZero -o {matrix_tss_out} -p {compute_cores} "
                    f"--samplesLabel {sample_labels}"
                )
                run_single_command(f"{step_key}_computeMatrix", cmd_compute, args.no_run, check_error=False)

                if args.no_run or matrix_tss_out.is_file():
                    # Plot Heatmap
                    cmd_heatmap = (
                        f"plotHeatmap -m {matrix_tss_out} -out {heatmap_tss_out} "
                        f"--colorMap Blues --sortRegions keep --whatToShow 'heatmap and colorbar' "
                        f"--heatmapHeight 15 --heatmapWidth 5"
                    )
                    run_single_command(f"{step_key}_plotHeatmap", cmd_heatmap, args.no_run, check_error=False)
                    # Plot Profile
                    cmd_profile = (
                        f"plotProfile -m {matrix_tss_out} -out {profile_tss_out} "
                        f"--plotTitle 'ATAC-seq Signal Profile around TSS'"
                    )
                    run_single_command(f"{step_key}_plotProfile", cmd_profile, args.no_run, check_error=False)
                else:
                     print(f"[WARN] Skipping DeepTools plotting: Matrix file '{matrix_tss_out}' not found or computeMatrix failed.")

    elif args.visualize:
         print("\n[INFO] Skipping DeepTools aggregate plots: --tss_bed argument not provided.")


    # --- Final Summary ---
    print("\n--------------------")
    print("ATAC-seq Pipeline (Step-wise Parallel with ParaFly) Finished.")
    print("--------------------")
    print("Key output directories:")
    # (Print directory structure using updated keys)
    print(f"  - {dir_structure['run0_scripts']}: Generated command files and ParaFly logs")
    print(f"  - {dir_structure['run1_trim']}: Trimmed FASTQ files and FastQC reports")
    print(f"  - {dir_structure['run2_align']}: Initial aligned BAM files")
    print(f"  - {dir_structure['run3_rmdup']}: Deduplicated BAM files")
    print(f"  - {dir_structure['run4_filter']}: Final filtered BAM files and index")
    print(f"  - {dir_structure['run5_stats']}: Alignment statistics (flagstat)")
    print(f"  - {dir_structure['run6_bed']}: Final BED file (from filtered BAM)")
    print(f"  - {dir_structure['run7_peaks']}: MACS2 peak calling results")
    print(f"  - {dir_structure['run8_bw']}: Normalized signal tracks (BigWig)")
    print(f"  - {dir_structure['run9_multiqc']}: Aggregate QC report (MultiQC)")
    if args.visualize:
        if args.homer_genome: print(f"  - {dir_structure['run10a_homer']}: HOMER Peak annotation results")
        if args.tss_bed: print(f"  - {dir_structure['run10b_plots']}: DeepTools aggregate plots")
    print("\nNext Steps Suggestions:")
    print(f"  - Check ParaFly logs in '{dir_structure['run0_scripts']}'.")
    print(f"  - Inspect MultiQC report in '{dir_structure['run9_multiqc']}/multiqc_report.html'.")
    print(f"  - Load final BAM ({dir_structure['run4_filter']}/<sample>_Run4.last.bam), BigWig ({dir_structure['run8_bw']}/<sample>_Run8.last.bw), and Peaks ({dir_structure['run7_peaks']}/<sample>_Run7_peaks.narrowPeak) into IGV.")
    if args.visualize and args.homer_genome: print(f"  - Check HOMER peak annotation results in {dir_structure['run10a_homer']}/<sample>_Run10a.peak_annotation.txt")
    if args.visualize and args.tss_bed: print(f"  - Check DeepTools aggregate plots in {dir_structure['run10b_plots']}/")


if __name__ == "__main__":
    main()
