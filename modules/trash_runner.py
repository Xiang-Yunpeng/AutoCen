import os
import shlex
import subprocess
import glob
import sys

def check_dependencies(genome_fa, trash_script, trash_dir, logger):
    logger.info("Checking dependencies...")
    tools = ["cd-hit-est", "blastn", "bedtools", "samtools"]
    for tool in tools:
        if subprocess.call(f"which {tool}", shell=True, stdout=subprocess.DEVNULL) != 0:
            logger.error(f"Tool not found in PATH: {tool}")
            sys.exit(1)

    if not os.path.exists(genome_fa):
        logger.error(f"Genome file not found: {genome_fa}")
        sys.exit(1)

    if trash_dir:
        if not os.path.exists(trash_dir):
            logger.error(f"Provided TRASH directory not found: {trash_dir}")
            sys.exit(1)
        logger.info("Running in 'Analyze Existing TRASH Results' mode.")
    else:
        if not trash_script or not os.path.exists(trash_script):
            logger.error("TRASH.R script not found. Please provide valid --trash_path.")
            sys.exit(1)
        if subprocess.call("which Rscript", shell=True, stdout=subprocess.DEVNULL) != 0:
            logger.error("Rscript not found in PATH.")
            sys.exit(1)
        logger.info("Running in 'Full Automation' mode.")
    logger.info("All dependencies checked successfully.")

def find_trash_files(search_dir, logger):
    arrays_csv = glob.glob(os.path.join(search_dir, "**", "*arrays.csv"), recursive=True)
    repeats_csv = glob.glob(os.path.join(search_dir, "**", "*repeats_with_seq.csv"), recursive=True)
    
    if not arrays_csv or not repeats_csv:
        logger.error("Could not find required TRASH output CSV files (arrays.csv and repeats_with_seq.csv).")
        sys.exit(1)
        
    logger.info(f"Found arrays.csv: {max(arrays_csv, key=os.path.getsize)}")
    logger.info(f"Found repeats_with_seq.csv: {max(repeats_csv, key=os.path.getsize)}")
    
    return max(arrays_csv, key=os.path.getsize), max(repeats_csv, key=os.path.getsize)

def run_trash(genome_fa, output_dir, trash_script, threads, logger):
    logger.info("--- Step 1: TRASH Analysis Phase ---")
    cmd = f"Rscript {shlex.quote(trash_script)} -f {shlex.quote(genome_fa)} -o {shlex.quote(output_dir)} -p {threads}"
    logger.info(f"Executing: {cmd}")
    
    try:
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        for line in process.stdout:
            line = line.strip()
            if line:
                logger.info(f"[TRASH] {line}")
        process.wait()
        if process.returncode != 0:
            logger.error("TRASH execution failed.")
            sys.exit(1)
    except Exception as e:
        logger.error(f"Failed to run TRASH: {e}")
        sys.exit(1)

    return find_trash_files(output_dir, logger)
