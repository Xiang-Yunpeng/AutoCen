import pandas as pd
import os
import shlex
import subprocess
import sys

# ================= Helper Functions =================
def _get_cdhit_word_length(cdhit_id):
    """Dynamically determine CD-HIT word length (-n) based on identity threshold (-c)."""
    if cdhit_id >= 0.95: return 10
    elif cdhit_id >= 0.90: return 8
    elif cdhit_id >= 0.88: return 7
    elif cdhit_id >= 0.85: return 6
    elif cdhit_id >= 0.80: return 5
    else: return 4

def _parse_clusters(clstr_file):
    """Parse a CD-HIT .clstr file. Returns (sorted_clusters, cluster_reps)."""
    cluster_sizes = {}
    cluster_reps = {}
    
    with open(clstr_file) as f:
        for line in f:
            if line.startswith(">Cluster"):
                curr_id = line.strip().split()[1]
                cluster_sizes[curr_id] = 0
            else:
                cluster_sizes[curr_id] += 1
                if "*" in line:
                    rep_name = line.split('>')[1].split('...')[0]
                    cluster_reps[curr_id] = rep_name
                    
    sorted_clusters = sorted(cluster_sizes.items(), key=lambda x: x[1], reverse=True)
    return sorted_clusters, cluster_reps

def _extract_representatives(clustered_fasta, sorted_clusters, cluster_reps, out_dir, top_monomers, prefix, logger, min_cluster_size=5):
    """Extract representative sequences for top N families from CD-HIT output.
    
    The prefix controls output naming: 
      - Global/Standard mode uses "" -> Family_1_size100.fasta
      - Refined mode uses chrom name -> chr1_Family_1_size100.fasta
    
    Args:
        min_cluster_size: minimum number of members in a cluster to be output.
                          Clusters smaller than this are skipped (default: 5).
    """
    output_fastas = []
    
    # Build name prefix: "chr1_Family" for refined, "Family" for standard
    name_prefix = f"{prefix}_Family" if prefix else "Family"
    log_prefix = prefix if prefix else "Global"
    
    with open(clustered_fasta, 'r') as f_in:
        fasta_lines = f_in.readlines()

    rank = 0
    for c_id, size in sorted_clusters:
        if size < min_cluster_size:
            # Sorted descending, so all remaining clusters are also too small
            break
        
        rank += 1
        if rank > top_monomers:
            break
        
        rep_seq_name = cluster_reps[c_id]
        logger.info(f"  {log_prefix} Rank {rank}: Cluster {c_id}, Size: {size}, Rep: {rep_seq_name}")
        
        fam_fasta = os.path.join(out_dir, f"{name_prefix}_{rank}_size{size}.fasta")
        output_fastas.append(fam_fasta)
        
        with open(fam_fasta, 'w') as f_out:
            write = False
            for line in fasta_lines:
                if line.startswith(">"):
                    if line.strip().lstrip(">") == rep_seq_name:
                        write = True
                        f_out.write(f">{name_prefix}_{rank}_Cluster_{c_id}_Size_{size}\n")
                    else:
                        write = False
                elif write:
                    f_out.write(line)
                    
    return output_fastas

# ================= Original Global Clustering =================
def extract_and_cluster(trash_repeats_csv, candidate_bed, out_dir, cdhit_id, cdhit_aS, cdhit_aL, threads, top_monomers, logger, min_monomer_len=50, min_cluster_size=5):
    logger.info("--- Step 3: Monomer Extraction and Clustering ---")
    
    # 1. Load candidate regions
    logger.info(f"Loading candidate regions from: {candidate_bed}")
    candidate_regions = {}
    try:
        with open(candidate_bed, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                seqID, start, end = parts[0], int(parts[1]), int(parts[2])
                if seqID not in candidate_regions:
                    candidate_regions[seqID] = []
                candidate_regions[seqID].append((start, end))
    except FileNotFoundError:
        logger.error(f"Candidate BED file not found: {candidate_bed}")
        sys.exit(1)

    logger.info(f"Loaded {sum(len(v) for v in candidate_regions.values())} candidate regions.")

    # 2. Extract monomers
    target_fasta = os.path.join(out_dir, "candidate_monomers.fasta")
    extracted_count = 0
    skipped_short = 0
    chunk_size = 100000
    
    logger.info(f"Scanning TRASH repeats file: {trash_repeats_csv}")
    logger.info(f"Minimum monomer length filter: {min_monomer_len} bp")
    try:
        with open(target_fasta, 'w') as f_out:
            for chunk in pd.read_csv(trash_repeats_csv, chunksize=chunk_size, low_memory=False):
                for monomer in chunk.itertuples(index=False):
                    if monomer.seqID in candidate_regions:
                        for array_start, array_end in candidate_regions[monomer.seqID]:
                            if monomer.start >= array_start and monomer.end <= array_end:
                                seq = str(monomer.sequence)
                                if len(seq) < min_monomer_len:
                                    skipped_short += 1
                                    break
                                extracted_count += 1
                                fasta_header = f">{monomer.seqID}_{monomer.start}_{monomer.end}"
                                f_out.write(f"{fasta_header}\n{seq}\n")
                                break 
    except Exception as e:
        logger.error(f"Error extracting monomers: {e}")
        sys.exit(1)
        
    logger.info(f"Extraction complete! Found {extracted_count} candidate monomers "
                f"({skipped_short} skipped below {min_monomer_len} bp).")
    
    if extracted_count == 0:
        logger.error("No monomers extracted. Check coordinate matching.")
        sys.exit(1)

    # 3. CD-HIT Clustering
    clustered_fasta = os.path.join(out_dir, "clustered_monomers.fasta")
    clstr_file = clustered_fasta + ".clstr"
    cdhit_n = _get_cdhit_word_length(cdhit_id)
        
    logger.info(f"Starting CD-HIT clustering with identity (-c) {cdhit_id} and word length (-n) {cdhit_n}...")
    cmd = (f"cd-hit-est -i {shlex.quote(target_fasta)} -o {shlex.quote(clustered_fasta)} -c {cdhit_id} "
           f"-aS {cdhit_aS} -aL {cdhit_aL} -n {cdhit_n} -d 0 -M 0 -T {threads}")
    subprocess.run(cmd, shell=True, check=True, stdout=subprocess.DEVNULL)
    logger.info("CD-HIT clustering completed.")

    # 4. Parse clusters and extract representatives (reusing shared helpers)
    sorted_clusters, cluster_reps = _parse_clusters(clstr_file)
    logger.info(f"Identified top {min(top_monomers, len(sorted_clusters))} monomer families.")
    
    output_fastas = _extract_representatives(
        clustered_fasta, sorted_clusters, cluster_reps,
        out_dir, top_monomers, "", logger, min_cluster_size=min_cluster_size
    )
                    
    return output_fastas

# ================= Refined Per-Chromosome Clustering =================
def extract_and_cluster_per_chrom(trash_repeats_csv, candidate_bed, out_dir, cdhit_id, cdhit_aS, cdhit_aL, threads, top_monomers, logger, min_monomer=10, min_monomer_len=50, min_cluster_size=5):
    """
    Refined mode: cluster monomers independently for each chromosome.
    
    Key difference from global mode: reads CSV once, distributes monomers to 
    per-chromosome FASTA files, then runs CD-HIT separately per chromosome.
    
    Args:
        min_monomer: minimum number of monomers required per chromosome to run 
                     clustering (default: 10). Chromosomes below this threshold 
                     are skipped with a warning.
        min_monomer_len: minimum monomer sequence length in bp (default: 50).
                         Shorter monomers (e.g. microsatellites) are excluded.
        min_cluster_size: minimum CD-HIT cluster size to output a family (default: 5).
    
    Returns:
        per_chrom_results: dict { chrom_name: [fasta_path_1, fasta_path_2, ...] }
    """
    logger.info("--- Step 3: Per-Chromosome Monomer Extraction and Clustering (Refined Mode) ---")
    
    # 1. Load candidate regions grouped by chromosome
    logger.info(f"Loading candidate regions from: {candidate_bed}")
    candidate_regions = {}
    try:
        with open(candidate_bed, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                seqID, start, end = parts[0], int(parts[1]), int(parts[2])
                if seqID not in candidate_regions:
                    candidate_regions[seqID] = []
                candidate_regions[seqID].append((start, end))
    except FileNotFoundError:
        logger.error(f"Candidate BED file not found: {candidate_bed}")
        sys.exit(1)

    chromosomes = list(candidate_regions.keys())
    logger.info(f"Refined mode: {len(chromosomes)} chromosomes to process: {', '.join(chromosomes)}")

    # 2. Single-pass extraction: read CSV once, write per-chromosome FASTA files
    logger.info(f"Scanning TRASH repeats file (single-pass): {trash_repeats_csv}")
    logger.info(f"Minimum monomer length filter: {min_monomer_len} bp")
    
    chrom_dirs = {}
    chrom_fasta_paths = {}
    for chrom in chromosomes:
        chrom_dir = os.path.join(out_dir, chrom)
        os.makedirs(chrom_dir, exist_ok=True)
        chrom_dirs[chrom] = chrom_dir
        chrom_fasta_paths[chrom] = os.path.join(chrom_dir, f"{chrom}_monomers.fasta")
    
    extracted_counts = {c: 0 for c in chromosomes}
    skipped_short = 0
    
    try:
        chrom_file_handles = {c: open(chrom_fasta_paths[c], 'w') for c in chromosomes}
        
        for chunk in pd.read_csv(trash_repeats_csv, chunksize=100000, low_memory=False):
            for monomer in chunk.itertuples(index=False):
                if monomer.seqID in candidate_regions:
                    for array_start, array_end in candidate_regions[monomer.seqID]:
                        if monomer.start >= array_start and monomer.end <= array_end:
                            seq = str(monomer.sequence)
                            if len(seq) < min_monomer_len:
                                skipped_short += 1
                                break
                            chrom = monomer.seqID
                            extracted_counts[chrom] += 1
                            header = f">{monomer.seqID}_{monomer.start}_{monomer.end}"
                            chrom_file_handles[chrom].write(f"{header}\n{seq}\n")
                            break
                            
        for fh in chrom_file_handles.values():
            fh.close()
            
    except Exception as e:
        logger.error(f"Error during monomer extraction: {e}")
        sys.exit(1)

    for chrom in chromosomes:
        logger.info(f"  {chrom}: {extracted_counts[chrom]} monomers extracted")
    if skipped_short > 0:
        logger.info(f"  Total {skipped_short} monomers skipped (below {min_monomer_len} bp)")

    # 3. Per-chromosome CD-HIT clustering
    cdhit_n = _get_cdhit_word_length(cdhit_id)
    per_chrom_results = {}
    
    for chrom in chromosomes:
        if extracted_counts[chrom] == 0:
            logger.warning(f"  {chrom}: No monomers found, skipping.")
            continue
        
        if extracted_counts[chrom] < min_monomer:
            logger.warning(f"  {chrom}: Only {extracted_counts[chrom]} monomers extracted "
                           f"(below --min_monomer threshold of {min_monomer}), skipping.")
            continue
            
        logger.info(f"  Clustering {chrom} (identity={cdhit_id}, n={cdhit_n})...")
        
        chrom_dir = chrom_dirs[chrom]
        target_fasta = chrom_fasta_paths[chrom]
        clustered_fasta = os.path.join(chrom_dir, f"{chrom}_clustered.fasta")
        clstr_file = clustered_fasta + ".clstr"
        
        cmd = (f"cd-hit-est -i {shlex.quote(target_fasta)} -o {shlex.quote(clustered_fasta)} -c {cdhit_id} "
               f"-aS {cdhit_aS} -aL {cdhit_aL} -n {cdhit_n} -d 0 -M 0 -T {threads}")
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.DEVNULL)
        
        sorted_clusters, cluster_reps = _parse_clusters(clstr_file)
        
        output_fastas = _extract_representatives(
            clustered_fasta, sorted_clusters, cluster_reps,
            chrom_dir, top_monomers, chrom, logger, min_cluster_size=min_cluster_size
        )
        
        if output_fastas:
            per_chrom_results[chrom] = output_fastas
    
    if not per_chrom_results:
        logger.error("No monomers found on any chromosome in refined mode!")
        sys.exit(1)

    total_families = sum(len(v) for v in per_chrom_results.values())
    logger.info(f"Refined clustering complete: {total_families} families across {len(per_chrom_results)} chromosomes.")
    
    return per_chrom_results
