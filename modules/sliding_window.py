import os
import shlex
import subprocess

def run_sliding_window(genome_fa, monomer_fastas, out_dir, window_size, threads, blast_task, blast_identity, blast_evalue, blast_qcov, logger):
    logger.info("--- Step 4: Sliding Window Analysis ---")
    
    genome_sizes = os.path.join(out_dir, "genome.sizes")
    windows_bed = os.path.join(out_dir, f"windows_{window_size}.bed")
    
    q_genome = shlex.quote(genome_fa)
    q_sizes = shlex.quote(genome_sizes)
    q_windows = shlex.quote(windows_bed)
    
    # Generate genome windows
    subprocess.run(f"samtools faidx {q_genome}", shell=True, check=True)
    subprocess.run(f"cut -f1,2 {shlex.quote(genome_fa + '.fai')} > {q_sizes}", shell=True, check=True)
    subprocess.run(f"bedtools makewindows -g {q_sizes} -w {window_size} > {q_windows}", shell=True, check=True)
    
    # Build BLAST database once (instead of using -subject per query)
    db_prefix = os.path.join(out_dir, "genome_blastdb")
    q_db = shlex.quote(db_prefix)
    cmd_makedb = f"makeblastdb -in {q_genome} -dbtype nucl -out {q_db}"
    logger.info(f"Building BLAST database: {cmd_makedb}")
    subprocess.run(cmd_makedb, shell=True, check=True, stdout=subprocess.DEVNULL)
    
    density_tsvs = {}
    
    for fasta in monomer_fastas:
        fam_name = os.path.basename(fasta).replace(".fasta", "")
        logger.info(f"Running BLAST for {fam_name} with task '{blast_task}'...")
        
        blast_out = os.path.join(out_dir, f"{fam_name}.blast")
        monomer_bed = os.path.join(out_dir, f"{fam_name}_hits.bed")
        density_tsv = os.path.join(out_dir, f"{fam_name}_density.tsv")
        
        q_fasta = shlex.quote(fasta)
        q_blast_out = shlex.quote(blast_out)
        
        # BLAST against pre-built database (much faster than -subject for large genomes)
        cmd_blast = (
            f"blastn -task {blast_task} -query {q_fasta} -db {q_db} "
            f"-out {q_blast_out} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' "
            f"-perc_identity {blast_identity} -evalue {blast_evalue} -qcov_hsp_perc {blast_qcov} -num_threads {threads}"
        )
        subprocess.run(cmd_blast, shell=True, check=True)
        
        # Convert BLAST output to BED format
        # Since we used -qcov_hsp_perc, all hits here are already pre-filtered for length/coverage!
        with open(blast_out, 'r') as f_in, open(monomer_bed, 'w') as f_out:
            for line in f_in:
                cols = line.strip().split('\t')
                chrom = cols[1]
                s, e = int(cols[8]), int(cols[9])
                f_out.write(f"{chrom}\t{min(s, e)-1}\t{max(s, e)}\n")
                
        # Intersect with windows to calculate density
        q_bed = shlex.quote(monomer_bed)
        q_tsv = shlex.quote(density_tsv)
        subprocess.run(f"bedtools intersect -a {q_windows} -b {q_bed} -c > {q_tsv}", shell=True, check=True)
        density_tsvs[fam_name] = density_tsv
        
    return density_tsvs, genome_sizes
