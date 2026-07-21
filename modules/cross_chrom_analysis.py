import os
import shlex
import subprocess
import sys
import numpy as np

from utils.chrom_utils import natural_sort_key

def shorten_family_label(label):
    """Shorten FASTA header to readable label.
    
    'chr1_Family_1_Cluster_0_Size_100' -> 'chr1_Fam1(100)'
    Splits on '_Family_' to safely handle chrom names containing underscores.
    """
    if '_Family_' in label:
        chrom, rest = label.split('_Family_', 1)
        parts = rest.split('_')
        fam_rank = parts[0]
        size = ""
        for i, p in enumerate(parts):
            if p == "Size" and i + 1 < len(parts):
                size = parts[i + 1]
        return f"{chrom}_Fam{fam_rank}({size})"
    return label

def _get_cluster_size(label):
    """Extract cluster size from a full family label for representative selection.
    
    'chr1_Family_1_Cluster_0_Size_100' -> 100
    Returns 0 if parsing fails.
    """
    parts = label.split('_')
    for i, p in enumerate(parts):
        if p == "Size" and i + 1 < len(parts):
            try:
                return int(parts[i + 1])
            except ValueError:
                return 0
    return 0

# ================= Family Merging =================
def merge_families(sim_matrix, family_labels, short_labels, per_chrom_results, 
                   out_dir, merge_threshold, logger):
    """
    Merge per-chromosome families into final unified families using the identity matrix.
    
    Algorithm: complete-linkage agglomerative clustering.
    Two clusters merge ONLY when the minimum pairwise identity between ALL members
    of both clusters >= threshold. This prevents the chaining effect where A~B~C~D
    leads to A and D being grouped despite low direct identity.
    
    At each iteration, the pair of clusters with the highest minimum pairwise identity
    is merged first (greedy best-first). This is O(n^3) but n is typically <100.
    
    The representative for each merged family is the member with the largest CD-HIT 
    cluster size (i.e., the most abundant monomer variant).
    
    Args:
        sim_matrix:       n x n identity matrix
        family_labels:    full FASTA header names
        short_labels:     shortened display labels
        per_chrom_results: dict { chrom: [fasta_path, ...] }
        out_dir:          output directory (cross_chrom_comparison/)
        merge_threshold:  identity cutoff for merging (default 80.0)
        logger:           logger instance
        
    Returns:
        merged_fastas:    list of final family FASTA paths
        merge_summary:    path to the merge summary TSV
    """
    n = len(family_labels)
    cross_dir = os.path.join(out_dir, "cross_chrom_comparison")
    final_dir = os.path.join(cross_dir, "final_families")
    os.makedirs(final_dir, exist_ok=True)
    
    logger.info(f"--- Merging Families (complete-linkage, threshold: {merge_threshold}%) ---")
    
    # 1. Complete-linkage agglomerative clustering
    #    Start with each family as its own cluster
    clusters = [[i] for i in range(n)]
    
    while True:
        best_merge = None
        best_min_id = -1.0
        
        for a in range(len(clusters)):
            for b in range(a + 1, len(clusters)):
                # Complete linkage: the score is the MINIMUM identity across all inter-cluster pairs
                min_id = min(sim_matrix[i][j] for i in clusters[a] for j in clusters[b])
                if min_id >= merge_threshold and min_id > best_min_id:
                    best_min_id = min_id
                    best_merge = (a, b)
        
        if best_merge is None:
            break
        
        # Merge cluster b into cluster a, remove b
        a, b = best_merge
        clusters[a].extend(clusters[b])
        clusters.pop(b)
        logger.info(f"  Merged two groups (min identity: {best_min_id:.1f}%), "
                     f"new group size: {len(clusters[a])}")
    
    logger.info(f"Complete-linkage clustering done: {n} per-chrom families -> {len(clusters)} merged families")
    
    # Convert to dict format compatible with downstream code
    groups = {idx: members for idx, members in enumerate(clusters)}
    
    # Sort groups by total monomer count (descending) for intuitive ranking
    def group_total_count(root):
        return sum(_get_cluster_size(family_labels[idx]) for idx in groups[root])
    
    sorted_roots = sorted(groups.keys(), key=group_total_count, reverse=True)
    
    # 2. Read all representative sequences into memory
    #    Build a lookup: family_label -> sequence
    seq_lookup = {}
    for chrom in per_chrom_results:
        for fasta in per_chrom_results[chrom]:
            current_header = None
            current_seq = []
            with open(fasta, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        if current_header is not None:
                            seq_lookup[current_header] = ''.join(current_seq)
                        current_header = line.strip().lstrip('>')
                        current_seq = []
                    else:
                        current_seq.append(line.strip())
                if current_header is not None:
                    seq_lookup[current_header] = ''.join(current_seq)
    
    # 3. For each group, pick the representative and write output
    merged_fastas = []
    summary_rows = []
    
    for final_rank, root in enumerate(sorted_roots, 1):
        members = groups[root]
        
        # Select representative: the member with the largest cluster size
        best_idx = max(members, key=lambda idx: _get_cluster_size(family_labels[idx]))
        rep_label = family_labels[best_idx]
        rep_short = short_labels[best_idx]
        rep_size = _get_cluster_size(rep_label)
        
        # Calculate total monomer count across all member chromosomes
        total_count = sum(_get_cluster_size(family_labels[idx]) for idx in members)
        
        # Collect member info
        member_chroms = []
        member_shorts = []
        for idx in sorted(members, key=lambda x: natural_sort_key(short_labels[x])):
            member_shorts.append(short_labels[idx])
            chrom = family_labels[idx].split('_Family_')[0] if '_Family_' in family_labels[idx] else "?"
            if chrom not in member_chroms:
                member_chroms.append(chrom)
        
        member_chroms_sorted = sorted(member_chroms, key=natural_sort_key)
        is_shared = len(member_chroms_sorted) > 1
        
        # Simplified naming: Fam1_n35001, Fam2_n5602, ... (details go into summary TSV)
        family_name = f"Fam{final_rank}_n{total_count}"
        
        # Write FASTA
        fasta_path = os.path.join(final_dir, f"{family_name}.fasta")
        seq = seq_lookup.get(rep_label, "")
        if seq:
            with open(fasta_path, 'w') as f:
                f.write(f">{family_name} rep={rep_short} members={len(members)}\n")
                # Write sequence in 80-char lines
                for i in range(0, len(seq), 80):
                    f.write(seq[i:i+80] + "\n")
            merged_fastas.append(fasta_path)
        else:
            logger.warning(f"  Could not find sequence for {rep_label}, skipping.")
            continue
        
        # Log
        shared_tag = "SHARED" if is_shared else "specific"
        logger.info(f"  Final Family {final_rank} [{shared_tag}]: rep={rep_short}, "
                     f"total_monomers={total_count}, members={len(members)}, "
                     f"chromosomes={','.join(member_chroms_sorted)}")
        
        summary_rows.append({
            'family_id': final_rank,
            'family_name': family_name,
            'representative': rep_short,
            'total_monomers': total_count,
            'member_count': len(members),
            'members': "; ".join(member_shorts),
            'chromosomes': ", ".join(member_chroms_sorted),
            'shared': "Yes" if is_shared else "No"
        })
    
    # 4. Write summary TSV
    merge_summary = os.path.join(cross_dir, "merged_families_summary.tsv")
    with open(merge_summary, 'w') as f:
        f.write("Family_ID\tFamily_Name\tRepresentative\tTotal_Monomers\tMembers\tMember_Count\tChromosomes\tShared\n")
        for row in summary_rows:
            f.write(f"{row['family_id']}\t{row['family_name']}\t{row['representative']}\t"
                    f"{row['total_monomers']}\t{row['members']}\t{row['member_count']}\t"
                    f"{row['chromosomes']}\t{row['shared']}\n")
    
    # 5. Log final statistics
    n_shared = sum(1 for r in summary_rows if r['shared'] == "Yes")
    n_specific = sum(1 for r in summary_rows if r['shared'] == "No")
    logger.info(f"Family merging complete: {len(summary_rows)} final families "
                f"({n_shared} shared across chromosomes, {n_specific} chromosome-specific)")
    logger.info(f"Final family FASTAs saved to: {final_dir}/")
    logger.info(f"Merge summary saved to: {merge_summary}")
    
    return merged_fastas, merge_summary


# ================= Main Cross-Chromosome Analysis =================
def cross_chromosome_analysis(per_chrom_results, out_dir, threads, logger, min_qcov=80.0):
    """
    Compare representative monomer sequences across chromosomes using all-vs-all BLAST,
    then merge same-family sequences into final unified families.
    
    Similarity metric: pure percent identity (pident) from the best HSP that passes
    a **bidirectional** coverage filter. Both query coverage and subject coverage must
    reach the threshold, ensuring symmetric treatment regardless of which sequence is
    query vs subject. This prevents length-biased asymmetry when comparing monomers 
    of different sizes.
    
    Args:
        per_chrom_results: dict { chrom: [fasta_path, ...] } from refined clustering
        out_dir: output directory (typically 3_monomer_clustering)
        threads: number of CPU threads
        logger: logger instance
        min_qcov: minimum bidirectional coverage threshold (default: 80.0%)
        
    Returns:
        sim_matrix:      np.ndarray of pairwise percent identity (n x n), values 0-100
        family_labels:   list of full FASTA header names (length n)
        short_labels:    list of shortened display labels (length n)
        summary_tsv:     path to the similarity matrix TSV file
        merged_fastas:   list of final merged family FASTA paths
        merge_summary:   path to the merge summary TSV
    """
    logger.info("--- Cross-Chromosome Monomer Family Comparison ---")
    
    cross_dir = os.path.join(out_dir, "cross_chrom_comparison")
    os.makedirs(cross_dir, exist_ok=True)
    
    # 1. Merge all representative sequences into one FASTA
    all_reps_fasta = os.path.join(cross_dir, "all_representatives.fasta")
    family_labels = []
    
    with open(all_reps_fasta, 'w') as f_out:
        for chrom in sorted(per_chrom_results.keys(), key=natural_sort_key):
            for fasta in per_chrom_results[chrom]:
                with open(fasta, 'r') as f_in:
                    for line in f_in:
                        if line.startswith('>'):
                            family_labels.append(line.strip().lstrip('>'))
                        f_out.write(line)
    
    n = len(family_labels)
    short_labels = [shorten_family_label(l) for l in family_labels]
    logger.info(f"Merged {n} representative sequences for cross-chromosome comparison.")
    logger.info(f"Bidirectional coverage threshold: {min_qcov}%")
    
    if n < 2:
        logger.warning("Less than 2 families found, skipping cross-chromosome comparison.")
        sim_matrix = np.array([[100.0]]) if n == 1 else np.array([])
        return sim_matrix, family_labels, short_labels, None, [], None

    # 2. All-vs-all BLAST
    blast_out = os.path.join(cross_dir, "all_vs_all.blast")
    q_reps = shlex.quote(all_reps_fasta)
    q_blast = shlex.quote(blast_out)
    cmd = (
        f"blastn -task blastn -query {q_reps} -subject {q_reps} "
        f"-out {q_blast} "
        f"-outfmt '6 qseqid sseqid pident length qlen slen evalue bitscore' "
        f"-num_threads {threads} -evalue 1e-3 -dust no"
    )
    subprocess.run(cmd, shell=True, check=True)
    logger.info("All-vs-all BLAST completed.")
    
    # 3. Build similarity matrix with bidirectional coverage filtering
    #    For each pair, we require BOTH query and subject coverage >= threshold.
    #    This ensures symmetric treatment regardless of sequence length differences.
    #    We aggregate all valid HSPs per pair and keep the max pident.
    label_to_idx = {label: i for i, label in enumerate(family_labels)}
    sim_matrix = np.zeros((n, n))
    np.fill_diagonal(sim_matrix, 100.0)
    
    # Collect best pident per ordered pair (considering bidirectional coverage)
    pair_best = {}  # (i, j) -> best pident where i < j
    
    with open(blast_out, 'r') as f:
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) < 8:
                continue
            q_name, s_name = cols[0], cols[1]
            pident = float(cols[2])
            aln_len = int(cols[3])
            q_len = int(cols[4])
            s_len = int(cols[5])
            
            if q_name == s_name:
                continue
            if q_name not in label_to_idx or s_name not in label_to_idx:
                continue
            
            # Bidirectional coverage: require both query and subject to be well-covered
            query_cov = (aln_len / q_len) * 100 if q_len > 0 else 0
            subject_cov = (aln_len / s_len) * 100 if s_len > 0 else 0
            min_cov = min(query_cov, subject_cov)
            
            if min_cov < min_qcov:
                continue
                
            i, j = label_to_idx[q_name], label_to_idx[s_name]
            # Normalize to ordered pair (smaller index first) for consistent aggregation
            pair_key = (min(i, j), max(i, j))
            pair_best[pair_key] = max(pair_best.get(pair_key, 0), pident)
    
    # Fill symmetric matrix from aggregated pairs
    for (i, j), best_pident in pair_best.items():
        sim_matrix[i][j] = best_pident
        sim_matrix[j][i] = best_pident
    
    # 4. Write similarity matrix TSV
    summary_tsv = os.path.join(cross_dir, "similarity_matrix.tsv")
    with open(summary_tsv, 'w') as f:
        f.write("Family\t" + "\t".join(short_labels) + "\n")
        for i, label in enumerate(short_labels):
            row = "\t".join(f"{sim_matrix[i][j]:.1f}" for j in range(n))
            f.write(f"{label}\t{row}\n")
    logger.info(f"Identity matrix saved to: {summary_tsv}")
    
    # 5. Write cross-chromosome pairwise summary
    pair_tsv = os.path.join(cross_dir, "cross_chrom_pairs.tsv")
    with open(pair_tsv, 'w') as f:
        f.write("Family_A\tFamily_B\tIdentity(%)\tRelationship\n")
        for i in range(n):
            chrom_i = family_labels[i].split('_Family_')[0] if '_Family_' in family_labels[i] else ""
            for j in range(i + 1, n):
                chrom_j = family_labels[j].split('_Family_')[0] if '_Family_' in family_labels[j] else ""
                if chrom_i == chrom_j:
                    continue
                sim = sim_matrix[i][j]
                if sim >= 80:
                    relationship = "Same_family"
                elif sim >= 60:
                    relationship = "Divergent"
                else:
                    relationship = "Different"
                f.write(f"{short_labels[i]}\t{short_labels[j]}\t{sim:.1f}\t{relationship}\n")
    
    logger.info(f"Cross-chromosome pairwise comparisons saved to: {pair_tsv}")
    
    # 6. Log summary
    logger.info("Cross-chromosome relationship summary:")
    reported = False
    for i in range(n):
        chrom_i = family_labels[i].split('_Family_')[0] if '_Family_' in family_labels[i] else ""
        for j in range(i + 1, n):
            chrom_j = family_labels[j].split('_Family_')[0] if '_Family_' in family_labels[j] else ""
            if chrom_i != chrom_j and sim_matrix[i][j] >= 80:
                logger.info(f"  {short_labels[i]} <-> {short_labels[j]}: {sim_matrix[i][j]:.1f}% identity (Same family)")
                reported = True
    
    if not reported:
        logger.info("  No cross-chromosome same-family pairs found (threshold: 80% identity).")
        logger.info("  All chromosomes appear to use distinct monomer families.")
    
    # 7. Merge families into final unified set
    merged_fastas, merge_summary = merge_families(
        sim_matrix, family_labels, short_labels, per_chrom_results,
        out_dir, merge_threshold=80.0, logger=logger
    )
    
    return sim_matrix, family_labels, short_labels, summary_tsv, merged_fastas, merge_summary
