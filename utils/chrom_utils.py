import re
import sys

def natural_sort_key(s):
    """Generate a sort key that orders strings in natural/human order.
    
    'chr1, chr2, chr10' instead of 'chr1, chr10, chr2'.
    Splits the string into text and numeric chunks, converting numeric parts to int.
    """
    return [int(c) if c.isdigit() else c.lower() for c in re.split(r'(\d+)', str(s))]

def get_target_chroms(genome_fa, chrom_num, logger):
    """
    Read FASTA headers in order and return the first `chrom_num` sequence names.
    
    In standard genome assemblies, chromosomes are listed first (usually sorted
    by size or ID), followed by unplaced contigs/scaffolds. By taking the first N
    sequences, we effectively select only chromosomal-level sequences.
    
    Returns None if chrom_num is not set (meaning: analyze all sequences).
    """
    if chrom_num is None or chrom_num <= 0:
        return None

    chrom_names = []
    try:
        with open(genome_fa, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    # Take the first whitespace-delimited token after '>'
                    seq_name = line[1:].strip().split()[0]
                    chrom_names.append(seq_name)
                    if len(chrom_names) >= chrom_num:
                        break
    except FileNotFoundError:
        logger.error(f"Genome FASTA not found: {genome_fa}")
        sys.exit(1)

    if len(chrom_names) < chrom_num:
        logger.warning(
            f"Requested {chrom_num} chromosomes, but genome only contains "
            f"{len(chrom_names)} sequences. Using all {len(chrom_names)}."
        )

    logger.info(f"Target chromosomes ({len(chrom_names)}): {', '.join(chrom_names)}")
    return chrom_names
