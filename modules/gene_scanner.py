import sys
from collections import defaultdict

def extract_gene_coordinates(gff_path, logger):
    """
    Parse a standard GFF3 file and extract coordinates for 'gene' features.
    Structure: { "chr1": [(start, end), ...], ... }
    """
    extracted_data = defaultdict(list)
    gene_count = 0
    
    try:
        with open(gff_path, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                    
                parts = line.strip().split('\t')
                # A standard GFF3 line should have 9 columns
                if len(parts) < 9:
                    continue
                    
                chrom = parts[0]
                feature_type = parts[2]
                
                # We only care about the main 'gene' feature for density
                if feature_type == 'gene':
                    try:
                        start = int(parts[3])
                        end = int(parts[4])
                        extracted_data[chrom].append((start, end))
                        gene_count += 1
                    except ValueError:
                        continue
                        
        logger.info(f"Successfully extracted {gene_count} genes from GFF file.")
        
    except FileNotFoundError:
        logger.error(f"GFF file not found: {gff_path}")
        sys.exit(1)
        
    return extracted_data
