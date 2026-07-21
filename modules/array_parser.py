import pandas as pd
import os
import sys

def extract_top_arrays(arrays_csv, output_bed, top_n, min_length, logger, target_chroms=None):
    logger.info("--- Step 2: Parsing Candidate Arrays ---")
    logger.info(f"Reading arrays file: {arrays_csv}")
    
    try:
        df = pd.read_csv(arrays_csv)
        df['array_length'] = df['end'] - df['start']
        logger.info(f"Total arrays found: {len(df)}")
        
        # Filter by target chromosomes (if specified)
        if target_chroms is not None:
            before = len(df)
            df = df[df['seqID'].isin(target_chroms)]
            logger.info(f"Chromosome filter applied: {before} -> {len(df)} arrays "
                        f"(kept {len(target_chroms)} chromosomes)")

        # Filter by minimum length
        if min_length > 0:
            df = df[df['array_length'] >= min_length]
            logger.info(f"Arrays passing length filter (>={min_length}bp): {len(df)}")
            
        if df.empty:
            logger.error("No candidate arrays found after filtering!")
            sys.exit(1)

        # Sort by length descending, then group by seqID and take top_n
        df_sorted = df.sort_values(by=['seqID', 'array_length'], ascending=[True, False])
        top_arrays_df = df_sorted.groupby('seqID').head(top_n)
        final_df = top_arrays_df.sort_values(by=['seqID', 'start'])
        
        logger.info(f"Extracted top {top_n} arrays per chromosome. Total candidates: {len(final_df)}")
        
        # Write to BED
        bed_output = final_df[['seqID', 'start', 'end', 'array_length']]
        bed_output.to_csv(output_bed, sep='\t', index=False, header=False)
        logger.info(f"Candidate arrays saved to: {output_bed}")
        
    except Exception as e:
        logger.error(f"Error parsing array file: {e}")
        sys.exit(1)
