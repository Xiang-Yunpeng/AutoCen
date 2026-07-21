#!/usr/bin/env python3
import os
import sys
import argparse
from utils.logger import setup_logger
from utils.chrom_utils import get_target_chroms, natural_sort_key
from modules.trash_runner import check_dependencies, run_trash, find_trash_files
from modules.array_parser import extract_top_arrays
from modules.monomer_core import extract_and_cluster, extract_and_cluster_per_chrom
from modules.sliding_window import run_sliding_window
from modules.visualization import plot_multi_tracks, plot_refined_report
from modules.cross_chrom_analysis import cross_chromosome_analysis
from modules.te_scanner import run_scan_te

__version__ = "1.0.0"

def main():
    parser = argparse.ArgumentParser(description="AutoCen: Automated Centromere Identification and Multi-omics Visualization Toolkit")
    parser.add_argument("--version", action="version", version=f"AutoCen v{__version__}")
    subparsers = parser.add_subparsers(dest="command", help="Available subcommands")

    # ================= Subcommand 1: scan_te =================
    parser_scan = subparsers.add_parser("scan_te", help="Pre-scan TE composition from RepeatMasker .out files")
    parser_scan.add_argument("-o", "--rm_out", required=True, help="Path to the RepeatMasker .out file")
    parser_scan.add_argument("-f", "--fai", required=True, help="Path to the genome .fai file")

    # ================= Subcommand 2: centromere =================
    parser_cen = subparsers.add_parser("centromere", help="Task 1: Identify centromere monomers and generate basic density plot")
    
    # Basic Arguments
    parser_cen.add_argument("--genome", required=True, help="Path to input genome FASTA")
    parser_cen.add_argument("--work_dir", required=True, help="Working directory for output")
    
    # TRASH Arguments
    parser_cen.add_argument("--trash_path", help="Path to TRASH.R script")
    parser_cen.add_argument("--trash_dir", help="Path to existing TRASH output directory")
    
    # Execution Arguments
    parser_cen.add_argument("--threads", type=int, default=20, help="Number of threads to use (default: 20)")
    parser_cen.add_argument("--window", type=int, default=10000, help="Sliding window size in bp for basic plot (default: 10000)")
    
    # Chromosome Selection
    parser_cen.add_argument("--chrom_num", type=int, default=None, 
                            help="Only analyze the first N sequences in the genome FASTA (e.g. --chrom_num 12 for rice). "
                                 "Chromosomes are typically listed before unplaced contigs/scaffolds in standard assemblies. "
                                 "If not set, all sequences are analyzed.")
    
    # Refined Mode
    parser_cen.add_argument("--refined", action='store_true', default=False,
                            help="Enable refined mode: cluster monomers independently per chromosome, "
                                 "then perform cross-chromosome family comparison. Recommended for species "
                                 "where different chromosomes may use distinct centromeric monomer families.")
    
    # Extraction Arguments
    parser_cen.add_argument("--top_arrays", type=int, default=1, help="Top N longest arrays to extract per chromosome (default: 1)")
    parser_cen.add_argument("--min_array_len", type=int, default=50000, help="Minimum array length in bp (default: 50000)")
    parser_cen.add_argument("--top_monomers", type=int, default=3, help="Top X monomer families to extract (default: 3)")
    parser_cen.add_argument("--min_monomer", type=int, default=10, 
                            help="[Refined mode only] Minimum number of monomers per chromosome to run clustering (default: 10). "
                                 "Chromosomes below this threshold are skipped with a warning. "
                                 "Ignored in Standard mode where all monomers are pooled globally.")
    parser_cen.add_argument("--min_monomer_len", type=int, default=50, 
                            help="Minimum monomer sequence length in bp (default: 50). "
                                 "Shorter repeats (e.g. microsatellites, telomeric repeats) are excluded before clustering.")
    parser_cen.add_argument("--min_cluster_size", type=int, default=5, 
                            help="Minimum CD-HIT cluster size to output a monomer family (default: 5). "
                                 "Clusters smaller than this are discarded.")
    
    # Alignment & Clustering Parameters
    parser_cen.add_argument("--cdhit_id", type=float, default=0.8, help="CD-HIT clustering identity threshold (default: 0.8)")
    parser_cen.add_argument("--cdhit_aS", type=float, default=0.0, help="CD-HIT alignment coverage for shorter sequence (default: 0.0)")
    parser_cen.add_argument("--cdhit_aL", type=float, default=0.0, help="CD-HIT alignment coverage for longer sequence (default: 0.0)")
    
    parser_cen.add_argument("--blast_task", type=str, choices=['megablast', 'blastn', 'dc-megablast'], default='megablast', help="BLAST task to execute (default: megablast)")
    parser_cen.add_argument("--blast_identity", type=float, default=80.0, help="BLAST minimum percent identity (default: 80)")
    parser_cen.add_argument("--blast_qcov", type=float, default=0, help="BLAST minimum query coverage per HSP (default: 0)")
    parser_cen.add_argument("--blast_evalue", type=str, default="1e-10", help="BLAST e-value threshold (default: 1e-10)")
    parser_cen.add_argument("--cross_qcov", type=float, default=80.0, 
                            help="Minimum bidirectional coverage (%%) for cross-chromosome comparison in refined mode (default: 80.0). "
                                 "Both query and subject coverage must exceed this threshold.")

    # ================= Subcommand 3: plot =================
    parser_plot = subparsers.add_parser("plot", help="Task 2: Generate highly customizable multi-omics visualization tracks")
    
    # Basic & Monomer Arguments
    parser_plot.add_argument("--genome", required=True, help="Path to input genome FASTA")
    parser_plot.add_argument("--work_dir", required=True, help="Working directory for output")
    parser_plot.add_argument("--monomers", nargs='+', required=True, 
                             help="One or more monomer FASTA files. Supports custom track names (eg. CEN1:Fam_1.fasta CEN2:Fam_2.fasta)")
    
    # Execution & Visualization Arguments
    parser_plot.add_argument("--threads", type=int, default=20, help="Number of threads to use (default: 20)")
    parser_plot.add_argument("--window", type=int, default=10000, help="Sliding window size in bp (default: 10000)")
    parser_plot.add_argument("--region", help="Specific regions to plot (eg. 'chr1,chr2' or 'chr1:16M-18M'). Overrides --chrom_num.")
    parser_plot.add_argument("--chrom_num", type=int, default=None, help="Number of top contigs to visualize if --region is not set.")
    parser_plot.add_argument("--fig_width", type=float, default=10.0, help="Width of the output figure in inches (default: 10.0)")
    parser_plot.add_argument("--track_height", type=float, default=1.8, help="Height of each standard track in inches (default: 1.8)")
    
    # Alignment Parameters for Plotting
    parser_plot.add_argument("--blast_task", type=str, choices=['megablast', 'blastn', 'dc-megablast'], default='megablast', help="BLAST task to execute (default: megablast)")
    parser_plot.add_argument("--blast_identity", type=float, default=80.0, help="BLAST minimum percent identity (default: 80)")
    parser_plot.add_argument("--blast_qcov", type=float, default=0, help="BLAST minimum query coverage per HSP (default: 0)")
    parser_plot.add_argument("--blast_evalue", type=str, default="1e-10", help="BLAST e-value threshold (default: 1e-10)")

    # Multi-omics Tracks
    parser_plot.add_argument("--rm_out", help="Optional: Path to RepeatMasker .out file for TE density track")
    parser_plot.add_argument("--te_mapping", nargs='+', 
                             help="Optional: TE mapping rules. Format: TrackName:TE1,TE2 (eg. LINE/Rex:Rex-Babar,LINE/L2 LTR/Gypsy:Gypsy)")
    parser_plot.add_argument("--gff", help="Optional: Path to GFF3 file for Gene density track")
    parser_plot.add_argument("--hic", help="Optional: Path to .hic file for Hi-C heatmap track")
    parser_plot.add_argument("--hic_res", type=int, default=250000, help="Resolution for Hi-C heatmap (default: 250000)")
    parser_plot.add_argument("--hic_vmax", type=int, default=20, help="Maximum value for Hi-C heatmap color scale (default: 20)")

    args = parser.parse_args()
    
    # ---------------- Execution Logic ----------------
    if args.command == "scan_te":
        run_scan_te(args.rm_out, args.fai)
        sys.exit(0)
        
    elif args.command == "centromere":
        work_dir = os.path.abspath(args.work_dir)
        dirs = {
            "step1": os.path.join(work_dir, "1_TRASH"),
            "step2": os.path.join(work_dir, "2_candidate_arrays"),
            "step3": os.path.join(work_dir, "3_monomer_clustering"),
            "step4": os.path.join(work_dir, "4_sliding_windows"),
            "step5": os.path.join(work_dir, "5_visualization")
        }
        for d in dirs.values(): os.makedirs(d, exist_ok=True)
            
        logger = setup_logger(work_dir)
        
        mode_label = "Refined" if args.refined else "Standard"
        logger.info(f"Starting AutoCen Task 1: Centromere Identification ({mode_label} Mode)...")
        
        check_dependencies(args.genome, args.trash_path, args.trash_dir, logger)
        
        # Resolve target chromosomes
        target_chroms = get_target_chroms(args.genome, args.chrom_num, logger)
        
        if args.trash_dir:
            trash_arrays_csv, trash_repeats_csv = find_trash_files(args.trash_dir, logger)
        else:
            trash_arrays_csv, trash_repeats_csv = run_trash(args.genome, dirs["step1"], args.trash_path, args.threads, logger)
            
        candidate_bed = os.path.join(dirs["step2"], "candidate_arrays.bed")
        extract_top_arrays(trash_arrays_csv, candidate_bed, args.top_arrays, args.min_array_len, logger, 
                           target_chroms=target_chroms)
        
        if args.refined:
            # ====== Refined Mode Pipeline ======
            
            # Step 3: Per-chromosome clustering
            per_chrom_results = extract_and_cluster_per_chrom(
                trash_repeats_csv, candidate_bed, dirs["step3"],
                args.cdhit_id, args.cdhit_aS, args.cdhit_aL, args.threads, args.top_monomers, logger,
                min_monomer=args.min_monomer, min_monomer_len=args.min_monomer_len,
                min_cluster_size=args.min_cluster_size
            )
            
            # Flatten all per-chromosome FASTAs for sliding window
            all_monomer_fastas = []
            for chrom in sorted(per_chrom_results.keys(), key=natural_sort_key):
                all_monomer_fastas.extend(per_chrom_results[chrom])
            
            # Step 4: Sliding window (runs on all families at once)
            density_tsvs, genome_sizes = run_sliding_window(
                args.genome, all_monomer_fastas, dirs["step4"], args.window, args.threads,
                args.blast_task, args.blast_identity, args.blast_evalue, args.blast_qcov, logger
            )
            
            # Step 4.5: Cross-chromosome comparison and family merging
            sim_matrix, family_labels, short_labels, summary_tsv, merged_fastas, merge_summary = cross_chromosome_analysis(
                per_chrom_results, dirs["step3"], args.threads, logger, min_qcov=args.cross_qcov
            )
            
            # Restructure density_tsvs into per-chromosome dict for visualization
            # Key format from sliding_window: "chr1_Family_1_size100" (basename of fasta without .fasta)
            per_chrom_density_tsvs = {}
            for fam_name, tsv_path in density_tsvs.items():
                if '_Family_' in fam_name:
                    chrom = fam_name.split('_Family_')[0]
                else:
                    chrom = "unknown"
                if chrom not in per_chrom_density_tsvs:
                    per_chrom_density_tsvs[chrom] = {}
                per_chrom_density_tsvs[chrom][fam_name] = tsv_path
            
            # Step 5: Refined visualization
            out_pdf = os.path.join(dirs["step5"], "autocen_refined_report.pdf")
            plot_refined_report(
                per_chrom_density_tsvs=per_chrom_density_tsvs,
                genome_sizes_file=genome_sizes,
                sim_matrix=sim_matrix,
                short_labels=short_labels,
                family_labels=family_labels,
                out_pdf=out_pdf,
                window_size=args.window,
                logger=logger
            )
            
            logger.info("Task 1 Completed (Refined Mode)! Per-chromosome clustering and cross-chromosome analysis done.")
            
        else:
            # ====== Standard Mode Pipeline (unchanged) ======
            
            monomer_fastas = extract_and_cluster(
                trash_repeats_csv, candidate_bed, dirs["step3"], 
                args.cdhit_id, args.cdhit_aS, args.cdhit_aL, args.threads, args.top_monomers, logger,
                min_monomer_len=args.min_monomer_len, min_cluster_size=args.min_cluster_size
            )
            
            density_tsvs, genome_sizes = run_sliding_window(
                args.genome, monomer_fastas, dirs["step4"], args.window, args.threads,
                args.blast_task, args.blast_identity, args.blast_evalue, args.blast_qcov, logger
            )
            
            out_pdf = os.path.join(dirs["step5"], "autocen_monomer_basic_report.pdf")
            plot_multi_tracks(
                density_tsvs=density_tsvs, genome_sizes_file=genome_sizes, out_pdf=out_pdf, 
                region_str=None, chrom_num=args.chrom_num, window_size=args.window, rm_out=None, 
                te_mapping=None, gff_file=None, hic_file=None, hic_res=None, logger=logger,
                separate_monomer_tracks=True
            )
            logger.info("Task 1 Completed! Centromere monomers extracted and basic plot generated.")

    elif args.command == "plot":
        if (args.rm_out and not args.te_mapping) or (args.te_mapping and not args.rm_out):
            print("Error: --rm_out and --te_mapping must be provided together.")
            sys.exit(1)
            
        work_dir = os.path.abspath(args.work_dir)
        dirs = {
            "step4": os.path.join(work_dir, "4_sliding_windows"),
            "step5": os.path.join(work_dir, "5_visualization")
        }
        for d in dirs.values(): os.makedirs(d, exist_ok=True)
            
        logger = setup_logger(work_dir)
        logger.info("Starting AutoCen Task 2: Multi-omics Plotting...")
        
        monomer_fastas = []
        monomer_names_map = {}
        
        for m in args.monomers:
            if ':' in m:
                custom_name, filepath = m.split(':', 1)
                filepath = os.path.abspath(filepath)
                monomer_fastas.append(filepath)
                base_key = os.path.basename(filepath).replace(".fasta", "")
                monomer_names_map[base_key] = custom_name
            else:
                filepath = os.path.abspath(m)
                monomer_fastas.append(filepath)
                base_key = os.path.basename(filepath).replace(".fasta", "")
                monomer_names_map[base_key] = base_key.split('_size')[0]
                
        for f in monomer_fastas:
            if not os.path.exists(f):
                logger.error(f"Provided monomer FASTA file not found: {f}")
                sys.exit(1)
        
        density_tsvs, genome_sizes = run_sliding_window(
            args.genome, monomer_fastas, dirs["step4"], args.window, args.threads,
            args.blast_task, args.blast_identity, args.blast_evalue, args.blast_qcov, logger
        )
        
        out_pdf = os.path.join(dirs["step5"], "autocen_multi_track_report.pdf")
        plot_multi_tracks(
            density_tsvs=density_tsvs, genome_sizes_file=genome_sizes, out_pdf=out_pdf, 
            region_str=args.region, chrom_num=args.chrom_num, window_size=args.window,
            rm_out=args.rm_out, te_mapping=args.te_mapping, gff_file=args.gff,           
            hic_file=args.hic, hic_res=args.hic_res, logger=logger,
            fig_width=args.fig_width, track_height=args.track_height, monomer_names_map=monomer_names_map,
            hic_vmax=args.hic_vmax
        )
        logger.info("Task 2 Completed! High-quality multi-omics plots generated successfully.")
        
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
