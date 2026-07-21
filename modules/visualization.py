import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, LinearSegmentedColormap
from matplotlib.backends.backend_pdf import PdfPages
from collections import defaultdict
from matplotlib.ticker import FuncFormatter

from modules.gene_scanner import extract_gene_coordinates
from utils.chrom_utils import natural_sort_key

# ================= Density Processing Functions =================
def parse_te_mapping(te_args, logger):
    reverse_mapping, track_names = {}, []
    for item in te_args:
        if ':' not in item:
            logger.warning(f"Skipping invalid TE mapping argument: '{item}'")
            continue
        track_name, te_str = item.split(':', 1)
        track_name = track_name.strip()
        track_names.append(track_name)
        for te in [t.strip() for t in te_str.split(',')]:
            reverse_mapping[te] = track_name
    return reverse_mapping, track_names

def merge_intervals(intervals):
    if not intervals: return []
    intervals.sort()
    merged = []
    curr_start, curr_end = intervals[0]
    for next_start, next_end in intervals[1:]:
        if next_start <= curr_end:
            curr_end = max(curr_end, next_end)
        else:
            merged.append((curr_start, curr_end))
            curr_start, curr_end = next_start, next_end
    merged.append((curr_start, curr_end))
    return merged

def extract_te_coordinates(out_path, reverse_mapping, logger):
    extracted_data = defaultdict(lambda: defaultdict(list))
    try:
        with open(out_path, 'r') as f:
            for _ in range(3): next(f)
            for line in f:
                parts = line.strip().split()
                if len(parts) < 11: continue
                try:
                    chrom, start, end, raw_label = parts[4], int(parts[5]), int(parts[6]), parts[10]
                except ValueError: continue
                track_name = reverse_mapping.get(raw_label)
                if track_name:
                    extracted_data[track_name][chrom].append((start, end))
    except FileNotFoundError:
        logger.error(f"RepeatMasker .out file not found: {out_path}")
        sys.exit(1)
    return extracted_data

def calculate_window_density(coords, chrom_length, window_size):
    """Generic function to calculate density for TEs or Genes."""
    merged_coords = merge_intervals(coords)
    num_bins = int(np.ceil(chrom_length / window_size))
    if num_bins == 0: return np.array([]), np.array([])
    coverage = np.zeros(num_bins)
    for start, end in merged_coords:
        start_bin = min(max(0, start // window_size), num_bins - 1)
        end_bin = min(max(0, end // window_size), num_bins - 1)
        if start_bin == end_bin:
            coverage[start_bin] += (end - start + 1)
        else:
            coverage[start_bin] += ((start_bin + 1) * window_size - start + 1)
            for b in range(start_bin + 1, end_bin): coverage[b] += window_size
            coverage[end_bin] += (end - end_bin * window_size)
    bin_sizes = np.full(num_bins, float(window_size))
    bin_sizes[-1] = chrom_length - (num_bins - 1) * window_size
    return np.arange(num_bins) * window_size + (bin_sizes / 2), (coverage / bin_sizes) * 100

# ================= Region Parsing =================
def parse_num(val_str):
    val_str = val_str.upper().strip()
    if val_str.endswith('K'): return int(float(val_str[:-1]) * 1_000)
    if val_str.endswith('M'): return int(float(val_str[:-1]) * 1_000_000)
    return int(val_str)

def parse_target_regions(region_str, chrom_sizes, chrom_num, logger):
    regions = []
    if region_str:
        parts = [p.strip() for p in region_str.split(',')]
        for p in parts:
            if ':' in p:
                chrom, coords = p.split(':', 1)
                coords = coords.replace(',', '-') 
                start_str, end_str = coords.split('-')
                start, end = parse_num(start_str), parse_num(end_str)
                max_len = chrom_sizes.get(chrom, end)
                start, end = max(0, start), min(max_len, end)
                regions.append({'chrom': chrom, 'start': start, 'end': end})
            else:
                chrom = p
                if chrom in chrom_sizes:
                    regions.append({'chrom': chrom, 'start': 0, 'end': chrom_sizes[chrom]})
                else:
                    logger.warning(f"Chromosome {chrom} not found in genome.sizes")
    else:
        ordered_chroms = list(chrom_sizes.keys())
        target = ordered_chroms[:chrom_num] if chrom_num else ordered_chroms
        for chrom in target:
            regions.append({'chrom': chrom, 'start': 0, 'end': chrom_sizes[chrom]})
    return regions

# ================= Formatters & Styling Helpers =================
def m_fmt(x, pos):
    """Format X-axis ticks to show 'M' suffix (e.g., 19.5M)"""
    return f"{x:.1f}M"

def style_axis(ax, label_text):
    """Apply common publication-style minimalist aesthetics to a track axis"""
    ax.set_ylabel(label_text, rotation=270, labelpad=30, va='center', fontsize=11)
    ax.yaxis.set_label_position("right")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.tick_params(bottom=False, labelbottom=False)
    ax.margins(x=0)
    ax.grid(False)

# ================= Main Multi-Track Plotter (Original) =================
def plot_multi_tracks(density_tsvs, genome_sizes_file, out_pdf, region_str, chrom_num, window_size, rm_out, te_mapping, gff_file, hic_file, hic_res, logger, fig_width=10.0, track_height=1.8, monomer_names_map=None, separate_monomer_tracks=False, hic_vmax=20):
    logger.info("--- Step 5: Multi-Track Visualization ---")
    
    # 1. Parse Chromosomes & Regions
    chrom_sizes = {}
    with open(genome_sizes_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            chrom_sizes[parts[0]] = int(parts[1])
            
    target_regions = parse_target_regions(region_str, chrom_sizes, chrom_num, logger)
    logger.info(f"Visualizing {len(target_regions)} specified region(s).")

    # 2. Process TE Data
    te_data = {}
    te_tracks = []
    if rm_out and te_mapping:
        logger.info("Processing TE annotations...")
        reverse_mapping, te_tracks = parse_te_mapping(te_mapping, logger)
        te_data = extract_te_coordinates(rm_out, reverse_mapping, logger)

    # 3. Process Gene Data
    gene_data = {}
    if gff_file:
        logger.info("Processing Gene annotations...")
        gene_data = extract_gene_coordinates(gff_file, logger)
        
    # 4. Check Hi-C
    hic_obj = None
    if hic_file:
        try:
            import hicstraw
            logger.info(f"Loading Hi-C file: {hic_file} at resolution {hic_res} bp")
            hic_obj = hicstraw.HiCFile(hic_file)
        except ImportError:
            logger.error("hicstraw is not installed! Cannot plot Hi-C track.")
            hic_file = None

    # 5. Plotting Setup
    if monomer_names_map is None:
        monomer_names_map = {}
        
    monomer_colors = ['#f8766d', '#00bfc4', '#7cae00', '#c77cff', '#ff7f00', '#ffff33']
    te_colors = ['#00bfc4', '#c77cff', '#e6ab02', '#a6761d'] 
    gene_color = '#7cae00' 
    
    monomer_dfs = {name: pd.read_csv(tsv, sep='\t', names=['chrom', 'start', 'end', 'count']) 
                   for name, tsv in density_tsvs.items()}

    with PdfPages(out_pdf) as pdf:
        for region in target_regions:
            chrom = region['chrom']
            r_start, r_end = region['start'], region['end']
            r_start_mb = r_start / 1_000_000.0
            r_end_mb = r_end / 1_000_000.0
            c_len_bp = chrom_sizes[chrom]
            
            title_text = f"{chrom}" if r_start == 0 and r_end == c_len_bp else f"{chrom}:{r_start}-{r_end}"
            logger.info(f"Plotting tracks for {title_text} ...")
            
            num_monomer_tracks = len(monomer_dfs) if separate_monomer_tracks else 1
            num_te_tracks = len(te_tracks)
            has_gene = 1 if gff_file else 0
            has_hic = 1 if hic_obj else 0
            
            total_tracks = num_monomer_tracks + num_te_tracks + has_gene + has_hic
            height_ratios = [1] * num_monomer_tracks + [1] * num_te_tracks + [1] * has_gene
            
            if has_hic: 
                height_ratios.append(5) 
                fig_height = track_height * (num_monomer_tracks + num_te_tracks + has_gene) + (track_height * 2.5)
            else:
                fig_height = track_height * (num_monomer_tracks + num_te_tracks + has_gene)
                
            fig, axes = plt.subplots(total_tracks, 1, figsize=(fig_width, fig_height), 
                                     sharex=True, gridspec_kw={'height_ratios': height_ratios})
            
            if total_tracks == 1: axes = [axes]
            axes_iter = iter(axes)
            
            # --- Monomer Track(s) ---
            if separate_monomer_tracks:
                # Each monomer family gets its own independent track
                for j, (fam_name, df) in enumerate(monomer_dfs.items()):
                    ax = next(axes_iter)
                    if j == 0:
                        ax.set_title(title_text, loc='left', fontsize=12, fontweight='bold')
                    
                    color = monomer_colors[j % len(monomer_colors)]
                    label_name = monomer_names_map.get(fam_name, fam_name.split('_size')[0])
                    
                    sub_df = df[df['chrom'] == chrom]
                    max_y = 0
                    if not sub_df.empty:
                        view_df = sub_df[(sub_df['start'] >= r_start) & (sub_df['end'] <= r_end)]
                        if not view_df.empty and view_df['count'].sum() > 0:
                            x_mb = sub_df['start'] / 1_000_000.0
                            ax.fill_between(x_mb, sub_df['count'], color=color, alpha=0.4)
                            ax.plot(x_mb, sub_df['count'], color=color, lw=1.2)
                            max_y = view_df['count'].max()
                    
                    style_axis(ax, label_name)
                    y_upper = max_y * 1.1 if max_y > 0 else 10
                    ax.set_ylim(-y_upper * 0.05, y_upper)
            else:
                # All monomer families overlaid on a single track
                ax_mono = next(axes_iter)
                ax_mono.set_title(title_text, loc='left', fontsize=12, fontweight='bold')
                max_mono_y = 0
                
                mono_labels = []
                for j, (fam_name, df) in enumerate(monomer_dfs.items()):
                    sub_df = df[df['chrom'] == chrom]
                    if not sub_df.empty:
                        view_df = sub_df[(sub_df['start'] >= r_start) & (sub_df['end'] <= r_end)]
                        if not view_df.empty and view_df['count'].sum() > 0:
                            x_mb = sub_df['start'] / 1_000_000.0
                            label_name = monomer_names_map.get(fam_name, fam_name.split('_size')[0])
                            
                            ax_mono.plot(x_mb, sub_df['count'], color=monomer_colors[j % len(monomer_colors)], 
                                         lw=1.2, label=label_name)
                            max_mono_y = max(max_mono_y, view_df['count'].max())
                            mono_labels.append(label_name)
                            
                track_label = "\n".join(mono_labels) if mono_labels else "Monomer"
                style_axis(ax_mono, track_label)
                
                y_upper_mono = max_mono_y * 1.1 if max_mono_y > 0 else 10
                y_lower_mono = -y_upper_mono * 0.05
                ax_mono.set_ylim(y_lower_mono, y_upper_mono)
                
                if ax_mono.get_legend_handles_labels()[0] and len(mono_labels) > 1:
                    ax_mono.legend(loc="upper left", fontsize=8, frameon=False)

            # --- Track 2..N: TEs ---
            for j, t_name in enumerate(te_tracks):
                ax_te = next(axes_iter)
                color = te_colors[j % len(te_colors)]
                coords = te_data[t_name].get(chrom, [])
                
                max_te_y = 0
                if coords:
                    x_pos, y_dens = calculate_window_density(coords, c_len_bp, window_size)
                    x_pos_mb = x_pos / 1_000_000.0
                    
                    mask = (x_pos_mb >= r_start_mb) & (x_pos_mb <= r_end_mb)
                    if any(mask): max_te_y = max(y_dens[mask])
                    
                    ax_te.plot(x_pos_mb, y_dens, color=color, linewidth=1.2)
                    
                style_axis(ax_te, t_name.replace('/', '/\n')) 
                
                y_upper_te = max_te_y * 1.1 if max_te_y > 0 else 100
                y_lower_te = -y_upper_te * 0.05
                ax_te.set_ylim(y_lower_te, y_upper_te)
                
            # --- Track N+1: Genes ---
            if has_gene:
                ax_gene = next(axes_iter)
                coords = gene_data.get(chrom, [])
                
                max_gene_y = 0
                if coords:
                    x_pos, y_dens = calculate_window_density(coords, c_len_bp, window_size)
                    x_pos_mb = x_pos / 1_000_000.0
                    
                    mask = (x_pos_mb >= r_start_mb) & (x_pos_mb <= r_end_mb)
                    if any(mask): max_gene_y = max(y_dens[mask])
                        
                    ax_gene.plot(x_pos_mb, y_dens, color=gene_color, linewidth=1.2)
                    
                style_axis(ax_gene, "Gene")
                
                y_upper_gene = max_gene_y * 1.1 if max_gene_y > 0 else 100
                y_lower_gene = -y_upper_gene * 0.05
                ax_gene.set_ylim(y_lower_gene, y_upper_gene)

            # --- Track Last: Hi-C ---
            if has_hic:
                ax_hic = next(axes_iter)
                style_axis(ax_hic, "Hi-C")
                try:
                    mzd = hic_obj.getMatrixZoomData(chrom, chrom, "observed", "NONE", "BP", hic_res)
                    matrix = mzd.getRecordsAsMatrix(r_start, r_end, r_start, r_end)
                    
                    if matrix is not None and matrix.size > 0:
                        N = matrix.shape[0]
                        edges_mb = np.linspace(r_start_mb, r_end_mb, N + 1)
                        J, I = np.meshgrid(edges_mb, edges_mb)
                        X = (I + J) / 2.0
                        Y = (J - I) / 2.0
                        
                        row_idx, col_idx = np.indices(matrix.shape)
                        masked_matrix = np.ma.masked_where((row_idx > col_idx) | (matrix < 1), matrix)
                        
                        cmap = plt.get_cmap('YlOrRd').copy()
                        cmap.set_bad('white')
                        
                        im = ax_hic.pcolormesh(X, Y, masked_matrix, cmap=cmap, 
                                               norm=LogNorm(vmin=1, vmax=hic_vmax), rasterized=True)
                        
                        ax_hic.set_ylim((r_end_mb - r_start_mb) / 2.0, 0) 
                        ax_hic.set_yticks([])
                        
                        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
                        cbaxes = inset_axes(ax_hic, width="2%", height="30%", loc='lower left')
                        plt.colorbar(im, cax=cbaxes, orientation='vertical')
                        cbaxes.tick_params(labelsize=6)
                except Exception as e:
                    logger.warning(f"Could not plot Hi-C for {title_text}: {e}")
                    ax_hic.text(0.5, 0.5, 'Hi-C Data Unavailable', ha='center', va='center', transform=ax_hic.transAxes)
                    ax_hic.axis('off')

            # --- Final Formatting for the bottom X axis ---
            last_ax = axes[-1]
            last_ax.spines['bottom'].set_visible(True)
            last_ax.tick_params(bottom=True, labelbottom=True)
            last_ax.set_xlabel("Genomic position(Mbp)", fontsize=11)
            last_ax.set_xlim(r_start_mb, r_end_mb) 
            last_ax.xaxis.set_major_formatter(FuncFormatter(m_fmt))
            
            plt.tight_layout()
            plt.subplots_adjust(hspace=0.08) 
            pdf.savefig(fig, dpi=300)
            plt.close()
            
    logger.info(f"Multi-track plotting completed. Saved to: {out_pdf}")


# ================= Refined Mode Visualization =================

def plot_cross_chrom_heatmap(sim_matrix, short_labels, family_labels, ax):
    """Draw an annotated cross-chromosome identity heatmap on the given axes.
    
    Args:
        sim_matrix:    n x n numpy array of percent identity values (0-100)
        short_labels:  list of short display labels for axes
        family_labels: list of full labels (used to determine chromosome grouping)
        ax:            matplotlib Axes to draw on
    """
    n = len(short_labels)
    
    # Custom colormap: cream/flesh (0) -> warm orange (50) -> deep red (100)
    cmap = LinearSegmentedColormap.from_list(
        'flesh_to_red',
        ['#FFF5EB', '#FDDCB5', '#FCAE6B', '#F16B13', '#CC2900', '#8B0000'],
        N=256
    )
    
    im = ax.imshow(sim_matrix, cmap=cmap, vmin=0, vmax=100, aspect='equal')
    
    # Annotate cells with values
    for i in range(n):
        for j in range(n):
            val = sim_matrix[i][j]
            # Light background (low values) -> black text; dark background (high) -> white text
            text_color = 'white' if val > 65 else '#333333'
            fontsize = 7 if n <= 15 else (5 if n <= 40 else 3.5)
            ax.text(j, i, f"{val:.0f}", ha='center', va='center', 
                    fontsize=fontsize, color=text_color, fontweight='bold' if i == j else 'normal')
    
    # Axis labels - fully vertical X labels
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    tick_fontsize = 8 if n <= 15 else (6 if n <= 40 else 4.5)
    ax.set_xticklabels(short_labels, rotation=90, ha='center', va='top', fontsize=tick_fontsize)
    ax.set_yticklabels(short_labels, fontsize=tick_fontsize)
    
    # Draw chromosome group separators
    chrom_groups = []
    for label in family_labels:
        chrom = label.split('_Family_')[0] if '_Family_' in label else label
        chrom_groups.append(chrom)
    
    boundaries = []
    for i in range(1, n):
        if chrom_groups[i] != chrom_groups[i - 1]:
            boundaries.append(i)
    
    for b in boundaries:
        ax.axhline(y=b - 0.5, color='#333333', linewidth=1.2)
        ax.axvline(x=b - 0.5, color='#333333', linewidth=1.2)
    
    # Colorbar
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('Sequence Identity (%)', fontsize=10)
    cbar.ax.tick_params(labelsize=8)
    
    ax.set_title("Cross-Chromosome Monomer Family Identity", fontsize=13, fontweight='bold', pad=15)


def plot_summary_table(sim_matrix, short_labels, family_labels, pdf, logger):
    """Draw a clean summary table page showing per-chromosome dominant families
    and notable cross-chromosome relationships.
    
    This provides an at-a-glance overview complementing the detailed heatmap.
    """
    n = len(short_labels)
    if n == 0:
        return
    
    # --- Part 1: Per-chromosome family summary ---
    # Group families by chromosome
    chrom_families = {}  # { chrom: [(short_label, size, rank)] }
    for i, label in enumerate(family_labels):
        if '_Family_' in label:
            chrom = label.split('_Family_')[0]
            rest = label.split('_Family_')[1]
            parts = rest.split('_')
            rank = parts[0]
            size = ""
            for k, p in enumerate(parts):
                if p == "Size" and k + 1 < len(parts):
                    size = parts[k + 1]
        else:
            chrom = label
            rank = "?"
            size = ""
        
        if chrom not in chrom_families:
            chrom_families[chrom] = []
        chrom_families[chrom].append((short_labels[i], size, rank, i))
    
    sorted_chroms = sorted(chrom_families.keys(), key=natural_sort_key)
    
    # --- Part 2: Cross-chromosome same-family pairs ---
    cross_pairs = []
    for i in range(n):
        chrom_i = family_labels[i].split('_Family_')[0] if '_Family_' in family_labels[i] else ""
        for j in range(i + 1, n):
            chrom_j = family_labels[j].split('_Family_')[0] if '_Family_' in family_labels[j] else ""
            if chrom_i != chrom_j and sim_matrix[i][j] >= 60:
                relationship = "Same family" if sim_matrix[i][j] >= 80 else "Divergent"
                cross_pairs.append((short_labels[i], short_labels[j], sim_matrix[i][j], relationship))
    
    cross_pairs.sort(key=lambda x: -x[2])
    
    # --- Draw table ---
    fig, axes = plt.subplots(2, 1, figsize=(11, max(8, len(sorted_chroms) * 0.45 + 4)),
                             gridspec_kw={'height_ratios': [max(1, len(sorted_chroms)), max(1, len(cross_pairs))]})
    
    # Table 1: Per-chromosome summary
    ax1 = axes[0]
    ax1.axis('off')
    ax1.set_title("Per-Chromosome Monomer Family Summary", fontsize=13, fontweight='bold', pad=10)
    
    table1_data = []
    for chrom in sorted_chroms:
        fams = chrom_families[chrom]
        fam_count = len(fams)
        dominant = fams[0]  # Rank 1 is always first
        fam_names = ", ".join(f.split('(')[0].split('_')[-1] for f, _, _, _ in fams)
        sizes = ", ".join(s if s else "?" for _, s, _, _ in fams)
        table1_data.append([chrom, str(fam_count), fam_names, sizes])
    
    if table1_data:
        t1 = ax1.table(cellText=table1_data,
                       colLabels=["Chromosome", "Families", "Family Ranks", "Cluster Sizes"],
                       cellLoc='center', loc='center', colColours=['#4472C4'] * 4)
        t1.auto_set_font_size(False)
        t1.set_fontsize(8 if len(sorted_chroms) <= 20 else 6.5)
        t1.scale(1, 1.3)
        
        # Style header
        for key, cell in t1.get_celld().items():
            row, col = key
            if row == 0:
                cell.set_text_props(color='white', fontweight='bold')
                cell.set_facecolor('#4472C4')
            else:
                cell.set_facecolor('#D9E2F3' if row % 2 == 0 else 'white')
            cell.set_edgecolor('#B4C7E7')
    
    # Table 2: Cross-chromosome relationships
    ax2 = axes[1]
    ax2.axis('off')
    ax2.set_title("Cross-Chromosome Relationships (Identity ≥ 60%)", fontsize=13, fontweight='bold', pad=10)
    
    if cross_pairs:
        table2_data = [[a, b, f"{ident:.1f}%", rel] for a, b, ident, rel in cross_pairs[:30]]  # Cap at 30 rows
        
        t2 = ax2.table(cellText=table2_data,
                       colLabels=["Family A", "Family B", "Identity", "Relationship"],
                       cellLoc='center', loc='center', colColours=['#4472C4'] * 4)
        t2.auto_set_font_size(False)
        t2.set_fontsize(8 if len(cross_pairs) <= 20 else 6.5)
        t2.scale(1, 1.3)
        
        for key, cell in t2.get_celld().items():
            row, col = key
            if row == 0:
                cell.set_text_props(color='white', fontweight='bold')
                cell.set_facecolor('#4472C4')
            else:
                cell.set_facecolor('#D9E2F3' if row % 2 == 0 else 'white')
                # Color-code relationship column
                if col == 3:
                    rel_text = cell.get_text().get_text()
                    if rel_text == "Same family":
                        cell.set_facecolor('#C6EFCE')
                    elif rel_text == "Divergent":
                        cell.set_facecolor('#FFEB9C')
            cell.set_edgecolor('#B4C7E7')
    else:
        ax2.text(0.5, 0.5, "No cross-chromosome pairs with ≥ 60% identity found.\n"
                 "All chromosomes appear to use distinct monomer families.",
                 ha='center', va='center', fontsize=11, style='italic', color='#666666',
                 transform=ax2.transAxes)
    
    plt.tight_layout()
    pdf.savefig(fig, dpi=300)
    plt.close()
    logger.info("Summary table page added to report.")


def plot_refined_report(per_chrom_density_tsvs, genome_sizes_file, sim_matrix, short_labels, 
                        family_labels, out_pdf, window_size, logger, fig_width=10.0, track_height=1.8):
    """
    Generate the refined mode PDF report:
      - One page per chromosome showing its own monomer families' density
      - One final summary page with the cross-chromosome identity heatmap
    
    Args:
        per_chrom_density_tsvs: dict { chrom: { family_display_name: tsv_path } }
        genome_sizes_file:      path to genome.sizes file
        sim_matrix:             n x n identity matrix from cross_chrom_analysis
        short_labels:           short display labels for heatmap
        family_labels:          full labels for heatmap grouping
        out_pdf:                output PDF path
        window_size:            window size used in density calculation
        logger:                 logger instance
    """
    logger.info("--- Step 5: Refined Mode Visualization ---")
    
    # Parse chromosome sizes
    chrom_sizes = {}
    with open(genome_sizes_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            chrom_sizes[parts[0]] = int(parts[1])
    
    monomer_colors = ['#f8766d', '#00bfc4', '#7cae00', '#c77cff', '#ff7f00', '#ffff33']
    
    with PdfPages(out_pdf) as pdf:
        
        # === Per-chromosome density pages ===
        for chrom in sorted(per_chrom_density_tsvs.keys(), key=natural_sort_key):
            fam_tsvs = per_chrom_density_tsvs[chrom]
            if not fam_tsvs:
                continue
                
            c_len_bp = chrom_sizes.get(chrom)
            if c_len_bp is None:
                logger.warning(f"Chromosome {chrom} not in genome.sizes, skipping plot.")
                continue
            
            c_len_mb = c_len_bp / 1_000_000.0
            num_tracks = len(fam_tsvs)
            
            logger.info(f"Plotting {chrom}: {num_tracks} monomer families ...")
            
            fig, axes = plt.subplots(num_tracks, 1, figsize=(fig_width, track_height * num_tracks), 
                                     sharex=True)
            if num_tracks == 1:
                axes = [axes]
            
            fig.suptitle(f"{chrom}  (Refined Mode)", fontsize=13, fontweight='bold', y=0.98)
            
            for j, (fam_name, tsv_path) in enumerate(fam_tsvs.items()):
                ax = axes[j]
                df = pd.read_csv(tsv_path, sep='\t', names=['chrom', 'start', 'end', 'count'])
                sub_df = df[df['chrom'] == chrom]
                
                color = monomer_colors[j % len(monomer_colors)]
                
                # Clean display name: "chr1_Family_1_size100" -> "Family_1"
                display_name = fam_name
                if '_Family_' in fam_name:
                    rest = fam_name.split('_Family_', 1)[1]
                    rank = rest.split('_size')[0]
                    size_part = rest.split('_size')[1] if '_size' in rest else ''
                    display_name = f"Family {rank} (n={size_part})" if size_part else f"Family {rank}"
                
                max_y = 0
                if not sub_df.empty and sub_df['count'].sum() > 0:
                    x_mb = sub_df['start'] / 1_000_000.0
                    ax.fill_between(x_mb, sub_df['count'], color=color, alpha=0.4)
                    ax.plot(x_mb, sub_df['count'], color=color, lw=1.2)
                    max_y = sub_df['count'].max()
                
                style_axis(ax, display_name)
                y_upper = max_y * 1.1 if max_y > 0 else 10
                ax.set_ylim(-y_upper * 0.05, y_upper)
            
            # Format bottom axis
            last_ax = axes[-1]
            last_ax.spines['bottom'].set_visible(True)
            last_ax.tick_params(bottom=True, labelbottom=True)
            last_ax.set_xlabel("Genomic position (Mbp)", fontsize=11)
            last_ax.set_xlim(0, c_len_mb)
            last_ax.xaxis.set_major_formatter(FuncFormatter(m_fmt))
            
            plt.tight_layout(rect=[0, 0, 1, 0.96])
            plt.subplots_adjust(hspace=0.08)
            pdf.savefig(fig, dpi=300)
            plt.close()
        
        # === Summary table page ===
        if sim_matrix is not None and sim_matrix.size > 1:
            plot_summary_table(sim_matrix, short_labels, family_labels, pdf, logger)
        
        # === Summary heatmap page ===
        if sim_matrix is not None and sim_matrix.size > 1:
            n = len(short_labels)
            # Dynamic figure size: generous for large matrices
            if n <= 15:
                heatmap_size = 8
            elif n <= 40:
                heatmap_size = max(10, n * 0.4 + 4)
            else:
                heatmap_size = max(16, n * 0.3 + 5)
            
            fig, ax = plt.subplots(1, 1, figsize=(heatmap_size, heatmap_size * 0.9))
            plot_cross_chrom_heatmap(sim_matrix, short_labels, family_labels, ax)
            
            # Extra bottom margin for vertical X labels
            plt.subplots_adjust(bottom=0.18, left=0.15)
            pdf.savefig(fig, dpi=300)
            plt.close()
            
            logger.info("Cross-chromosome identity heatmap page added to report.")
    
    logger.info(f"Refined mode report saved to: {out_pdf}")
