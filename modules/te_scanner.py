import sys
from collections import defaultdict

def get_genome_size(fai_path):
    total_size = 0
    try:
        with open(fai_path, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    total_size += int(parts[1])
    except FileNotFoundError:
        print(f"Error: Genome index file not found at {fai_path}")
        sys.exit(1)
    return total_size

def scan_te_out(out_path):
    stats = defaultdict(lambda: {"count": 0, "length": 0})
    try:
        with open(out_path, 'r') as f:
            for _ in range(3): next(f) # Skip header
            for line in f:
                parts = line.strip().split()
                if len(parts) < 11: continue
                try:
                    start, end = int(parts[5]), int(parts[6])
                except ValueError: continue
                
                raw_class = parts[10]
                length = end - start + 1
                stats[raw_class]["count"] += 1
                stats[raw_class]["length"] += length
    except FileNotFoundError:
        print(f"Error: RepeatMasker .out file not found at {out_path}")
        sys.exit(1)
    return stats

def run_scan_te(out_path, fai_path):
    """Scan entry point called by the AutoCen.py main program."""
    genome_size = get_genome_size(fai_path)
    stats = scan_te_out(out_path)
    
    print("\n" + "="*75)
    print(f"{'AutoCen TE Distribution Summary':^75}")
    print("="*75)
    print(f"{'TE Class/Family (Raw)':<30} | {'Count':<10} | {'Total Length (bp)':<18} | {'Genome %':<10}")
    print("-" * 75)
    
    sorted_stats = sorted(stats.items(), key=lambda x: x[1]['length'], reverse=True)
    for te_name, data in sorted_stats:
        count = data['count']
        length = data['length']
        percentage = (length / genome_size) * 100
        if percentage < 0.01: continue
        print(f"{te_name:<30} | {count:<10,} | {length:<18,} | {percentage:>6.2f}%")
        
    print("="*75)
    print("Tip: Use these exact names for AutoCen's --te_mapping argument.\n")
