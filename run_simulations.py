#!/usr/bin/env python3
"""
Script to run codoff simulations across random regions of a genome
and generate histograms of the resulting discordance percentiles.
"""

import argparse
import random
import sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from tqdm import tqdm

from src.codoff import util
from src.codoff.codoff import codoff_main_coords


def run_one_simulation(args_tuple):
    """Helper function for parallel execution - must be at module level for pickling."""
    genome_file, scaffold_id, start, end, num_sims, seed = args_tuple
    try:
        result = codoff_main_coords(
            full_genome_file=str(genome_file),
            focal_scaffold=scaffold_id,
            focal_start_coord=start,
            focal_end_coord=end,
            verbose=False,
            num_sims=num_sims,
            seed=seed
        )
        return result
    except SystemExit as e:
        return None
    except Exception as e:
        return None


def main():
    parser = argparse.ArgumentParser(
        description="Run codoff simulations on random genome regions and plot discordance percentile distributions."
    )
    parser.add_argument("-g", "--genome-file", type=Path, required=True, help="Path to genome file (GenBank or FASTA)")
    parser.add_argument("-ns","--num-sims", type=int, default=1000, help="Number of random genomic regions to test [Default: 1000]")
    parser.add_argument("-rs","--region-size", type=int, default=100000, help="Size of each focal region in bp [Default: 100000]")
    parser.add_argument("-o","--output-file", type=Path, default="discordance_percentile_histogram.png", help="Output file for histogram [Default: discordance_percentile_histogram.png]")
    parser.add_argument("-d","--data-file", type=Path, default=None, help="Optional: Save raw percentile data to CSV file")
    parser.add_argument("-ncs","--num-codoff-sims", type=int, default=1000, help="Number of simulations for each codoff run [Default: 1000]")
    parser.add_argument("-mj", "--max-jobs", type=int, default=1, help="Number of parallel jobs [Default: 1]")
    parser.add_argument('-x', '--seed', type=int, default=None, help='Random seed for reproducible region selection')
    
    args = parser.parse_args()
    
    # Set the master seed for the simulation initiator if provided
    if args.seed is not None:
        random.seed(args.seed)
        
    # Parse genome file
    sys.stderr.write(f"Reading genome file: {args.genome_file}\n")
    is_genbank = util.checkIsGenBankWithCDS(str(args.genome_file))
    
    if is_genbank:
        records = list(SeqIO.parse(str(args.genome_file), "genbank"))
    else:
        records = list(SeqIO.parse(str(args.genome_file), "fasta"))
    
    if not records:
        sys.stderr.write("ERROR: No sequences found in genome file\n")
        sys.exit(1)
    
    # Filter scaffolds large enough for the region size
    valid_scaffolds = [(rec.id, len(rec.seq)) for rec in records if len(rec.seq) > args.region_size]
    
    if not valid_scaffolds:
        sys.stderr.write(f"ERROR: No scaffolds larger than {args.region_size} bp\n")
        sys.exit(1)
    
    sys.stderr.write(f"Found {len(valid_scaffolds)} valid scaffolds\n")
    sys.stderr.write(f"Running {args.num_sims} simulations with {args.max_jobs} parallel jobs...\n")
    
    # Prepare simulation tasks (now include all args needed by the worker function)
    tasks = []
    for i in range(args.num_sims):
        scaffold_id, scaffold_len = random.choice(valid_scaffolds)
        start = random.randint(1, scaffold_len - args.region_size)
        end = start + args.region_size
        
        # Generate a unique seed for each simulation run
        # This makes each run reproducible, but ensures they are different
        run_seed = random.randint(0, 2**32 - 1)
        tasks.append((args.genome_file, scaffold_id, start, end, args.num_codoff_sims, run_seed))
    
    # Run simulations (with or without parallelization)
    percentiles = []
    
    if args.max_jobs == 1:
        # Sequential execution with simple progress
        for i, task in enumerate(tasks):
            if (i + 1) % 100 == 0:
                sys.stderr.write(f"  Completed {i + 1}/{args.num_sims} simulations\n")
            result = run_one_simulation(task)
            if result:
                empirical_freq = result.get('empirical_freq')
                if empirical_freq is not None:
                    percentile = empirical_freq * 100.0
                    percentiles.append(percentile)
    else:
        # Parallel execution with progress bar
        with ProcessPoolExecutor(max_workers=args.max_jobs) as executor:
            futures = [executor.submit(run_one_simulation, task) for task in tasks]
            for future in tqdm(as_completed(futures), total=len(futures), desc="Simulations"):
                result = future.result()
                if result:
                    empirical_freq = result.get('empirical_freq')
                    if empirical_freq is not None:
                        percentile = empirical_freq * 100.0
                        percentiles.append(percentile)
    
    sys.stderr.write(f"Completed all simulations\n")
    sys.stderr.write(f"  Discordance percentiles collected: {len(percentiles)}\n")
    
    # Generate plot
    if percentiles:
        plt.figure(figsize=(10, 6))
        sns.histplot(percentiles, bins=50, kde=True)
        plt.xlabel("Discordance Percentile (%)")
        plt.ylabel("Frequency")
        plt.title(f"Discordance Percentile Distribution\n({len(percentiles)} random genomic regions of {args.region_size:,} bp)")
        plt.xlim(0, 100)
        plt.axvline(x=5, color='r', linestyle='--', alpha=0.5, label='5% threshold')
        plt.axvline(x=10, color='orange', linestyle='--', alpha=0.5, label='10% threshold')
        plt.legend()
        plt.savefig(args.output_file, dpi=300, bbox_inches='tight')
        plt.close()
        sys.stderr.write(f"Saved histogram to {args.output_file}\n")
        
        # Print summary statistics
        sys.stderr.write(f"\nSummary Statistics:\n")
        sys.stderr.write(f"  Mean: {np.mean(percentiles):.2f}%\n")
        sys.stderr.write(f"  Median: {np.median(percentiles):.2f}%\n")
        sys.stderr.write(f"  Std Dev: {np.std(percentiles):.2f}%\n")
        sys.stderr.write(f"  Min: {np.min(percentiles):.2f}%\n")
        sys.stderr.write(f"  Max: {np.max(percentiles):.2f}%\n")
        sys.stderr.write(f"  < 5%: {sum(1 for p in percentiles if p < 5)} ({100*sum(1 for p in percentiles if p < 5)/len(percentiles):.1f}%)\n")
        sys.stderr.write(f"  < 10%: {sum(1 for p in percentiles if p < 10)} ({100*sum(1 for p in percentiles if p < 10)/len(percentiles):.1f}%)\n")
        
        # Save raw data if requested
        if args.data_file:
            with open(args.data_file, 'w') as f:
                f.write("discordance_percentile\n")
                for p in percentiles:
                    f.write(f"{p:.4f}\n")
            sys.stderr.write(f"Saved raw data to {args.data_file}\n")
    else:
        sys.stderr.write("No discordance percentiles collected\n")
    
    sys.stderr.write("\nDone!\n")


if __name__ == "__main__":
    main()
