#!/usr/bin/env python3
"""
Script to run codoff simulations across random regions of a genome
and generate histograms of the resulting p-values.
"""

import argparse
import random
import sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from tqdm import tqdm

from src.codoff import util
from src.codoff.codoff import codoff_main_coords


def run_one_simulation(args_tuple):
    """Helper function for parallel execution - must be at module level for pickling."""
    genome_file, scaffold_id, start, end, num_sims, sequential_sampling, seed = args_tuple
    try:
        result = codoff_main_coords(
            full_genome_file=str(genome_file),
            focal_scaffold=scaffold_id,
            focal_start_coord=start,
            focal_end_coord=end,
            verbose=False,
            num_sims=num_sims,
            sequential_sampling=sequential_sampling,
            seed=seed
        )
        return result
    except SystemExit as e:
        return None
    except Exception as e:
        return None


def main():
    parser = argparse.ArgumentParser(
        description="Run codoff simulations on random genome regions and plot p-value distributions."
    )
    parser.add_argument("-g", "--genome-file", type=Path, required=True, help="Path to genome file (GenBank or FASTA)")
    parser.add_argument("-ns","--num-sims", type=int, default=1000, help="Number of random genomic regions to test [Default: 1000]")
    parser.add_argument("-rs","--region-size", type=int, default=100000, help="Size of each focal region in bp [Default: 100000]")
    parser.add_argument("-o","--output-file", type=Path, default="pvalue_histogram.png", help="Output file for histogram [Default: pvalue_histogram.png]")
    parser.add_argument("-r", "--random-sampling", action="store_true", help="Use legacy random gene sampling instead of sequential contiguous-window sampling")
    parser.add_argument("-ncs","--num-codoff-sims", type=int, default=1000, help="Number of simulations for each codoff run [Default: 1000]")
    parser.add_argument("-mj", "--max-jobs", type=int, default=1, help="Number of parallel jobs [Default: 1]")
    parser.add_argument('-x', '--seed', type=int, default=None, help='Random seed for reproducible region selection')
    
    args = parser.parse_args()
    
    proportional_sampling = args.random_sampling
    sequential_sampling = not proportional_sampling
    
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
        tasks.append((args.genome_file, scaffold_id, start, end, args.num_codoff_sims, sequential_sampling, run_seed))
    
    # Run simulations (with or without parallelization)
    p_values = []
    
    if args.max_jobs == 1:
        # Sequential execution with simple progress
        for i, task in enumerate(tasks):
            if (i + 1) % 100 == 0:
                sys.stderr.write(f"  Completed {i + 1}/{args.num_sims} simulations\n")
            result = run_one_simulation(task)
            if result:
                p_val = result.get('emp_pval_freq')
                if p_val is not None:
                    p_values.append(p_val)
    else:
        # Parallel execution with progress bar
        with ProcessPoolExecutor(max_workers=args.max_jobs) as executor:
            futures = [executor.submit(run_one_simulation, task) for task in tasks]
            for future in tqdm(as_completed(futures), total=len(futures), desc="Simulations"):
                result = future.result()
                if result:
                    p_val = result.get('emp_pval_freq')
                    if p_val is not None:
                        p_values.append(p_val)
    
    sys.stderr.write(f"Completed all simulations\n")
    sys.stderr.write(f"  P-values collected: {len(p_values)}\n")
    
    # Generate plot
    if p_values:
        plt.figure(figsize=(8, 6))
        sns.histplot(p_values, bins=50, kde=True)
        plt.xlabel("Empirical P-value")
        plt.ylabel("Frequency")
        plt.title(f"P-value Distribution ({len(p_values)} simulations)")
        plt.savefig(args.output_file, dpi=300, bbox_inches='tight')
        plt.close()
        sys.stderr.write(f"Saved histogram to {args.output_file}\n")
    else:
        sys.stderr.write("No p-values collected\n")
    
    sys.stderr.write("Done!\n")


if __name__ == "__main__":
    main()
