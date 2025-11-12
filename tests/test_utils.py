#!/usr/bin/env python3

"""Test utilities for codoff testing"""

import sys
import os
from collections import defaultdict

# Add the src directory to the path so we can import codoff
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))


def _create_test_genome_data():
    """Create test genome data for testing purposes"""
    return {
        'locus_tag_sequences': {
            'gene1': 'ATGAAACCCGGGTTT',  # 5 codons
            'gene2': 'ATGCCCAAATTT',     # 4 codons  
            'gene3': 'ATGGGGCCCTTT'      # 4 codons
        },
        'gene_codons': {
            'gene1': defaultdict(int, {'ATG': 1, 'AAA': 1, 'CCC': 1, 'GGG': 1, 'TTT': 1}),
            'gene2': defaultdict(int, {'ATG': 1, 'CCC': 1, 'AAA': 1, 'TTT': 1}),
            'gene3': defaultdict(int, {'ATG': 1, 'GGG': 1, 'CCC': 1, 'TTT': 1})
        },
        'gene_list': ['gene1', 'gene2', 'gene3'],
        'all_cods': {'ATG', 'AAA', 'CCC', 'GGG', 'TTT'},
        'all_codon_counts': defaultdict(int, {'ATG': 3, 'AAA': 2, 'CCC': 3, 'GGG': 2, 'TTT': 3}),
        'total_cds_length': 39
    }


def _create_test_focal_data():
    """Create test focal region data for testing"""
    return ['gene1']  # Focal region contains gene1


def _run_direct_analysis(genome_data, focal_data, num_sims=100, verbose=False):
    """Run direct analysis for comparison with cached version"""
    # Simulate the non-cached path for testing
    focal_lts = set(focal_data)
    
    cod_freq_dict_focal = defaultdict(int)
    cod_freq_dict_background = defaultdict(int)
    foc_codon_count = 0
    
    # Calculate focal and background counts
    for gene, codons in genome_data['gene_codons'].items():
        for codon, count in codons.items():
            if gene in focal_lts:
                cod_freq_dict_focal[codon] += count
                foc_codon_count += count
            else:
                cod_freq_dict_background[codon] += count
    
    # Return simplified result for testing
    return {
        'empirical_freq': 0.5,    # Dummy value for testing
        'cosine_distance': 0.1,   # Dummy value for testing
        'rho': 0.8,               # Dummy value for testing
        'focal_region_codons': [cod_freq_dict_focal[cod] for cod in sorted(genome_data['all_cods'])]
    }
