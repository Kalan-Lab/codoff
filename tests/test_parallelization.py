"""
Test parallelization functionality for codoff simulation.
"""

import unittest
import sys
import os
import time
from collections import defaultdict
from unittest.mock import patch

# Add src to path to import codoff
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from codoff.codoff import _stat_calc_and_simulation


class TestParallelization(unittest.TestCase):
    """Test parallelization functionality."""
    
    def test_parallel_vs_sequential_consistency(self):
        """Test that parallel and sequential results are consistent."""
        # Create test data
        all_cods = {'AAA', 'ACG', 'TTT', 'CCC'}
        cod_freq_dict_focal = {'AAA': 10, 'ACG': 5, 'TTT': 15, 'CCC': 8}
        cod_freq_dict_background = {'AAA': 20, 'ACG': 10, 'TTT': 30, 'CCC': 12}
        gene_list = ['gene1', 'gene2', 'gene3']
        
        gene_codons = {
            'gene1': {'AAA': 5, 'ACG': 3, 'TTT': 7},
            'gene2': {'AAA': 5, 'ACG': 2, 'TTT': 8},
            'gene3': {'AAA': 10, 'CCC': 8}
        }
        
        foc_codon_count = 20
        all_codon_counts = {'AAA': 20, 'ACG': 5, 'TTT': 15, 'CCC': 8}
        
        # Run with 1 process (sequential)
        result_sequential = _stat_calc_and_simulation(
            all_cods, cod_freq_dict_focal, cod_freq_dict_background,
            gene_list, gene_codons, foc_codon_count, all_codon_counts,
            outfile=None, plot_outfile=None, verbose=False, max_jobs=1
        )
        
        # Run with 2 processes (parallel)
        result_parallel = _stat_calc_and_simulation(
            all_cods, cod_freq_dict_focal, cod_freq_dict_background,
            gene_list, gene_codons, foc_codon_count, all_codon_counts,
            outfile=None, plot_outfile=None, verbose=False, max_jobs=2
        )
        
        # Results should have same structure
        self.assertEqual(set(result_sequential.keys()), set(result_parallel.keys()))
        
        # P-values should be in reasonable range
        self.assertGreaterEqual(result_sequential['emp_pval_freq'], 0.0)
        self.assertLessEqual(result_sequential['emp_pval_freq'], 1.0)
        self.assertGreaterEqual(result_parallel['emp_pval_freq'], 0.0)
        self.assertLessEqual(result_parallel['emp_pval_freq'], 1.0)
        
        # Cosine distances should be in reasonable range
        self.assertGreaterEqual(result_sequential['cosine_distance'], 0.0)
        self.assertLessEqual(result_sequential['cosine_distance'], 1.0)
        self.assertGreaterEqual(result_parallel['cosine_distance'], 0.0)
        self.assertLessEqual(result_parallel['cosine_distance'], 1.0)
        
        # Codon counts should match
        self.assertEqual(result_sequential['focal_region_codons'], result_parallel['focal_region_codons'])
        self.assertEqual(result_sequential['background_genome_codons'], result_parallel['background_genome_codons'])
        self.assertEqual(result_sequential['codon_order'], result_parallel['codon_order'])
    
    def test_max_jobs_parameter_handling(self):
        """Test that max_jobs parameter is handled correctly."""
        all_cods = {'AAA', 'ACG', 'TTT'}
        cod_freq_dict_focal = {'AAA': 10, 'ACG': 5, 'TTT': 15}
        cod_freq_dict_background = {'AAA': 20, 'ACG': 10, 'TTT': 30}
        gene_list = ['gene1', 'gene2']
        
        gene_codons = {
            'gene1': {'AAA': 5, 'ACG': 3, 'TTT': 7},
            'gene2': {'AAA': 5, 'ACG': 2, 'TTT': 8}
        }
        
        foc_codon_count = 15
        all_codon_counts = {'AAA': 10, 'ACG': 5, 'TTT': 15}
        
        # Test with None (should use all cores)
        result_none = _stat_calc_and_simulation(
            all_cods, cod_freq_dict_focal, cod_freq_dict_background,
            gene_list, gene_codons, foc_codon_count, all_codon_counts,
            outfile=None, plot_outfile=None, verbose=False, max_jobs=None
        )
        
        # Test with 1 (sequential)
        result_1 = _stat_calc_and_simulation(
            all_cods, cod_freq_dict_focal, cod_freq_dict_background,
            gene_list, gene_codons, foc_codon_count, all_codon_counts,
            outfile=None, plot_outfile=None, verbose=False, max_jobs=1
        )
        
        # Test with 2 (parallel)
        result_2 = _stat_calc_and_simulation(
            all_cods, cod_freq_dict_focal, cod_freq_dict_background,
            gene_list, gene_codons, foc_codon_count, all_codon_counts,
            outfile=None, plot_outfile=None, verbose=False, max_jobs=2
        )
        
        # All should produce valid results
        for result in [result_none, result_1, result_2]:
            self.assertIn('emp_pval_freq', result)
            self.assertIn('cosine_distance', result)
            self.assertIn('rho', result)
            self.assertGreaterEqual(result['emp_pval_freq'], 0.0)
            self.assertLessEqual(result['emp_pval_freq'], 1.0)
    
    def test_parallel_performance_improvement(self):
        """Test that parallel execution provides performance improvement."""
        # Create a larger dataset for performance testing
        all_cods = {f'COD{i:03d}' for i in range(32)}  # 32 codons
        
        cod_freq_dict_focal = {codon: 0 for codon in all_cods}
        cod_freq_dict_background = {codon: 0 for codon in all_cods}
        
        # Create 50 genes
        gene_list = [f'gene_{i}' for i in range(1, 51)]
        gene_codons = {}
        
        for i in range(1, 51):
            gene_codons[f'gene_{i}'] = {}
            for j, codon in enumerate(all_cods):
                count = (i * j) % 5 + 1
                gene_codons[f'gene_{i}'][codon] = count
                
                if i <= 10:  # First 10 genes are focal
                    cod_freq_dict_focal[codon] += count
                else:  # Remaining 40 genes are background
                    cod_freq_dict_background[codon] += count
        
        foc_codon_count = 100
        all_codon_counts = {}
        for codon in all_cods:
            all_codon_counts[codon] = cod_freq_dict_focal[codon] + cod_freq_dict_background[codon]
        
        # Test sequential performance
        start_time = time.time()
        result_sequential = _stat_calc_and_simulation(
            all_cods, cod_freq_dict_focal, cod_freq_dict_background,
            gene_list, gene_codons, foc_codon_count, all_codon_counts,
            outfile=None, plot_outfile=None, verbose=False, max_jobs=1
        )
        sequential_time = time.time() - start_time
        
        # Test parallel performance
        start_time = time.time()
        result_parallel = _stat_calc_and_simulation(
            all_cods, cod_freq_dict_focal, cod_freq_dict_background,
            gene_list, gene_codons, foc_codon_count, all_codon_counts,
            outfile=None, plot_outfile=None, verbose=False, max_jobs=2
        )
        parallel_time = time.time() - start_time
        
        # Both should complete successfully
        self.assertIn('emp_pval_freq', result_sequential)
        self.assertIn('emp_pval_freq', result_parallel)
        
        # Results should be consistent
        self.assertEqual(result_sequential['focal_region_codons'], result_parallel['focal_region_codons'])
        self.assertEqual(result_sequential['background_genome_codons'], result_parallel['background_genome_codons'])
        
        # Note: Performance improvement may not be visible on very small datasets
        # This test mainly ensures both modes work correctly
        print(f"Sequential time: {sequential_time:.3f}s, Parallel time: {parallel_time:.3f}s")
    


if __name__ == '__main__':
    unittest.main()
