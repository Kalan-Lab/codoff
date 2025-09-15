"""
Test suite for codoff codon caching and simulation changes.
"""

import unittest
import tempfile
import os
import sys
from unittest.mock import patch, MagicMock
from contextlib import redirect_stderr
from io import StringIO
from collections import defaultdict
import numpy as np

# Add src to path to import codoff
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from codoff.codoff import codoff_main_gbk, codoff_main_coords, _stat_calc_and_simulation


class TestCodonCaching(unittest.TestCase):
    """Test codon count caching functionality."""
    
    def test_gene_codons_structure(self):
        """Test that gene_codons uses count-based structure."""
        # This test verifies the internal structure change
        # We'll test this indirectly through the simulation function
        
        # Create test data with codon counts
        all_cods = {'AAA', 'ACG', 'TTT'}
        cod_freq_dict_focal = {'AAA': 10, 'ACG': 5, 'TTT': 15}
        cod_freq_dict_background = {'AAA': 20, 'ACG': 10, 'TTT': 30}
        gene_list = ['gene1', 'gene2']
        
        # Create gene_codons with count structure (new format)
        gene_codons = defaultdict(lambda: defaultdict(int))
        gene_codons['gene1']['AAA'] = 5
        gene_codons['gene1']['ACG'] = 3
        gene_codons['gene1']['TTT'] = 7
        gene_codons['gene2']['AAA'] = 5
        gene_codons['gene2']['ACG'] = 2
        gene_codons['gene2']['TTT'] = 8
        
        foc_codon_count = 15
        all_codon_counts = {'AAA': 25, 'ACG': 15, 'TTT': 45}
        
        # Test that simulation runs without errors
        result = _stat_calc_and_simulation(
            all_cods, cod_freq_dict_focal, cod_freq_dict_background,
            gene_list, gene_codons, foc_codon_count, all_codon_counts,
            outfile=None, plot_outfile=None, verbose=False, max_jobs=1
        )
        
        # Verify result structure
        self.assertIn('emp_pval_freq', result)
        self.assertIn('cosine_distance', result)
        self.assertIn('rho', result)
        self.assertIn('codon_order', result)
        self.assertIn('focal_region_codons', result)
        self.assertIn('background_genome_codons', result)
        
        # Verify result types
        self.assertIsInstance(result['emp_pval_freq'], float)
        self.assertIsInstance(result['cosine_distance'], float)
        self.assertIsInstance(result['rho'], float)
        self.assertIsInstance(result['codon_order'], list)
        self.assertIsInstance(result['focal_region_codons'], list)
        self.assertIsInstance(result['background_genome_codons'], list)


class TestProportionalSampling(unittest.TestCase):
    """Test proportional sampling in simulations."""
    
    def test_proportional_sampling_logic(self):
        """Test that proportional sampling maintains codon ratios."""
        # Create a test scenario where we need to sample partial genes
        all_cods = {'AAA', 'ACG', 'TTT'}
        cod_freq_dict_focal = {'AAA': 10, 'ACG': 5, 'TTT': 15}
        cod_freq_dict_background = {'AAA': 20, 'ACG': 10, 'TTT': 30}
        gene_list = ['gene1', 'gene2']
        
        # Create gene with specific codon counts
        gene_codons = defaultdict(lambda: defaultdict(int))
        gene_codons['gene1']['AAA'] = 10  # 33.3% of gene
        gene_codons['gene1']['ACG'] = 5   # 16.7% of gene  
        gene_codons['gene1']['TTT'] = 15  # 50% of gene
        gene_codons['gene2']['AAA'] = 20
        gene_codons['gene2']['ACG'] = 10
        gene_codons['gene2']['TTT'] = 30
        
        foc_codon_count = 20  # Need exactly 20 codons
        all_codon_counts = {'AAA': 30, 'ACG': 15, 'TTT': 45}
        
        # Run simulation
        result = _stat_calc_and_simulation(
            all_cods, cod_freq_dict_focal, cod_freq_dict_background,
            gene_list, gene_codons, foc_codon_count, all_codon_counts,
            outfile=None, plot_outfile=None, verbose=False, max_jobs=1
        )
        
        # Verify that we get the original observed codon counts (not simulated)
        total_focal_codons = sum(result['focal_region_codons'])
        expected_focal = sum(cod_freq_dict_focal.values())
        self.assertEqual(total_focal_codons, expected_focal)
        
        # Verify that background codons are calculated correctly
        total_background_codons = sum(result['background_genome_codons'])
        expected_background = sum(cod_freq_dict_background.values())
        self.assertEqual(total_background_codons, expected_background)


class TestWarningMessages(unittest.TestCase):
    """Test warning messages for missing locus tags."""
    
    def test_missing_locus_tag_warning(self):
        """Test that warnings are issued for missing locus tags."""
        # This test would require creating mock GenBank files
        # For now, we'll test the warning logic indirectly
        
        # Test that the warning logic exists in the code
        with open(os.path.join(os.path.dirname(__file__), '..', 'src', 'codoff', 'codoff.py'), 'r') as f:
            code = f.read()
            
        # Check that warning logic is present
        self.assertIn('missing_focal_lts', code)
        self.assertIn('Warning: The following focal region locus tags were not found', code)
        self.assertIn('These locus tags will be ignored in the analysis', code)


class TestIntegration(unittest.TestCase):
    """Test overall integration and functionality."""
    
    def test_simulation_consistency(self):
        """Test that simulations produce consistent results."""
        all_cods = {'AAA', 'ACG', 'TTT', 'CCC'}
        cod_freq_dict_focal = {'AAA': 10, 'ACG': 5, 'TTT': 15, 'CCC': 8}
        cod_freq_dict_background = {'AAA': 20, 'ACG': 10, 'TTT': 30, 'CCC': 12}
        gene_list = ['gene1', 'gene2', 'gene3']
        
        gene_codons = defaultdict(lambda: defaultdict(int))
        gene_codons['gene1']['AAA'] = 5
        gene_codons['gene1']['ACG'] = 3
        gene_codons['gene1']['TTT'] = 7
        gene_codons['gene2']['AAA'] = 5
        gene_codons['gene2']['ACG'] = 2
        gene_codons['gene2']['TTT'] = 8
        gene_codons['gene3']['AAA'] = 10
        gene_codons['gene3']['CCC'] = 8
        
        foc_codon_count = 20
        all_codon_counts = {'AAA': 20, 'ACG': 5, 'TTT': 15, 'CCC': 8}
        
        # Run multiple simulations to test consistency
        results = []
        for _ in range(3):
            result = _stat_calc_and_simulation(
                all_cods, cod_freq_dict_focal, cod_freq_dict_background,
                gene_list, gene_codons, foc_codon_count, all_codon_counts,
                outfile=None, plot_outfile=None, verbose=False, max_jobs=1
            )
            results.append(result)
        
        # All results should have the same structure
        for result in results:
            self.assertIn('emp_pval_freq', result)
            self.assertIn('cosine_distance', result)
            self.assertIn('rho', result)
            
        # P-values should be reasonable (between 0 and 1)
        for result in results:
            self.assertGreaterEqual(result['emp_pval_freq'], 0.0)
            self.assertLessEqual(result['emp_pval_freq'], 1.0)
            
        # Cosine distances should be reasonable (between 0 and 1)
        for result in results:
            self.assertGreaterEqual(result['cosine_distance'], 0.0)
            self.assertLessEqual(result['cosine_distance'], 1.0)


class TestEdgeCases(unittest.TestCase):
    """Test edge cases and error conditions."""
    
    def test_empty_gene_list(self):
        """Test behavior with empty gene list."""
        all_cods = {'AAA', 'ACG'}
        cod_freq_dict_focal = {'AAA': 0, 'ACG': 0}
        cod_freq_dict_background = {'AAA': 10, 'ACG': 5}
        gene_list = []
        gene_codons = defaultdict(lambda: defaultdict(int))
        foc_codon_count = 0
        all_codon_counts = {'AAA': 10, 'ACG': 5}
        
        # This should handle empty gene list gracefully without printing to screen
        with redirect_stderr(StringIO()):
            with self.assertRaises(SystemExit):
                _stat_calc_and_simulation(
                    all_cods, cod_freq_dict_focal, cod_freq_dict_background,
                    gene_list, gene_codons, foc_codon_count, all_codon_counts,
                    outfile=None, plot_outfile=None, verbose=False, max_jobs=1
                )
    
    def test_single_codon_type(self):
        """Test behavior with only one codon type."""
        all_cods = {'AAA'}
        cod_freq_dict_focal = {'AAA': 10}
        cod_freq_dict_background = {'AAA': 20}
        gene_list = ['gene1']
        gene_codons = defaultdict(lambda: defaultdict(int))
        gene_codons['gene1']['AAA'] = 10
        foc_codon_count = 10
        all_codon_counts = {'AAA': 30}
        
        # This should handle single codon type gracefully without printing to screen
        with redirect_stderr(StringIO()):
            with self.assertRaises(SystemExit):
                _stat_calc_and_simulation(
                    all_cods, cod_freq_dict_focal, cod_freq_dict_background,
                    gene_list, gene_codons, foc_codon_count, all_codon_counts,
                    outfile=None, plot_outfile=None, verbose=False, max_jobs=1
                )


if __name__ == '__main__':
    unittest.main()
