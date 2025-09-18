"""
Integration tests for the complete codoff workflow.
"""

import unittest
import sys
import os
import tempfile
from unittest.mock import patch, MagicMock
from contextlib import redirect_stderr
from io import StringIO
from collections import defaultdict

# Add src to path to import codoff
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from codoff.codoff import _stat_calc_and_simulation


class TestIntegrationWorkflow(unittest.TestCase):
    """Test the complete workflow with realistic data."""
    
    def test_realistic_genome_simulation(self):
        """Test simulation with realistic genome data."""
        # Create realistic codon data
        all_cods = {
            'AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT',
            'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT',
            'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT',
            'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT',
            'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT',
            'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT',
            'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT',
            'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT'
        }
        
        # Create focal region with specific codon bias
        cod_freq_dict_focal = {codon: 0 for codon in all_cods}
        cod_freq_dict_focal['AAA'] = 50  # High AAA usage
        cod_freq_dict_focal['TTT'] = 30  # High TTT usage
        cod_freq_dict_focal['CCC'] = 20  # Moderate CCC usage
        
        # Create background genome with different codon bias
        cod_freq_dict_background = {codon: 0 for codon in all_cods}
        cod_freq_dict_background['AAA'] = 100  # Lower AAA usage
        cod_freq_dict_background['TTT'] = 200  # Lower TTT usage
        cod_freq_dict_background['CCC'] = 150  # Higher CCC usage
        
        # Create gene list with multiple genes
        gene_list = [f'gene_{i}' for i in range(1, 21)]  # 20 genes
        
        # Create gene_codons with realistic distributions
        gene_codons = defaultdict(lambda: defaultdict(int))
        
        # Focal genes (first 5) - biased towards AAA and TTT
        for i in range(1, 6):
            gene_codons[f'gene_{i}']['AAA'] = 10
            gene_codons[f'gene_{i}']['TTT'] = 8
            gene_codons[f'gene_{i}']['CCC'] = 2
            gene_codons[f'gene_{i}']['ACG'] = 5
            gene_codons[f'gene_{i}']['GCA'] = 3
        
        # Background genes (remaining 15) - biased towards CCC
        for i in range(6, 21):
            gene_codons[f'gene_{i}']['AAA'] = 5
            gene_codons[f'gene_{i}']['TTT'] = 4
            gene_codons[f'gene_{i}']['CCC'] = 8
            gene_codons[f'gene_{i}']['ACG'] = 3
            gene_codons[f'gene_{i}']['GCA'] = 2
        
        foc_codon_count = 100  # 100 codons from focal region
        
        # Calculate total codon counts
        all_codon_counts = {}
        for codon in all_cods:
            all_codon_counts[codon] = cod_freq_dict_focal[codon] + cod_freq_dict_background[codon]
        
        # Run simulation
        result = _stat_calc_and_simulation(
            all_cods, cod_freq_dict_focal, cod_freq_dict_background,
            gene_list, gene_codons, foc_codon_count, all_codon_counts,
            outfile=None, plot_outfile=None, verbose=False
        )
        
        # Verify results
        self.assertIn('emp_pval_freq', result)
        self.assertIn('cosine_distance', result)
        self.assertIn('rho', result)
        self.assertIn('codon_order', result)
        self.assertIn('focal_region_codons', result)
        self.assertIn('background_genome_codons', result)
        
        # Verify data types
        self.assertIsInstance(result['emp_pval_freq'], float)
        self.assertIsInstance(result['cosine_distance'], float)
        self.assertIsInstance(result['rho'], float)
        self.assertIsInstance(result['codon_order'], list)
        self.assertIsInstance(result['focal_region_codons'], list)
        self.assertIsInstance(result['background_genome_codons'], list)
        
        # Verify reasonable values
        self.assertGreaterEqual(result['emp_pval_freq'], 0.0)
        self.assertLessEqual(result['emp_pval_freq'], 1.0)
        self.assertGreaterEqual(result['cosine_distance'], 0.0)
        self.assertLessEqual(result['cosine_distance'], 1.0)
        self.assertGreaterEqual(result['rho'], -1.0)
        self.assertLessEqual(result['rho'], 1.0)
        
        # Verify codon counts (should be original observed values)
        total_focal = sum(result['focal_region_codons'])
        total_background = sum(result['background_genome_codons'])
        expected_focal = sum(cod_freq_dict_focal.values())
        expected_background = sum(cod_freq_dict_background.values())
        self.assertEqual(total_focal, expected_focal)
        self.assertEqual(total_background, expected_background)
        
        # Verify codon order
        self.assertEqual(len(result['codon_order']), len(all_cods))
        self.assertEqual(set(result['codon_order']), all_cods)
    
    def test_performance_with_large_dataset(self):
        """Test performance with a large dataset."""
        # Create a large dataset
        all_cods = {f'COD{i:03d}' for i in range(64)}  # 64 codons
        
        cod_freq_dict_focal = {codon: 0 for codon in all_cods}
        cod_freq_dict_background = {codon: 0 for codon in all_cods}
        
        # Create 100 genes
        gene_list = [f'gene_{i}' for i in range(1, 101)]
        gene_codons = defaultdict(lambda: defaultdict(int))
        
        for i in range(1, 101):
            for j, codon in enumerate(all_cods):
                # Each gene has different codon usage - create variation
                count = (i * j) % 10 + 1
                gene_codons[f'gene_{i}'][codon] = count
                
                # Update frequency dictionaries
                if i <= 20:  # First 20 genes are focal
                    cod_freq_dict_focal[codon] += count
                else:  # Remaining 80 genes are background
                    cod_freq_dict_background[codon] += count
        
        foc_codon_count = 500  # 500 codons from focal region
        all_codon_counts = {}
        for codon in all_cods:
            all_codon_counts[codon] = cod_freq_dict_focal[codon] + cod_freq_dict_background[codon]
        
        # Run simulation and measure time
        import time
        start_time = time.time()
        
        result = _stat_calc_and_simulation(
            all_cods, cod_freq_dict_focal, cod_freq_dict_background,
            gene_list, gene_codons, foc_codon_count, all_codon_counts,
            outfile=None, plot_outfile=None, verbose=False
        )
        
        end_time = time.time()
        execution_time = end_time - start_time
        
        # Verify results
        self.assertIn('emp_pval_freq', result)
        self.assertIn('cosine_distance', result)
        self.assertIn('rho', result)
        
        # Verify performance (should complete in reasonable time)
        self.assertLess(execution_time, 10.0)  # Should complete within 10 seconds
        
        # Verify codon counts (should be original observed values)
        total_focal = sum(result['focal_region_codons'])
        total_background = sum(result['background_genome_codons'])
        expected_focal = sum(cod_freq_dict_focal.values())
        expected_background = sum(cod_freq_dict_background.values())
        self.assertEqual(total_focal, expected_focal)
        self.assertEqual(total_background, expected_background)
    
    def test_edge_case_empty_focal_region(self):
        """Test behavior when focal region has no codons."""
        all_cods = {'AAA', 'ACG', 'TTT'}
        cod_freq_dict_focal = {'AAA': 0, 'ACG': 0, 'TTT': 0}
        cod_freq_dict_background = {'AAA': 10, 'ACG': 5, 'TTT': 15}
        gene_list = ['gene1']
        gene_codons = defaultdict(lambda: defaultdict(int))
        foc_codon_count = 0
        all_codon_counts = {'AAA': 10, 'ACG': 5, 'TTT': 15}
        
        # This should handle empty focal region gracefully without printing to screen
        with redirect_stderr(StringIO()):
            with self.assertRaises(SystemExit):
                _stat_calc_and_simulation(
                    all_cods, cod_freq_dict_focal, cod_freq_dict_background,
                    gene_list, gene_codons, foc_codon_count, all_codon_counts,
                    outfile=None, plot_outfile=None, verbose=False
                )
    
    def test_edge_case_single_gene(self):
        """Test behavior with only one gene."""
        all_cods = {'AAA', 'ACG', 'TTT'}
        cod_freq_dict_focal = {'AAA': 5, 'ACG': 3, 'TTT': 7}
        cod_freq_dict_background = {'AAA': 10, 'ACG': 6, 'TTT': 14}
        gene_list = ['gene1', 'gene2']
        gene_codons = defaultdict(lambda: defaultdict(int))
        gene_codons['gene1']['AAA'] = 5
        gene_codons['gene1']['ACG'] = 3
        gene_codons['gene1']['TTT'] = 7
        gene_codons['gene2']['AAA'] = 10
        gene_codons['gene2']['ACG'] = 6
        gene_codons['gene2']['TTT'] = 14
        foc_codon_count = 15
        all_codon_counts = {'AAA': 15, 'ACG': 9, 'TTT': 21}
        
        result = _stat_calc_and_simulation(
            all_cods, cod_freq_dict_focal, cod_freq_dict_background,
            gene_list, gene_codons, foc_codon_count, all_codon_counts,
            outfile=None, plot_outfile=None, verbose=False
        )
        
        # Verify results
        self.assertIn('emp_pval_freq', result)
        self.assertIn('cosine_distance', result)
        self.assertIn('rho', result)
        
        # Verify codon counts (should be original observed values)
        total_focal = sum(result['focal_region_codons'])
        total_background = sum(result['background_genome_codons'])
        expected_focal = sum(cod_freq_dict_focal.values())
        expected_background = sum(cod_freq_dict_background.values())
        self.assertEqual(total_focal, expected_focal)
        self.assertEqual(total_background, expected_background)


if __name__ == '__main__':
    unittest.main()
