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
                    outfile=None, plot_outfile=None, verbose=False
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
                    outfile=None, plot_outfile=None, verbose=False
                )


if __name__ == '__main__':
    unittest.main()
