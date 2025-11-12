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


if __name__ == '__main__':
    unittest.main()
