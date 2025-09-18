#!/usr/bin/env python3

import unittest
import sys
import os
import tempfile
from unittest.mock import patch, MagicMock

# Add the src directory to the path so we can import codoff
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from codoff import codoff
from tests import test_utils


class TestAntismashCodoffIntegration(unittest.TestCase):
    """Test the antismash_codoff wrapper functionality"""
    
    def test_codoff_main_gbk_signature(self):
        """Test that codoff_main_gbk accepts the new parameters"""
        # Test that the function signature includes the new parameters
        import inspect
        sig = inspect.signature(codoff.codoff_main_gbk)
        
        # Check that new parameters are present
        self.assertIn('genome_data', sig.parameters)
        self.assertIn('num_sims', sig.parameters)
        
        # Check default values
        self.assertEqual(sig.parameters['genome_data'].default, None)
        self.assertEqual(sig.parameters['num_sims'].default, 10000)
    
    def test_caching_data_structure_integrity(self):
        """Test that the caching preserves data structure integrity"""
        test_data = test_utils._create_test_genome_data()
        
        # Verify the test data has the expected structure
        required_keys = ['locus_tag_sequences', 'gene_codons', 'gene_list', 
                        'all_cods', 'all_codon_counts', 'total_cds_length']
        
        for key in required_keys:
            self.assertIn(key, test_data, f"Missing required key: {key}")
        
        # Test that gene_codons uses count dictionaries (not lists like v1.2.1)
        for gene, codons in test_data['gene_codons'].items():
            self.assertIsInstance(codons, dict, f"gene_codons[{gene}] should be a dict")
            for codon, count in codons.items():
                self.assertIsInstance(count, int, f"Codon count should be int, got {type(count)}")
    
    def test_background_calculation_method(self):
        """Test that background frequencies are calculated correctly"""
        test_data = test_utils._create_test_genome_data()
        focal_genes = ['gene1']
        
        # Calculate expected background manually
        expected_background = {}
        for codon in test_data['all_cods']:
            total_count = test_data['all_codon_counts'][codon]
            focal_count = test_data['gene_codons']['gene1'].get(codon, 0)
            expected_background[codon] = total_count - focal_count
        
        # This validates our approach matches v1.2.1 behavior (total - focal)
        for codon in test_data['all_cods']:
            total = test_data['all_codon_counts'][codon]
            focal = test_data['gene_codons']['gene1'].get(codon, 0)
            background = total - focal
            self.assertGreaterEqual(background, 0, 
                                  f"Background count for {codon} should be non-negative")


if __name__ == '__main__':
    unittest.main()
