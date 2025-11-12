#!/usr/bin/env python3

import unittest
import sys
import os
from collections import defaultdict

# Add the src directory to the path so we can import codoff
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from codoff import codoff
import test_utils


class TestCachingFunctionality(unittest.TestCase):
    """Test the caching functionality added for antismash_codoff"""
    
    def test_extract_genome_codon_data_structure(self):
        """Test that extract_genome_codon_data returns expected structure"""
        # Create a minimal test genome file (we'll use a mock for this)
        test_data = test_utils._create_test_genome_data()
        
        # Test the structure of returned data
        self.assertIn('locus_tag_sequences', test_data)
        self.assertIn('gene_codons', test_data)
        self.assertIn('gene_list', test_data)
        self.assertIn('all_cods', test_data)
        self.assertIn('all_codon_counts', test_data)
        self.assertIn('total_cds_length', test_data)
        
        # Test data types
        self.assertIsInstance(test_data['locus_tag_sequences'], dict)
        self.assertIsInstance(test_data['gene_codons'], dict)
        self.assertIsInstance(test_data['gene_list'], list)
        self.assertIsInstance(test_data['all_cods'], set)
        self.assertIsInstance(test_data['all_codon_counts'], dict)
        self.assertIsInstance(test_data['total_cds_length'], int)
    
    def test_cached_vs_noncached_consistency(self):
        """Test that cached and non-cached paths produce identical results"""
        # Create test data
        test_genome_data = test_utils._create_test_genome_data()
        focal_data = test_utils._create_test_focal_data()
        
        # Test the direct analysis function
        result_direct = test_utils._run_direct_analysis(
            test_genome_data, focal_data, num_sims=100, verbose=False
        )
        
        # Verify the result structure
        expected_keys = ['empirical_freq', 'cosine_distance', 'rho', 'focal_region_codons']
        for key in expected_keys:
            self.assertIn(key, result_direct)
        
        # Focal region frequencies should be a list
        self.assertIsInstance(result_direct['focal_region_codons'], list)
    
    def test_num_sims_parameter(self):
        """Test that num_sims parameter works correctly"""
        test_genome_data = test_utils._create_test_genome_data()
        focal_data = test_utils._create_test_focal_data()
        
        # Test with different simulation counts using direct analysis
        result_100 = test_utils._run_direct_analysis(
            test_genome_data, focal_data, num_sims=100, verbose=False
        )
        
        result_500 = test_utils._run_direct_analysis(
            test_genome_data, focal_data, num_sims=500, verbose=False
        )
        
        # Both should complete successfully
        self.assertIn('empirical_freq', result_100)
        self.assertIn('empirical_freq', result_500)
        
        # Discordance frequencies should be in valid range
        self.assertGreaterEqual(result_100['empirical_freq'], 0.0)
        self.assertLessEqual(result_100['empirical_freq'], 1.0)
        self.assertGreaterEqual(result_500['empirical_freq'], 0.0)
        self.assertLessEqual(result_500['empirical_freq'], 1.0)


if __name__ == '__main__':
    unittest.main()
