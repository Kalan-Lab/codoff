#!/usr/bin/env python3

import unittest
import sys
import os

# Add the src directory to the path so we can import codoff
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from codoff import codoff


class TestCachingConsistency(unittest.TestCase):
    """Test that cached and non-cached paths produce identical results for codoff_main_gbk"""
    
    @classmethod
    def setUpClass(cls):
        """Set up test data paths"""
        cls.test_dir = os.path.join(os.path.dirname(__file__), '..', 'Csimulans_Data', 'Coryne_simulans_PES1')
        cls.full_genome_gbk = os.path.join(cls.test_dir, 'Coryne_simulans_PES1.gbk')
        cls.focal_gbk = os.path.join(cls.test_dir, 'NZ_CP014634.1.region001.gbk')
        
        # Check if test files exist
        if not os.path.exists(cls.full_genome_gbk):
            raise unittest.SkipTest(f"Test genome file not found: {cls.full_genome_gbk}")
        if not os.path.exists(cls.focal_gbk):
            raise unittest.SkipTest(f"Test focal file not found: {cls.focal_gbk}")
    
    def test_cached_vs_noncached_identical_results(self):
        """Test that using genome_data cache produces identical results to non-cached"""
        
        # Run WITHOUT caching (normal path)
        result_nocache = codoff.codoff_main_gbk(
            full_genome_file=self.full_genome_gbk,
            focal_genbank_files=[self.focal_gbk],
            verbose=False,
            num_sims=1000,
            seed=42
        )
        
        # Extract genome data once
        genome_data = codoff.extract_genome_codon_data(
            full_genome_file=self.full_genome_gbk,
            verbose=False
        )
        
        # Run WITH caching
        result_cached = codoff.codoff_main_gbk(
            full_genome_file=self.full_genome_gbk,
            focal_genbank_files=[self.focal_gbk],
            genome_data=genome_data,
            verbose=False,
            num_sims=1000,
            seed=42
        )
        
        # Compare results - they should be identical
        self.assertEqual(result_nocache['empirical_freq'], result_cached['empirical_freq'],
                        "Empirical frequencies should be identical")
        self.assertEqual(result_nocache['cosine_distance'], result_cached['cosine_distance'],
                        "Cosine distances should be identical")
        self.assertEqual(result_nocache['rho'], result_cached['rho'],
                        "Rho values should be identical")
        self.assertEqual(result_nocache['codon_order'], result_cached['codon_order'],
                        "Codon orders should be identical")
        self.assertEqual(result_nocache['focal_region_codons'], result_cached['focal_region_codons'],
                        "Focal region codons should be identical")
        self.assertEqual(result_nocache['background_genome_codons'], result_cached['background_genome_codons'],
                        "Background genome codons should be identical")
    
    def test_extract_genome_data_contains_required_keys(self):
        """Test that extract_genome_codon_data returns all required keys"""
        genome_data = codoff.extract_genome_codon_data(
            full_genome_file=self.full_genome_gbk,
            verbose=False
        )
        
        required_keys = [
            'locus_tag_sequences', 'gene_codons', 'gene_list', 'all_cods',
            'all_codon_counts', 'total_cds_length', 'gene_coords', 'scaffold_lengths'
        ]
        
        for key in required_keys:
            self.assertIn(key, genome_data, f"genome_data should contain '{key}'")


if __name__ == '__main__':
    unittest.main()

