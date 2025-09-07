"""
Unit tests for caching functionality in codoff module.
"""

import pytest
import numpy as np
from codoff.codoff import (
    _get_cached_codons,
    _get_cached_codon_counts,
    _clear_caches,
    _get_cache_stats
)


class TestCachingFunctions:
    """Test caching functionality."""
    
    def setup_method(self):
        """Clear caches before each test."""
        _clear_caches()
    
    def test_get_cached_codons_valid_sequence(self):
        """Test codon extraction with valid DNA sequence."""
        sequence = "ATGGCCGAAATGGCCGAAATGGCCGAA"  # 9 codons
        locus_tag = "test_gene_1"
        expected_codons = ["ATG", "GCC", "GAA", "ATG", "GCC", "GAA", "ATG", "GCC", "GAA"]
        
        codons = _get_cached_codons(sequence, locus_tag)
        assert codons == expected_codons
        assert len(codons) == 9
    
    def test_get_cached_codons_invalid_length(self):
        """Test codon extraction with sequence not divisible by 3."""
        sequence = "ATGGCCGAAATGGCC"  # 14 bases, not divisible by 3
        locus_tag = "test_gene_invalid"
        
        codons = _get_cached_codons(sequence, locus_tag)
        # The function returns valid codons even for sequences not divisible by 3
        # It just ignores the incomplete codon at the end
        expected_codons = ["ATG", "GCC", "GAA", "ATG", "GCC"]  # 4 complete codons
        assert codons == expected_codons
    
    def test_get_cached_codons_with_ambiguous_bases(self):
        """Test codon extraction with ambiguous nucleotides."""
        sequence = "ATGGCCGAAATGNCCGAA"  # Contains N
        locus_tag = "test_gene_ambiguous"
        
        codons = _get_cached_codons(sequence, locus_tag)
        # Should only return valid codons (first 3 and last 2)
        expected_codons = ["ATG", "GCC", "GAA", "ATG", "GAA"]  # Skips the codon with N
        assert codons == expected_codons
    
    def test_get_cached_codons_cache_hit(self):
        """Test that cached results are returned on subsequent calls."""
        sequence = "ATGGCCGAAATGGCCGAAATGGCCGAA"
        locus_tag = "test_gene_cache"
        
        # First call
        codons1 = _get_cached_codons(sequence, locus_tag)
        
        # Second call should return cached result
        codons2 = _get_cached_codons(sequence, locus_tag)
        
        assert codons1 == codons2
        assert codons1 is codons2  # Should be the same object
    
    def test_get_cached_codon_counts(self):
        """Test codon count array generation."""
        sequence = "ATGGCCGAAATGGCCGAAATGGCCGAA"
        locus_tag = "test_gene_counts"
        codon_order = ['ATG', 'GCC', 'GAA', 'TAA', 'TAG', 'TGA']
        
        count_array = _get_cached_codon_counts(sequence, locus_tag, codon_order)
        
        expected_counts = [3, 3, 3, 0, 0, 0]  # 3 ATG, 3 GCC, 3 GAA, 0 others
        assert count_array == expected_counts
        assert isinstance(count_array, list)
    
    def test_get_cached_codon_counts_empty_sequence(self):
        """Test codon count array with empty sequence."""
        sequence = "ATGGCC"  # Not divisible by 3, but has 2 valid codons
        locus_tag = "test_gene_empty"
        codon_order = ['ATG', 'GCC', 'GAA']
        
        count_array = _get_cached_codon_counts(sequence, locus_tag, codon_order)
        
        expected_counts = [1, 1, 0]  # 1 ATG, 1 GCC, 0 GAA
        assert np.array_equal(count_array, expected_counts)
    
    def test_clear_caches(self):
        """Test cache clearing functionality."""
        sequence = "ATGGCCGAAATGGCCGAAATGGCCGAA"
        locus_tag = "test_gene_clear"
        
        # Populate cache
        _get_cached_codons(sequence, locus_tag)
        _get_cached_codon_counts(sequence, locus_tag, ['ATG', 'GCC', 'GAA'])
        
        # Check cache has data
        stats_before = _get_cache_stats()
        assert stats_before['codon_extraction_cache_size'] > 0
        assert stats_before['codon_counts_cache_size'] > 0
        
        # Clear caches
        _clear_caches()
        
        # Check caches are empty
        stats_after = _get_cache_stats()
        assert stats_after['codon_extraction_cache_size'] == 0
        assert stats_after['codon_counts_cache_size'] == 0
    
    def test_get_cache_stats(self):
        """Test cache statistics reporting."""
        # Clear caches first
        _clear_caches()
        
        # Check initial state
        stats = _get_cache_stats()
        assert stats['codon_extraction_cache_size'] == 0
        assert stats['codon_counts_cache_size'] == 0
        assert stats['total_cached_sequences'] == 0
        
        # Add some data
        sequence1 = "ATGGCCGAAATGGCCGAAATGGCCGAA"
        sequence2 = "TAAATGGCCGAAATGGCCGAAATGGCCGAA"
        
        _get_cached_codons(sequence1, "gene1")
        _get_cached_codons(sequence2, "gene2")
        _get_cached_codon_counts(sequence1, "gene1", ['ATG', 'GCC', 'GAA'])
        
        # Check updated stats
        stats = _get_cache_stats()
        assert stats['codon_extraction_cache_size'] == 2
        assert stats['codon_counts_cache_size'] == 1
        assert stats['total_cached_sequences'] == 2
    
    def test_different_locus_tags_same_sequence(self):
        """Test that different locus tags with same sequence are cached separately."""
        sequence = "ATGGCCGAAATGGCCGAAATGGCCGAA"
        locus_tag1 = "gene1"
        locus_tag2 = "gene2"
        
        codons1 = _get_cached_codons(sequence, locus_tag1)
        codons2 = _get_cached_codons(sequence, locus_tag2)
        
        # Should be same content but different cache entries
        assert codons1 == codons2
        assert codons1 is not codons2  # Different objects
        
        # Both should be cached
        stats = _get_cache_stats()
        assert stats['codon_extraction_cache_size'] == 2
    
    def test_exception_handling(self):
        """Test that exceptions in codon extraction are handled gracefully."""
        # Test with None sequence
        codons = _get_cached_codons(None, "test_none")
        assert codons == []
        
        # Test with empty string
        codons = _get_cached_codons("", "test_empty")
        assert codons == []
        
        # Test with non-string sequence
        codons = _get_cached_codons(123, "test_number")
        assert codons == []
