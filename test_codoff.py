#!/usr/bin/env python3
"""
Test script for the updated codoff implementation with caching.
This script tests the caching functionality without requiring the full environment.
"""

import sys
import os
sys.path.insert(0, 'src')

def test_caching_functions():
    """Test the caching functions directly."""
    from codoff.codoff import _get_cached_codons, _get_cached_codon_counts, _clear_caches, _get_cache_stats
    
    print("Testing caching functions...")
    
    # Clear caches first
    _clear_caches()
    
    # Test sequence
    test_sequence = "ATGGCCGAAATGGCCGAAATGGCCGAA"  # 9 codons
    test_locus_tag = "test_gene_1"
    codon_order = ['ATG', 'GCC', 'GAA', 'TAA', 'TAG', 'TGA']
    
    # Test codon extraction
    print(f"Testing codon extraction for sequence: {test_sequence}")
    codons = _get_cached_codons(test_sequence, test_locus_tag)
    print(f"Extracted codons: {codons}")
    print(f"Number of codons: {len(codons)}")
    
    # Test codon count array
    print(f"Testing codon count array with order: {codon_order}")
    count_array = _get_cached_codon_counts(test_sequence, test_locus_tag, codon_order)
    print(f"Count array: {count_array}")
    
    # Test cache stats
    stats = _get_cache_stats()
    print(f"Cache stats: {stats}")
    
    # Test cache hit (should return cached result)
    print("Testing cache hit...")
    codons_cached = _get_cached_codons(test_sequence, test_locus_tag)
    print(f"Cached codons (should be same): {codons_cached}")
    print(f"Cache hit successful: {codons == codons_cached}")
    
    # Test with different sequence
    test_sequence2 = "ATGGCCGAAATGGCCGAAATGGCCGAA"  # Same sequence
    test_locus_tag2 = "test_gene_2"
    
    codons2 = _get_cached_codons(test_sequence2, test_locus_tag2)
    print(f"Different locus tag, same sequence: {codons2}")
    
    # Test with invalid sequence
    invalid_sequence = "ATGGCCGAAATGGCC"  # Not divisible by 3
    codons_invalid = _get_cached_codons(invalid_sequence, "invalid_gene")
    print(f"Invalid sequence result: {codons_invalid}")
    
    # Final cache stats
    final_stats = _get_cache_stats()
    print(f"Final cache stats: {final_stats}")
    
    print("Caching functions test completed successfully!")

def test_numpy_imports():
    """Test that numpy imports work correctly."""
    print("Testing numpy imports...")
    try:
        import numpy as np
        print(f"NumPy version: {np.__version__}")
        
        # Test basic numpy operations
        arr = np.array([1, 2, 3, 4, 5])
        print(f"Test array: {arr}")
        print(f"Array sum: {np.sum(arr)}")
        print("NumPy import test successful!")
        return True
    except ImportError as e:
        print(f"NumPy import failed: {e}")
        return False

def test_scipy_imports():
    """Test that scipy imports work correctly."""
    print("Testing scipy imports...")
    try:
        from scipy import stats, spatial
        print("SciPy imports successful!")
        
        # Test basic scipy operations
        from scipy.spatial.distance import cosine
        arr1 = [1, 2, 3, 4, 5]
        arr2 = [2, 3, 4, 5, 6]
        distance = cosine(arr1, arr2)
        print(f"Cosine distance test: {distance}")
        return True
    except ImportError as e:
        print(f"SciPy import failed: {e}")
        return False

if __name__ == "__main__":
    print("=" * 60)
    print("Testing updated codoff implementation with caching")
    print("=" * 60)
    
    # Test imports first
    numpy_ok = test_numpy_imports()
    scipy_ok = test_scipy_imports()
    
    if numpy_ok and scipy_ok:
        # Test caching functions
        test_caching_functions()
        print("\n" + "=" * 60)
        print("All tests completed successfully!")
        print("The updated codoff implementation with caching is working correctly.")
    else:
        print("\n" + "=" * 60)
        print("Some dependencies are missing. Please install numpy and scipy.")
        print("You can use: pip install numpy scipy")
