"""
Unit tests for codon analysis functionality in codoff module.
"""

import pytest
import numpy as np
from scipy import stats, spatial
from codoff.codoff import check_data_type


class TestCodonAnalysis:
    """Test codon analysis functionality."""
    
    def test_check_data_type_valid(self):
        """Test data type checking with valid inputs."""
        # Test string
        assert check_data_type("test_string", str) == True
        assert check_data_type("", str) == True
        
        # Test integer
        assert check_data_type(123, int) == True
        assert check_data_type(0, int) == True
        
        # Test list
        assert check_data_type([1, 2, 3], list) == True
        assert check_data_type([], list) == True
        
        # Test boolean
        assert check_data_type(True, bool) == True
        assert check_data_type(False, bool) == True
    
    def test_check_data_type_invalid(self):
        """Test data type checking with invalid inputs."""
        # Test wrong types
        assert check_data_type(123, str) == False
        assert check_data_type("test", int) == False
        assert check_data_type([1, 2, 3], str) == False
        # Note: In Python, bool is a subclass of int, so True is considered an int
        # assert check_data_type(True, int) == False  # This would fail
        
        # Test None
        assert check_data_type(None, str) == False
        assert check_data_type(None, int) == False
    
    def test_cosine_distance_calculation(self):
        """Test cosine distance calculation between codon frequencies."""
        # Test identical vectors
        vec1 = np.array([1, 2, 3, 4, 5], dtype=np.float64)
        vec2 = np.array([1, 2, 3, 4, 5], dtype=np.float64)
        distance = spatial.distance.cosine(vec1, vec2)
        assert distance == 0.0
        
        # Test orthogonal vectors
        vec1 = np.array([1, 0, 0], dtype=np.float64)
        vec2 = np.array([0, 1, 0], dtype=np.float64)
        distance = spatial.distance.cosine(vec1, vec2)
        assert abs(distance - 1.0) < 1e-10
        
        # Test different vectors
        vec1 = np.array([1, 2, 3], dtype=np.float64)
        vec2 = np.array([2, 4, 6], dtype=np.float64)
        distance = spatial.distance.cosine(vec1, vec2)
        assert distance == 0.0  # Should be identical (proportional)
    
    def test_spearman_correlation_calculation(self):
        """Test Spearman correlation calculation."""
        # Test perfect positive correlation
        x = np.array([1, 2, 3, 4, 5])
        y = np.array([2, 4, 6, 8, 10])
        rho, p_value = stats.spearmanr(x, y)
        assert abs(rho - 1.0) < 1e-10
        
        # Test perfect negative correlation
        x = np.array([1, 2, 3, 4, 5])
        y = np.array([5, 4, 3, 2, 1])
        rho, p_value = stats.spearmanr(x, y)
        assert abs(rho - (-1.0)) < 1e-10
        
        # Test no correlation (constant array)
        x = np.array([1, 2, 3, 4, 5])
        y = np.array([1, 1, 1, 1, 1])
        with pytest.warns(stats.ConstantInputWarning):
            rho, p_value = stats.spearmanr(x, y)
        # Should be NaN for constant input
        assert np.isnan(rho)
    
    def test_codon_frequency_arrays(self):
        """Test codon frequency array operations."""
        # Test array creation and operations
        codon_order = ['ATG', 'GCC', 'GAA', 'TAA', 'TAG', 'TGA']
        focal_freqs = np.array([10, 5, 8, 2, 1, 0], dtype=np.float64)
        background_freqs = np.array([20, 15, 12, 5, 3, 2], dtype=np.float64)
        
        # Test array properties
        assert len(focal_freqs) == len(codon_order)
        assert len(background_freqs) == len(codon_order)
        assert focal_freqs.dtype == np.float64
        assert background_freqs.dtype == np.float64
        
        # Test array operations
        total_freqs = focal_freqs + background_freqs
        expected_total = np.array([30, 20, 20, 7, 4, 2])
        assert np.array_equal(total_freqs, expected_total)
        
        # Test background calculation
        background_calc = total_freqs - focal_freqs
        assert np.array_equal(background_calc, background_freqs)
    
    def test_codon_frequency_validation(self):
        """Test codon frequency validation logic."""
        # Test arrays with variation
        focal_freqs = np.array([10, 5, 8, 2, 1, 0], dtype=np.float64)
        background_freqs = np.array([20, 15, 12, 5, 3, 2], dtype=np.float64)
        
        # Should have variation (more than 1 unique value)
        assert len(set(focal_freqs)) > 1
        assert len(set(background_freqs)) > 1
        
        # Test arrays without variation (should fail validation)
        focal_no_var = np.array([5, 5, 5, 5, 5, 5], dtype=np.float64)
        background_no_var = np.array([10, 10, 10, 10, 10, 10], dtype=np.float64)
        
        assert len(set(focal_no_var)) == 1
        assert len(set(background_no_var)) == 1
    
    def test_numpy_array_operations(self):
        """Test numpy array operations used in codon analysis."""
        # Test array creation
        codon_counts = np.array([10, 5, 8, 2, 1, 0], dtype=np.int32)
        assert codon_counts.dtype == np.int32
        
        # Test array summation
        total_codons = np.sum(codon_counts)
        assert total_codons == 26
        
        # Test array indexing
        assert codon_counts[0] == 10
        assert codon_counts[-1] == 0
        
        # Test array slicing
        first_three = codon_counts[:3]
        assert np.array_equal(first_three, np.array([10, 5, 8]))
        
        # Test array broadcasting
        proportions = codon_counts / total_codons
        expected_proportions = np.array([10/26, 5/26, 8/26, 2/26, 1/26, 0/26])
        assert np.allclose(proportions, expected_proportions)
    
    def test_codon_order_consistency(self):
        """Test codon order consistency in analysis."""
        codon_order = ['ATG', 'GCC', 'GAA', 'TAA', 'TAG', 'TGA']
        
        # Test that codon order is sorted (this specific order is not sorted)
        # assert codon_order == sorted(codon_order)  # This would fail
        
        # Test that all codons are 3 characters
        for codon in codon_order:
            assert len(codon) == 3
            assert all(base in 'ATCG' for base in codon)
        
        # Test that codon order is unique
        assert len(codon_order) == len(set(codon_order))
        
        # Test that we can sort the codon order
        sorted_order = sorted(codon_order)
        assert sorted_order == ['ATG', 'GAA', 'GCC', 'TAA', 'TAG', 'TGA']
    
    def test_edge_cases(self):
        """Test edge cases in codon analysis."""
        # Test empty arrays
        empty_array = np.array([], dtype=np.float64)
        assert len(empty_array) == 0
        
        # Test single element arrays
        single_array = np.array([5], dtype=np.float64)
        assert len(single_array) == 1
        assert np.sum(single_array) == 5
        
        # Test arrays with zeros
        zero_array = np.array([0, 0, 0, 0], dtype=np.float64)
        assert np.sum(zero_array) == 0
        
        # Test arrays with very small numbers
        small_array = np.array([1e-10, 2e-10, 3e-10], dtype=np.float64)
        assert np.sum(small_array) > 0
        assert np.all(small_array > 0)
    
    def test_statistical_calculations(self):
        """Test statistical calculations used in codoff."""
        # Test with realistic codon frequency data
        focal_freqs = np.array([100, 50, 80, 20, 10, 5], dtype=np.float64)
        background_freqs = np.array([200, 150, 120, 50, 30, 20], dtype=np.float64)
        
        # Test Spearman correlation
        rho, p_value = stats.spearmanr(focal_freqs, background_freqs)
        assert -1.0 <= rho <= 1.0
        assert 0.0 <= p_value <= 1.0
        
        # Test cosine distance
        cosine_dist = spatial.distance.cosine(focal_freqs, background_freqs)
        assert 0.0 <= cosine_dist <= 1.0
        
        # Test that both calculations work with the same data
        assert not np.isnan(rho)
        assert not np.isnan(cosine_dist)
        assert np.isfinite(rho)
        assert np.isfinite(cosine_dist)
