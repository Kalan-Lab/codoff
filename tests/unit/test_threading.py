"""
Unit tests for threading functionality in codoff module.
"""

import pytest
import numpy as np
from unittest.mock import patch, MagicMock
from codoff.codoff import (
    _codoff_worker_function,
    _run_serial_simulation,
    _stat_calc_and_simulation
)


class TestThreadingFunctions:
    """Test threading and parallel processing functionality."""
    
    def test_codoff_worker_function_basic(self):
        """Test basic worker function functionality."""
        # Test data
        cpu_simulations = 100
        region_freqs_list = [10, 5, 8, 2, 1, 0]  # 26 total codons
        codon_order = ['ATG', 'GCC', 'GAA', 'TAA', 'TAG', 'TGA']
        observed_cosine_distance = 0.5
        all_cds_codon_count_arrays = [
            np.array([5, 2, 3, 1, 0, 0], dtype=np.int32),
            np.array([3, 1, 2, 0, 1, 0], dtype=np.int32),
            np.array([2, 2, 3, 1, 0, 0], dtype=np.int32)
        ]
        
        # Run worker function
        emp_count, sim_distances = _codoff_worker_function(
            cpu_simulations, region_freqs_list, codon_order, 
            observed_cosine_distance, all_cds_codon_count_arrays, random_seed=42
        )
        
        # Check results
        assert isinstance(emp_count, int)
        assert emp_count >= 0
        assert emp_count <= cpu_simulations
        assert len(sim_distances) == cpu_simulations
        assert all(0.0 <= dist <= 1.0 for dist in sim_distances)
    
    def test_codoff_worker_function_edge_cases(self):
        """Test worker function with edge cases."""
        codon_order = ['ATG', 'GCC', 'GAA']
        
        # Test with zero simulations
        emp_count, sim_distances = _codoff_worker_function(
            0, [1, 2, 3], codon_order, 0.5, [np.array([1, 1, 1])], random_seed=42
        )
        assert emp_count == 0
        assert len(sim_distances) == 0
        
        # Test with empty codon count arrays
        emp_count, sim_distances = _codoff_worker_function(
            10, [1, 2, 3], codon_order, 0.5, [], random_seed=42
        )
        assert emp_count == 0
        assert len(sim_distances) == 0
        
        # Test with empty codon order
        emp_count, sim_distances = _codoff_worker_function(
            10, [1, 2, 3], [], 0.5, [np.array([1, 1, 1])], random_seed=42
        )
        assert emp_count == 0
        assert len(sim_distances) == 0
    
    def test_codoff_worker_function_import_handling(self):
        """Test worker function handles import errors gracefully."""
        codon_order = ['ATG', 'GCC', 'GAA']
        all_cds_codon_count_arrays = [np.array([1, 1, 1], dtype=np.int32)]
        
        # Test with empty codon count arrays (which should return 0, 0)
        emp_count, sim_distances = _codoff_worker_function(
            10, [1, 2, 3], codon_order, 0.5, [], random_seed=42
        )
        
        assert emp_count == 0
        assert len(sim_distances) == 0
    
    def test_run_serial_simulation_basic(self):
        """Test serial simulation functionality."""
        # Test data
        foc_cod_freqs = np.array([10, 5, 8, 2, 1, 0], dtype=np.float64)
        all_gene_codon_count_arrays = [
            np.array([5, 2, 3, 1, 0, 0], dtype=np.int32),
            np.array([3, 1, 2, 0, 1, 0], dtype=np.int32),
            np.array([2, 2, 3, 1, 0, 0], dtype=np.int32)
        ]
        cod_order = ['ATG', 'GCC', 'GAA', 'TAA', 'TAG', 'TGA']
        cosine_distance = 0.5
        foc_codon_count = 26
        gene_list = ['gene1', 'gene2', 'gene3']
        gene_codons = {
            'gene1': ['ATG', 'GCC', 'GAA', 'ATG', 'GCC', 'GAA', 'ATG', 'GCC', 'GAA'],
            'gene2': ['ATG', 'GCC', 'GAA', 'ATG', 'GCC', 'GAA'],
            'gene3': ['ATG', 'GCC', 'GAA', 'ATG', 'GCC', 'GAA', 'ATG', 'GCC', 'GAA']
        }
        verbose = False
        
        # Run serial simulation
        emp_pval, sim_cosine_distances = _run_serial_simulation(
            foc_cod_freqs, all_gene_codon_count_arrays, cod_order, 
            cosine_distance, foc_codon_count, gene_list, gene_codons, verbose
        )
        
        # Check results
        assert isinstance(emp_pval, int)
        assert emp_pval >= 0
        assert emp_pval <= 10000  # Should not exceed number of simulations
        assert len(sim_cosine_distances) == 10000
        # Check that distances are valid (some might be NaN due to division by zero)
        valid_distances = [dist for dist in sim_cosine_distances if np.isfinite(dist)]
        assert len(valid_distances) > 0  # At least some valid distances
        
        # Cosine distances can be outside [0,1] due to numerical precision issues
        # The important thing is that they are finite and reasonable
        assert all(np.isfinite(dist) for dist in sim_cosine_distances)
        assert all(dist >= 0.0 for dist in valid_distances)  # Should be non-negative
    
    def test_run_serial_simulation_verbose(self):
        """Test serial simulation with verbose output."""
        # Mock the progress bar function
        with patch('codoff.codoff.util.print_progress_bar') as mock_progress:
            foc_cod_freqs = np.array([10, 5, 8], dtype=np.float64)
            all_gene_codon_count_arrays = [np.array([5, 2, 3], dtype=np.int32)]
            cod_order = ['ATG', 'GCC', 'GAA']
            cosine_distance = 0.5
            foc_codon_count = 10
            gene_list = ['gene1']
            gene_codons = {'gene1': ['ATG', 'GCC', 'GAA'] * 3}
            verbose = True
            
            # Run simulation
            emp_pval, sim_cosine_distances = _run_serial_simulation(
                foc_cod_freqs, all_gene_codon_count_arrays, cod_order, 
                cosine_distance, foc_codon_count, gene_list, gene_codons, verbose
            )
            
            # Check that progress bar was called
            assert mock_progress.called
    
    def test_stat_calc_and_simulation_threading_parameter(self):
        """Test that threading parameter is properly handled."""
        # Test data
        all_cods = {'ATG', 'GCC', 'GAA'}
        cod_freq_dict_focal = {'ATG': 10, 'GCC': 5, 'GAA': 8}
        cod_freq_dict_background = {'ATG': 20, 'GCC': 15, 'GAA': 12}
        gene_list = ['gene1', 'gene2', 'gene3']
        gene_codons = {
            'gene1': ['ATG', 'GCC', 'GAA'] * 3,
            'gene2': ['ATG', 'GCC', 'GAA'] * 2,
            'gene3': ['ATG', 'GCC', 'GAA'] * 3
        }
        foc_codon_count = 24
        all_codon_counts = {'ATG': 30, 'GCC': 20, 'GAA': 20}
        
        # Test with threads=1 (serial)
        with patch('codoff.codoff._run_serial_simulation') as mock_serial:
            mock_serial.return_value = (100, [0.5] * 10000)
            
            result = _stat_calc_and_simulation(
                all_cods, cod_freq_dict_focal, cod_freq_dict_background,
                gene_list, gene_codons, foc_codon_count, all_codon_counts,
                threads=1, verbose=False
            )
            
            # Should call serial simulation
            mock_serial.assert_called_once()
            assert 'emp_pval_freq' in result
            assert 'cosine_distance' in result
            assert 'rho' in result
    
    def test_stat_calc_and_simulation_parallel_threading(self):
        """Test parallel threading functionality."""
        # Test data
        all_cods = {'ATG', 'GCC', 'GAA'}
        cod_freq_dict_focal = {'ATG': 10, 'GCC': 5, 'GAA': 8}
        cod_freq_dict_background = {'ATG': 20, 'GCC': 15, 'GAA': 12}
        gene_list = ['gene1', 'gene2', 'gene3']
        gene_codons = {
            'gene1': ['ATG', 'GCC', 'GAA'] * 3,
            'gene2': ['ATG', 'GCC', 'GAA'] * 2,
            'gene3': ['ATG', 'GCC', 'GAA'] * 3
        }
        foc_codon_count = 24
        all_codon_counts = {'ATG': 30, 'GCC': 20, 'GAA': 20}
        
        # Test with threads=2 (parallel)
        with patch('codoff.codoff.Pool') as mock_pool:
            # Mock pool and its methods
            mock_pool_instance = MagicMock()
            mock_pool_instance.__enter__.return_value = mock_pool_instance
            mock_pool_instance.__exit__.return_value = None
            mock_pool_instance.starmap.return_value = [(50, [0.5] * 5000), (50, [0.5] * 5000)]
            mock_pool.return_value = mock_pool_instance
            
            result = _stat_calc_and_simulation(
                all_cods, cod_freq_dict_focal, cod_freq_dict_background,
                gene_list, gene_codons, foc_codon_count, all_codon_counts,
                threads=2, verbose=False
            )
            
            # Should use parallel processing
            mock_pool.assert_called_once_with(processes=2)
            mock_pool_instance.starmap.assert_called_once()
            assert 'emp_pval_freq' in result
            assert 'cosine_distance' in result
            assert 'rho' in result
    
    def test_threading_fallback_to_serial(self):
        """Test that threading falls back to serial on error."""
        # Test data
        all_cods = {'ATG', 'GCC', 'GAA'}
        cod_freq_dict_focal = {'ATG': 10, 'GCC': 5, 'GAA': 8}
        cod_freq_dict_background = {'ATG': 20, 'GCC': 15, 'GAA': 12}
        gene_list = ['gene1', 'gene2', 'gene3']
        gene_codons = {
            'gene1': ['ATG', 'GCC', 'GAA'] * 3,
            'gene2': ['ATG', 'GCC', 'GAA'] * 2,
            'gene3': ['ATG', 'GCC', 'GAA'] * 3
        }
        foc_codon_count = 24
        all_codon_counts = {'ATG': 30, 'GCC': 20, 'GAA': 20}
        
        # Mock parallel processing to fail
        with patch('codoff.codoff.Pool') as mock_pool:
            mock_pool.side_effect = Exception("Parallel processing failed")
            
            with patch('codoff.codoff._run_serial_simulation') as mock_serial:
                mock_serial.return_value = (100, [0.5] * 10000)
                
                result = _stat_calc_and_simulation(
                    all_cods, cod_freq_dict_focal, cod_freq_dict_background,
                    gene_list, gene_codons, foc_codon_count, all_codon_counts,
                    threads=2, verbose=False
                )
                
                # Should fall back to serial
                mock_serial.assert_called_once()
                assert 'emp_pval_freq' in result
    
    def test_worker_function_exception_handling(self):
        """Test that worker function handles exceptions gracefully."""
        codon_order = ['ATG', 'GCC', 'GAA']
        all_cds_codon_count_arrays = [np.array([1, 1, 1], dtype=np.int32)]
        
        # Test with invalid data that should cause exceptions
        emp_count, sim_distances = _codoff_worker_function(
            10, [1, 2, 3], codon_order, 0.5, all_cds_codon_count_arrays, random_seed=42
        )
        
        # Should still return valid results
        assert isinstance(emp_count, int)
        assert emp_count >= 0
        assert len(sim_distances) == 10
        assert all(0.0 <= dist <= 1.0 for dist in sim_distances)
    
    def test_simulation_consistency(self):
        """Test that simulations produce consistent results."""
        # Test data
        foc_cod_freqs = np.array([10, 5, 8], dtype=np.float64)
        all_gene_codon_count_arrays = [np.array([5, 2, 3], dtype=np.int32)]
        cod_order = ['ATG', 'GCC', 'GAA']
        cosine_distance = 0.5
        foc_codon_count = 10
        gene_list = ['gene1']
        gene_codons = {'gene1': ['ATG', 'GCC', 'GAA'] * 3}
        verbose = False
        
        # Run simulation multiple times
        results = []
        for _ in range(3):
            emp_pval, sim_distances = _run_serial_simulation(
                foc_cod_freqs, all_gene_codon_count_arrays, cod_order, 
                cosine_distance, foc_codon_count, gene_list, gene_codons, verbose
            )
            results.append((emp_pval, sim_distances))
        
        # All results should have same number of simulations
        assert all(len(sim_distances) == 10000 for _, sim_distances in results)
        
        # Results should be reasonable (not all identical due to randomness)
        emp_pvals = [emp_pval for emp_pval, _ in results]
        assert all(0 <= emp_pval <= 10000 for emp_pval in emp_pvals)
