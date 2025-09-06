"""
Sample test data and fixtures for codoff testing.
"""

import numpy as np
from collections import defaultdict


class SampleData:
    """Sample data for testing codoff functionality."""
    
    # Sample DNA sequences
    SAMPLE_DNA_SEQUENCES = {
        'short_sequence': 'ATGGCCGAAATGGCCGAAATGGCCGAA',  # 9 codons
        'medium_sequence': 'ATGGCCGAAATGGCCGAAATGGCCGAAATGGCCGAAATGGCCGAAATGGCCGAA',  # 18 codons
        'long_sequence': 'ATGGCCGAAATGGCCGAAATGGCCGAAATGGCCGAAATGGCCGAAATGGCCGAAATGGCCGAAATGGCCGAAATGGCCGAA',  # 27 codons
        'invalid_length': 'ATGGCCGAAATGGCC',  # Not divisible by 3
        'with_ambiguous': 'ATGGCCGAAATGNCCGAAATGGCCGAA',  # Contains N
        'empty_sequence': '',  # Empty sequence
        'single_codon': 'ATG',  # Single codon
        'two_codons': 'ATGGCC',  # Two codons
    }
    
    # Sample codon orders
    SAMPLE_CODON_ORDERS = {
        'standard_6': ['ATG', 'GCC', 'GAA', 'TAA', 'TAG', 'TGA'],
        'standard_3': ['ATG', 'GCC', 'GAA'],
        'all_64': [
            'AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT',
            'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT',
            'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT',
            'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT',
            'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT',
            'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT',
            'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT',
            'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT'
        ]
    }
    
    # Sample codon frequency data
    SAMPLE_CODON_FREQUENCIES = {
        'focal_region': {
            'ATG': 10, 'GCC': 5, 'GAA': 8, 'TAA': 2, 'TAG': 1, 'TGA': 0
        },
        'background_genome': {
            'ATG': 20, 'GCC': 15, 'GAA': 12, 'TAA': 5, 'TAG': 3, 'TGA': 2
        },
        'all_codons': {
            'ATG': 30, 'GCC': 20, 'GAA': 20, 'TAA': 7, 'TAG': 4, 'TGA': 2
        }
    }
    
    # Sample gene data
    SAMPLE_GENE_DATA = {
        'gene_list': ['gene1', 'gene2', 'gene3', 'gene4', 'gene5'],
        'gene_codons': {
            'gene1': ['ATG', 'GCC', 'GAA'] * 3,  # 9 codons
            'gene2': ['ATG', 'GCC', 'GAA'] * 2,  # 6 codons
            'gene3': ['ATG', 'GCC', 'GAA'] * 3,  # 9 codons
            'gene4': ['TAA', 'TAG', 'TGA'] * 1,  # 3 codons
            'gene5': ['ATG', 'GCC', 'GAA'] * 2   # 6 codons
        },
        'focal_codon_count': 26,
        'total_codon_count': 33
    }
    
    # Sample statistical data
    SAMPLE_STATISTICAL_DATA = {
        'cosine_distance': 0.3,
        'spearman_rho': 0.7,
        'empirical_pvalue': 0.05,
        'simulation_count': 10000
    }
    
    @classmethod
    def get_sample_sequence(cls, name):
        """Get a sample DNA sequence by name."""
        return cls.SAMPLE_DNA_SEQUENCES.get(name, '')
    
    @classmethod
    def get_sample_codon_order(cls, name):
        """Get a sample codon order by name."""
        return cls.SAMPLE_CODON_ORDERS.get(name, [])
    
    @classmethod
    def get_sample_frequencies(cls, region):
        """Get sample codon frequencies for a region."""
        return cls.SAMPLE_CODON_FREQUENCIES.get(region, {})
    
    @classmethod
    def get_sample_gene_data(cls):
        """Get sample gene data."""
        return cls.SAMPLE_GENE_DATA
    
    @classmethod
    def get_sample_statistical_data(cls):
        """Get sample statistical data."""
        return cls.SAMPLE_STATISTICAL_DATA
    
    @classmethod
    def create_codon_count_array(cls, codon_order, frequencies):
        """Create a numpy array of codon counts."""
        return np.array([frequencies.get(codon, 0) for codon in codon_order], dtype=np.int32)
    
    @classmethod
    def create_frequency_array(cls, codon_order, frequencies):
        """Create a numpy array of codon frequencies."""
        return np.array([frequencies.get(codon, 0) for codon in codon_order], dtype=np.float64)
    
    @classmethod
    def create_mock_genbank_record(cls, sequence, features=None):
        """Create a mock GenBank record for testing."""
        if features is None:
            features = []
        
        mock_record = type('MockRecord', (), {})()
        mock_record.seq = type('MockSeq', (), {})()
        mock_record.seq.__str__ = lambda: sequence
        mock_record.features = features
        mock_record.id = "test_scaffold"
        
        return mock_record
    
    @classmethod
    def create_mock_cds_feature(cls, locus_tag, start, end, strand=1):
        """Create a mock CDS feature for testing."""
        mock_feature = type('MockFeature', (), {})()
        mock_feature.type = 'CDS'
        mock_feature.qualifiers = {'locus_tag': [locus_tag]}
        mock_feature.location = f"{start}..{end}" if strand == 1 else f"complement({start}..{end})"
        
        return mock_feature


class TestFixtures:
    """Test fixtures for pytest."""
    
    @staticmethod
    def sample_workflow_data():
        """Create sample data for workflow testing."""
        return {
            'all_cods': {'ATG', 'GCC', 'GAA', 'TAA', 'TAG', 'TGA'},
            'cod_freq_dict_focal': SampleData.get_sample_frequencies('focal_region'),
            'cod_freq_dict_background': SampleData.get_sample_frequencies('background_genome'),
            'gene_list': SampleData.get_sample_gene_data()['gene_list'],
            'gene_codons': SampleData.get_sample_gene_data()['gene_codons'],
            'foc_codon_count': SampleData.get_sample_gene_data()['focal_codon_count'],
            'all_codon_counts': SampleData.get_sample_frequencies('all_codons')
        }
    
    @staticmethod
    def sample_caching_data():
        """Create sample data for caching tests."""
        return {
            'sequence': SampleData.get_sample_sequence('short_sequence'),
            'locus_tag': 'test_gene_1',
            'codon_order': SampleData.get_sample_codon_order('standard_6')
        }
    
    @staticmethod
    def sample_threading_data():
        """Create sample data for threading tests."""
        codon_order = SampleData.get_sample_codon_order('standard_6')
        all_cds_codon_count_arrays = [
            SampleData.create_codon_count_array(codon_order, {'ATG': 5, 'GCC': 2, 'GAA': 3, 'TAA': 1, 'TAG': 0, 'TGA': 0}),
            SampleData.create_codon_count_array(codon_order, {'ATG': 3, 'GCC': 1, 'GAA': 2, 'TAA': 0, 'TAG': 1, 'TGA': 0}),
            SampleData.create_codon_count_array(codon_order, {'ATG': 2, 'GCC': 2, 'GAA': 3, 'TAA': 1, 'TAG': 0, 'TGA': 0})
        ]
        
        return {
            'cpu_simulations': 100,
            'region_freqs_list': [10, 5, 8, 2, 1, 0],
            'codon_order': codon_order,
            'observed_cosine_distance': 0.5,
            'all_cds_codon_count_arrays': all_cds_codon_count_arrays
        }
