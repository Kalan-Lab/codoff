import pytest
import sys
import os
import numpy as np
from collections import defaultdict

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from tests.fixtures.sample_data import SampleData, TestFixtures


@pytest.fixture
def sample_dna_sequence():
    """Sample DNA sequence for testing."""
    return "ATGGCCGAAATGGCCGAAATGGCCGAA"


@pytest.fixture
def sample_codon_order():
    """Sample codon order for testing."""
    return ['ATG', 'GCC', 'GAA', 'TAA', 'TAG', 'TGA']


@pytest.fixture
def sample_codon_frequencies():
    """Sample codon frequencies for testing."""
    return {
        'focal': {'ATG': 10, 'GCC': 5, 'GAA': 8, 'TAA': 2, 'TAG': 1, 'TGA': 0},
        'background': {'ATG': 20, 'GCC': 15, 'GAA': 12, 'TAA': 5, 'TAG': 3, 'TGA': 2}
    }


@pytest.fixture
def sample_gene_data():
    """Sample gene data for testing."""
    return {
        'gene_list': ['gene1', 'gene2', 'gene3'],
        'gene_codons': {
            'gene1': ['ATG', 'GCC', 'GAA'] * 3,
            'gene2': ['ATG', 'GCC', 'GAA'] * 2,
            'gene3': ['ATG', 'GCC', 'GAA'] * 3
        },
        'foc_codon_count': 24
    }


@pytest.fixture
def sample_workflow_data():
    """Sample workflow data for integration testing."""
    return TestFixtures.sample_workflow_data()


@pytest.fixture
def sample_caching_data():
    """Sample caching data for testing."""
    return TestFixtures.sample_caching_data()


@pytest.fixture
def sample_threading_data():
    """Sample threading data for testing."""
    return TestFixtures.sample_threading_data()


@pytest.fixture
def mock_genbank_record():
    """Mock GenBank record for testing."""
    return SampleData.create_mock_genbank_record("ATGGCCGAAATGGCCGAAATGGCCGAA")


@pytest.fixture
def mock_cds_feature():
    """Mock CDS feature for testing."""
    return SampleData.create_mock_cds_feature("test_gene_1", 1, 27, 1)


@pytest.fixture
def temp_output_file():
    """Temporary output file for testing."""
    import tempfile
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as tmp_file:
        yield tmp_file.name
    os.unlink(tmp_file.name)


@pytest.fixture
def temp_plot_file():
    """Temporary plot file for testing."""
    import tempfile
    with tempfile.NamedTemporaryFile(suffix='.svg', delete=False) as tmp_file:
        yield tmp_file.name
    os.unlink(tmp_file.name)


@pytest.fixture(autouse=True)
def clear_caches():
    """Clear caches before each test."""
    from codoff.codoff import _clear_caches
    _clear_caches()
    yield
    _clear_caches()