import pytest
import pandas as pd
from unittest.mock import Mock


@pytest.fixture
def sample_gene_list():
    """Sample gene list for testing."""
    return ["TP53", "BRCA1", "EGFR", "MYC"]


@pytest.fixture
def sample_dataframe(sample_gene_list):
    """Sample DataFrame with gene names for testing."""
    return pd.DataFrame({"gene_name": sample_gene_list})


@pytest.fixture
def mock_interaction_data():
    """Mock interaction DataFrame for testing."""
    return pd.DataFrame({
        'stringId_A': ['9606.ENSP00000269305', '9606.ENSP00000471181', '9606.ENSP00000275493'],
        'stringId_B': ['9606.ENSP00000471181', '9606.ENSP00000275493', '9606.ENSP00000418960'],
        'preferredName_A': ['TP53', 'BRCA1', 'EGFR'],
        'preferredName_B': ['BRCA1', 'EGFR', 'MYC'],
        'score': [0.9, 0.8, 0.7]
    })


@pytest.fixture
def mock_enrichment_data():
    """Mock enrichment analysis DataFrame for testing."""
    return pd.DataFrame({
        'term': ['GO:0006915', 'GO:0007049'],
        'description': ['apoptotic process', 'cell cycle'],
        'FDR': [1.23e-15, 2.45e-12],
        'signal': [45.2, 32.1],
        'strength': [78.9, 65.4],
        'gene_count': [15, 12]
    })
