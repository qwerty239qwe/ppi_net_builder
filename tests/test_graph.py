import pytest
import pandas as pd
from unittest.mock import Mock, patch, MagicMock
from ppi_net_builder.src.graph import NetworkBuilder, require_attribute
from ppi_net_builder.src.data import DataManager


class TestRequireAttribute:
    """Test cases for the require_attribute decorator."""

    def test_require_attribute_success(self):
        """Test that decorator allows access when attribute exists and is not None."""
        class TestClass:
            def __init__(self):
                self.network = "test_network"

            @require_attribute("network")
            def test_method(self):
                return "success"

        obj = TestClass()
        assert obj.test_method() == "success"

    def test_require_attribute_missing_attribute(self):
        """Test that decorator raises AttributeError when attribute doesn't exist."""
        class TestClass:
            @require_attribute("network")
            def test_method(self):
                return "success"

        obj = TestClass()
        with pytest.raises(AttributeError, match="missing the required attribute 'network'"):
            obj.test_method()

    def test_require_attribute_none_value(self):
        """Test that decorator raises ValueError when attribute is None."""
        class TestClass:
            def __init__(self):
                self.network = None

            @require_attribute("network")
            def test_method(self):
                return "success"

        obj = TestClass()
        with pytest.raises(ValueError, match="Please call .construct_network first to get the main network"):
            obj.test_method()

    def test_require_attribute_none_value_generic(self):
        """Test that decorator raises ValueError with generic message for non-network attributes."""
        class TestClass:
            def __init__(self):
                self.other_attr = None

            @require_attribute("other_attr")
            def test_method(self):
                return "success"

        obj = TestClass()
        with pytest.raises(ValueError, match="The attribute 'other_attr' must not be None"):
            obj.test_method()



def test_network_builder_construct_network(monkeypatch, mock_string_ids_data, mock_interaction_data_simple):
    """Test network construction from interaction data."""
    # Monkeypatch the fetch functions
    monkeypatch.setattr(
        'ppi_net_builder.src.data.fetch_string_ids',
        lambda *args, **kwargs: mock_string_ids_data
    )
    monkeypatch.setattr(
        'ppi_net_builder.src.graph.fetch_stringdb',
        lambda *args, **kwargs: mock_interaction_data_simple
    )

    nb = NetworkBuilder(["Gene1", "Gene2", "Gene3"])
    nb.construct_network()

    assert nb.network is not None
    assert nb.network_vert_dic is not None
    assert len(nb.network.vs) == 3  # Should have 3 unique vertices
    assert 'gene_name' in nb.network.vs.attributes()
    assert nb.network.vs['gene_name'] == ['A', 'B', 'C']


def test_network_builder_init_with_gene_list(monkeypatch, mock_string_ids_data, mock_interaction_data_full):
    """Test NetworkBuilder initialization with gene list."""
    # Monkeypatch the fetch functions
    monkeypatch.setattr(
        'ppi_net_builder.src.data.fetch_string_ids',
        lambda *args, **kwargs: mock_string_ids_data
    )
    monkeypatch.setattr(
        'ppi_net_builder.src.graph.fetch_stringdb',
        lambda *args, **kwargs: mock_interaction_data_full
    )

    genes = ["TP53", "BRCA1"]
    nb = NetworkBuilder(genes, species="human")

    assert isinstance(nb.data, DataManager)
    assert nb.network is None
    assert nb.network_vert_dic is None
    assert nb.subnetworks == {}


def test_network_builder_init_with_dataframe(monkeypatch, mock_string_ids_data, mock_interaction_data_full):
    """Test NetworkBuilder initialization with DataFrame."""
    # Monkeypatch the fetch functions
    monkeypatch.setattr(
        'ppi_net_builder.src.data.fetch_string_ids',
        lambda *args, **kwargs: mock_string_ids_data
    )
    monkeypatch.setattr(
        'ppi_net_builder.src.graph.fetch_stringdb',
        lambda *args, **kwargs: mock_interaction_data_full
    )

    df = pd.DataFrame({"gene_name": ["TP53", "BRCA1"]})
    nb = NetworkBuilder(df, species="human")

    assert isinstance(nb.data, DataManager)
    assert nb.network is None


def test_network_builder_construct_network_vertex_mapping(monkeypatch, mock_string_ids_data):
    """Test that vertex mapping is created correctly."""
    # Monkeypatch the fetch functions
    monkeypatch.setattr(
        'ppi_net_builder.src.data.fetch_string_ids',
        lambda *args, **kwargs: mock_string_ids_data
    )
    monkeypatch.setattr(
        'ppi_net_builder.src.graph.fetch_stringdb',
        lambda *args, **kwargs: pd.DataFrame({
            'stringId_A': ['A', 'B'],
            'stringId_B': ['B', 'C'],
            'preferredName_A': ['Gene1', 'Gene2'],
            'preferredName_B': ['Gene2', 'Gene3'],
            'score': [0.9, 0.8]
        })
    )

    nb = NetworkBuilder(["Gene1", "Gene2", "Gene3"])
    nb.construct_network()

    # Check that vertex dictionary maps indices to gene names
    expected_mapping = {0: 'A', 1: 'B', 2: 'C'}
    assert nb.network_vert_dic == expected_mapping


def test_network_builder_extract_subnets(monkeypatch, mock_string_ids_data):
    """Test subnetwork extraction using community detection."""
    # Monkeypatch the fetch functions
    monkeypatch.setattr(
        'ppi_net_builder.src.data.fetch_string_ids',
        lambda *args, **kwargs: mock_string_ids_data
    )
    monkeypatch.setattr(
        'ppi_net_builder.src.graph.fetch_stringdb',
        lambda *args, **kwargs: pd.DataFrame({
            'stringId_A': ['A', 'A', 'B', 'B', 'C', 'D'],
            'stringId_B': ['B', 'C', 'C', 'D', 'D', 'E'],
            'preferredName_A': ['Gene1', 'Gene1', 'Gene2', 'Gene2', 'Gene3', 'Gene4'],
            'preferredName_B': ['Gene2', 'Gene3', 'Gene3', 'Gene4', 'Gene4', 'Gene5'],
            'score': [0.9, 0.8, 0.9, 0.8, 0.7, 0.6]
        })
    )

    nb = NetworkBuilder(["Gene1", "Gene2", "Gene3", "Gene4", "Gene5"])
    nb.construct_network()

    # Mock the community detection method
    from unittest.mock import Mock, patch
    mock_community = Mock()
    mock_clustering = Mock()
    mock_subgraphs = [Mock(), Mock()]  # Two communities
    mock_clustering.subgraphs.return_value = mock_subgraphs
    mock_community.return_value.as_clustering.return_value = mock_clustering

    with patch.object(nb.network, 'community_fastgreedy', mock_community):
        nb.extract_subnets(method="fastgreedy")

    assert len(nb.subnetworks) == 2
    assert 0 in nb.subnetworks
    assert 1 in nb.subnetworks


def test_network_builder_get_enrichment_table_main_network(monkeypatch, mock_string_ids_data, mock_enrichment_data):
    """Test enrichment table retrieval for main network."""
    # Monkeypatch the fetch functions
    monkeypatch.setattr(
        'ppi_net_builder.src.data.fetch_string_ids',
        lambda *args, **kwargs: mock_string_ids_data
    )
    monkeypatch.setattr(
        'ppi_net_builder.src.graph.fetch_stringdb',
        lambda *args, **kwargs: pd.DataFrame({
            'stringId_A': ['A', 'B'],
            'stringId_B': ['B', 'C'],
            'preferredName_A': ['Gene1', 'Gene2'],
            'preferredName_B': ['Gene2', 'Gene3'],
            'score': [0.9, 0.8]
        })
    )
    monkeypatch.setattr(
        'ppi_net_builder.src.fetch.fetch_enrichment',
        lambda *args, **kwargs: mock_enrichment_data
    )

    nb = NetworkBuilder(["Gene1", "Gene2", "Gene3"])
    nb.construct_network()

    result = nb.get_enrichment_table(use_main_network=True)

    assert result.equals(mock_enrichment_data)


def test_network_builder_get_enrichment_table_subnetwork(monkeypatch, mock_string_ids_data, mock_enrichment_data):
    """Test enrichment table retrieval for specific subnetwork."""
    # Monkeypatch the fetch functions
    monkeypatch.setattr(
        'ppi_net_builder.src.data.fetch_string_ids',
        lambda *args, **kwargs: mock_string_ids_data
    )
    monkeypatch.setattr(
        'ppi_net_builder.src.graph.fetch_stringdb',
        lambda *args, **kwargs: pd.DataFrame({
            'stringId_A': ['A', 'B'],
            'stringId_B': ['B', 'C'],
            'preferredName_A': ['Gene1', 'Gene2'],
            'preferredName_B': ['Gene2', 'Gene3'],
            'score': [0.9, 0.8]
        })
    )
    monkeypatch.setattr(
        'ppi_net_builder.src.graph.fetch_enrichment',
        lambda *args, **kwargs: mock_enrichment_data
    )

    nb = NetworkBuilder(["Gene1", "Gene2", "Gene3"])
    nb.construct_network()

    # Mock a subnetwork
    from unittest.mock import Mock
    mock_subnetwork = Mock()
    mock_subnetwork_df = pd.DataFrame({'gene_name': ['Gene1', 'Gene2']})
    mock_subnetwork.get_vertex_dataframe.return_value = mock_subnetwork_df
    nb.subnetworks[0] = mock_subnetwork

    result = nb.get_enrichment_table(use_main_network=False, subnetwork_id=0)

    assert result.equals(mock_enrichment_data)


def test_network_builder_get_enrichment_table_without_network_raises_error(monkeypatch, mock_string_ids_data):
    """Test that enrichment table retrieval fails without constructed network."""
    # Monkeypatch the fetch functions
    monkeypatch.setattr(
        'ppi_net_builder.src.data.fetch_string_ids',
        lambda *args, **kwargs: mock_string_ids_data
    )
    monkeypatch.setattr(
        'ppi_net_builder.src.graph.fetch_stringdb',
        lambda *args, **kwargs: pd.DataFrame()
    )

    nb = NetworkBuilder(["Gene1"])

    with pytest.raises(ValueError, match="Please call .construct_network first to get the main network"):
        nb.get_enrichment_table()


def test_network_builder_save_enrichment_plot(monkeypatch, mock_string_ids_data):
    """Test enrichment plot saving."""
    # Monkeypatch the fetch functions
    monkeypatch.setattr(
        'ppi_net_builder.src.data.fetch_string_ids',
        lambda *args, **kwargs: mock_string_ids_data
    )
    monkeypatch.setattr(
        'ppi_net_builder.src.graph.fetch_stringdb',
        lambda *args, **kwargs: pd.DataFrame({
            'stringId_A': ['A', 'B'],
            'stringId_B': ['B', 'C'],
            'preferredName_A': ['Gene1', 'Gene2'],
            'preferredName_B': ['Gene2', 'Gene3'],
            'score': [0.9, 0.8]
        })
    )
    monkeypatch.setattr(
        'ppi_net_builder.src.graph.fetch_enrichment_figure',
        lambda *args, **kwargs: None
    )

    nb = NetworkBuilder(["Gene1", "Gene2", "Gene3"])
    nb.construct_network()

    nb.save_enrichment_plot("test_plot.png", category="Process", n_terms=10)


def test_network_builder_get_interaction_table_calls_fetch(monkeypatch, mock_string_ids_data, mock_interaction_data_full):
    """Test that get_interaction_table calls fetch_stringdb correctly."""
    # Monkeypatch the fetch functions
    monkeypatch.setattr(
        'ppi_net_builder.src.data.fetch_string_ids',
        lambda *args, **kwargs: mock_string_ids_data
    )
    monkeypatch.setattr(
        'ppi_net_builder.src.graph.fetch_stringdb',
        lambda *args, **kwargs: mock_interaction_data_full
    )

    nb = NetworkBuilder(["Gene1"])
    result = nb.get_interaction_table(required_score=700)

    assert len(result) == 3  # From mock_interaction_data_full
