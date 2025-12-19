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


class TestNetworkBuilder:
    """Test cases for the NetworkBuilder class."""

    @patch('ppi_net_builder.src.graph.fetch_stringdb')
    @patch('ppi_net_builder.src.graph.fetch_string_ids')
    def test_init_with_gene_list(self, mock_fetch_ids, mock_fetch_interactions):
        """Test NetworkBuilder initialization with gene list."""
        # Mock fetch_string_ids
        mock_ids_df = pd.DataFrame({
            'gene_name': ['TP53', 'BRCA1'],
            'STRING_ID': ['9606.ENSP00000269305', '9606.ENSP00000471181'],
            'species': ['9606', '9606'],
            'species_name': ['Homo sapiens', 'Homo sapiens'],
            'preferred_name': ['TP53', 'BRCA1'],
            'annotation': ['annotation1', 'annotation2']
        })
        mock_fetch_ids.return_value = mock_ids_df

        # Mock the fetch_stringdb response
        mock_interactions_df = pd.DataFrame({
            'stringId_A': ['9606.ENSP00000269305', '9606.ENSP00000471181'],
            'stringId_B': ['9606.ENSP00000471181', '9606.ENSP00000269305'],
            'preferredName_A': ['TP53', 'BRCA1'],
            'preferredName_B': ['BRCA1', 'TP53'],
            'score': [0.9, 0.8]
        })
        mock_fetch_interactions.return_value = mock_interactions_df

        genes = ["TP53", "BRCA1"]
        nb = NetworkBuilder(genes, species="human")

        assert isinstance(nb.data, DataManager)
        assert nb.network is None
        assert nb.network_vert_dic is None
        assert nb.subnetworks == {}
        mock_fetch_interactions.assert_called_once()

    @patch('ppi_net_builder.src.graph.fetch_stringdb')
    @patch('ppi_net_builder.src.graph.fetch_string_ids')
    def test_init_with_dataframe(self, mock_fetch_ids, mock_fetch_interactions):
        """Test NetworkBuilder initialization with DataFrame."""
        # Mock fetch_string_ids
        mock_ids_df = pd.DataFrame({
            'gene_name': ['TP53', 'BRCA1'],
            'STRING_ID': ['9606.ENSP00000269305', '9606.ENSP00000471181'],
            'species': ['9606', '9606'],
            'species_name': ['Homo sapiens', 'Homo sapiens'],
            'preferred_name': ['TP53', 'BRCA1'],
            'annotation': ['annotation1', 'annotation2']
        })
        mock_fetch_ids.return_value = mock_ids_df

        mock_interactions_df = pd.DataFrame({
            'stringId_A': ['9606.ENSP00000269305'],
            'stringId_B': ['9606.ENSP00000471181'],
            'preferredName_A': ['TP53'],
            'preferredName_B': ['BRCA1'],
            'score': [0.9]
        })
        mock_fetch_interactions.return_value = mock_interactions_df

        df = pd.DataFrame({"gene_name": ["TP53", "BRCA1"]})
        nb = NetworkBuilder(df, species="human")

        assert isinstance(nb.data, DataManager)
        assert nb.network is None
        mock_fetch_interactions.assert_called_once()

    @patch('ppi_net_builder.src.graph.fetch_stringdb')
    @patch('ppi_net_builder.src.graph.fetch_string_ids')
    def test_construct_network(self, mock_fetch_ids, mock_fetch_interactions):
        """Test network construction from interaction data."""
        # Mock fetch_string_ids to return gene-to-ID mapping
        mock_ids_df = pd.DataFrame({
            'gene_name': ['Gene1', 'Gene2', 'Gene3'],
            'STRING_ID': ['A', 'B', 'C'],
            'species': ['9606', '9606', '9606'],
            'species_name': ['Homo sapiens', 'Homo sapiens', 'Homo sapiens'],
            'preferred_name': ['Gene1', 'Gene2', 'Gene3'],
            'annotation': ['annotation1', 'annotation2', 'annotation3']
        })
        mock_fetch_ids.return_value = mock_ids_df

        # Mock fetch_stringdb to return interaction data
        mock_interactions_df = pd.DataFrame({
            'stringId_A': ['A', 'B', 'A'],
            'stringId_B': ['B', 'C', 'C'],
            'preferredName_A': ['Gene1', 'Gene2', 'Gene1'],
            'preferredName_B': ['Gene2', 'Gene3', 'Gene3'],
            'score': [0.9, 0.8, 0.7]
        })
        mock_fetch_interactions.return_value = mock_interactions_df

        nb = NetworkBuilder(["Gene1", "Gene2", "Gene3"])
        nb.construct_network()

        assert nb.network is not None
        assert nb.network_vert_dic is not None
        assert len(nb.network.vs) == 3  # Should have 3 unique vertices
        assert 'gene_name' in nb.network.vs.attributes()
        assert nb.network.vs['gene_name'] == ['A', 'B', 'C']

    @patch('ppi_net_builder.src.graph.fetch_stringdb')
    @patch('ppi_net_builder.src.graph.fetch_string_ids')
    def test_construct_network_vertex_mapping(self, mock_fetch_ids, mock_fetch_interactions):
        """Test that vertex mapping is created correctly."""
        # Mock fetch_string_ids
        mock_ids_df = pd.DataFrame({
            'gene_name': ['Gene1', 'Gene2', 'Gene3'],
            'STRING_ID': ['X', 'Y', 'Z'],
            'species': ['9606', '9606', '9606'],
            'species_name': ['Homo sapiens', 'Homo sapiens', 'Homo sapiens'],
            'preferred_name': ['Gene1', 'Gene2', 'Gene3'],
            'annotation': ['annotation1', 'annotation2', 'annotation3']
        })
        mock_fetch_ids.return_value = mock_ids_df

        mock_interactions_df = pd.DataFrame({
            'stringId_A': ['X', 'Y'],
            'stringId_B': ['Y', 'Z'],
            'preferredName_A': ['Gene1', 'Gene2'],
            'preferredName_B': ['Gene2', 'Gene3'],
            'score': [0.9, 0.8]
        })
        mock_fetch_interactions.return_value = mock_interactions_df

        nb = NetworkBuilder(["Gene1", "Gene2", "Gene3"])
        nb.construct_network()

        # Check that vertex dictionary maps indices to gene names
        expected_mapping = {0: 'X', 1: 'Y', 2: 'Z'}
        assert nb.network_vert_dic == expected_mapping

    @patch('ppi_net_builder.src.graph.fetch_stringdb')
    @patch('ppi_net_builder.src.graph.fetch_string_ids')
    def test_extract_subnets(self, mock_fetch_ids, mock_fetch_interactions):
        """Test subnetwork extraction using community detection."""
        # Mock fetch_string_ids
        mock_ids_df = pd.DataFrame({
            'gene_name': ['Gene1', 'Gene2', 'Gene3', 'Gene4', 'Gene5'],
            'STRING_ID': ['A', 'B', 'C', 'D', 'E'],
            'species': ['9606'] * 5,
            'species_name': ['Homo sapiens'] * 5,
            'preferred_name': ['Gene1', 'Gene2', 'Gene3', 'Gene4', 'Gene5'],
            'annotation': ['ann1', 'ann2', 'ann3', 'ann4', 'ann5']
        })
        mock_fetch_ids.return_value = mock_ids_df

        # Create mock interaction data that forms two clear communities
        mock_interactions_df = pd.DataFrame({
            'stringId_A': ['A', 'A', 'B', 'B', 'C', 'D'],
            'stringId_B': ['B', 'C', 'C', 'D', 'D', 'E'],
            'preferredName_A': ['Gene1', 'Gene1', 'Gene2', 'Gene2', 'Gene3', 'Gene4'],
            'preferredName_B': ['Gene2', 'Gene3', 'Gene3', 'Gene4', 'Gene4', 'Gene5'],
            'score': [0.9, 0.8, 0.9, 0.8, 0.7, 0.6]
        })
        mock_fetch_interactions.return_value = mock_interactions_df

        nb = NetworkBuilder(["Gene1", "Gene2", "Gene3", "Gene4", "Gene5"])
        nb.construct_network()

        # Mock the community detection method
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

    @patch('ppi_net_builder.src.graph.fetch_stringdb')
    @patch('ppi_net_builder.src.graph.fetch_string_ids')
    @patch('ppi_net_builder.src.graph.fetch_enrichment')
    def test_get_enrichment_table_main_network(self, mock_enrichment, mock_fetch_ids, mock_fetch_interactions):
        """Test enrichment table retrieval for main network."""
        # Mock fetch_string_ids
        mock_ids_df = pd.DataFrame({
            'gene_name': ['Gene1', 'Gene2', 'Gene3'],
            'STRING_ID': ['A', 'B', 'C'],
            'species': ['9606'] * 3,
            'species_name': ['Homo sapiens'] * 3,
            'preferred_name': ['Gene1', 'Gene2', 'Gene3'],
            'annotation': ['ann1', 'ann2', 'ann3']
        })
        mock_fetch_ids.return_value = mock_ids_df

        # Setup mock data
        mock_interaction_df = pd.DataFrame({
            'stringId_A': ['A', 'B'],
            'stringId_B': ['B', 'C'],
            'preferredName_A': ['Gene1', 'Gene2'],
            'preferredName_B': ['Gene2', 'Gene3'],
            'score': [0.9, 0.8]
        })
        mock_fetch_interactions.return_value = mock_interaction_df

        mock_enrichment_df = pd.DataFrame({
            'term': ['GO:12345', 'GO:67890'],
            'description': ['Test term 1', 'Test term 2'],
            'FDR': [0.01, 0.05]
        })
        mock_enrichment.return_value = mock_enrichment_df

        nb = NetworkBuilder(["Gene1", "Gene2", "Gene3"])
        nb.construct_network()

        result = nb.get_enrichment_table(use_main_network=True)

        assert result.equals(mock_enrichment_df)
        mock_enrichment.assert_called_once_with(
            genes=['A', 'B', 'C'],
            version='12.0',
            species_id=9606
        )

    @patch('ppi_net_builder.src.graph.fetch_stringdb')
    @patch('ppi_net_builder.src.graph.fetch_string_ids')
    @patch('ppi_net_builder.src.graph.fetch_enrichment')
    def test_get_enrichment_table_subnetwork(self, mock_enrichment, mock_fetch_ids, mock_fetch_interactions):
        """Test enrichment table retrieval for specific subnetwork."""
        # Mock fetch_string_ids
        mock_ids_df = pd.DataFrame({
            'gene_name': ['Gene1', 'Gene2', 'Gene3'],
            'STRING_ID': ['A', 'B', 'C'],
            'species': ['9606'] * 3,
            'species_name': ['Homo sapiens'] * 3,
            'preferred_name': ['Gene1', 'Gene2', 'Gene3'],
            'annotation': ['ann1', 'ann2', 'ann3']
        })
        mock_fetch_ids.return_value = mock_ids_df

        # Setup mock data
        mock_interaction_df = pd.DataFrame({
            'stringId_A': ['A', 'B'],
            'stringId_B': ['B', 'C'],
            'preferredName_A': ['Gene1', 'Gene2'],
            'preferredName_B': ['Gene2', 'Gene3'],
            'score': [0.9, 0.8]
        })
        mock_fetch_interactions.return_value = mock_interaction_df

        nb = NetworkBuilder(["Gene1", "Gene2", "Gene3"])
        nb.construct_network()

        # Mock a subnetwork
        mock_subnetwork = Mock()
        mock_subnetwork_df = pd.DataFrame({'gene_name': ['Gene1', 'Gene2']})
        mock_subnetwork.get_vertex_dataframe.return_value = mock_subnetwork_df
        nb.subnetworks[0] = mock_subnetwork

        mock_enrichment_df = pd.DataFrame({
            'term': ['GO:12345'],
            'description': ['Test term'],
            'FDR': [0.01]
        })
        mock_enrichment.return_value = mock_enrichment_df

        result = nb.get_enrichment_table(use_main_network=False, subnetwork_id=0)

        assert result.equals(mock_enrichment_df)
        mock_enrichment.assert_called_once_with(
            genes=['Gene1', 'Gene2'],
            version='12.0',
            species_id=9606
        )

    @patch('ppi_net_builder.src.graph.fetch_stringdb')
    @patch('ppi_net_builder.src.graph.fetch_string_ids')
    def test_get_enrichment_table_without_network_raises_error(self, mock_fetch_ids, mock_fetch_interactions):
        """Test that enrichment table retrieval fails without constructed network."""
        mock_ids_df = pd.DataFrame({
            'gene_name': ['Gene1'],
            'STRING_ID': ['A'],
            'species': ['9606'],
            'species_name': ['Homo sapiens'],
            'preferred_name': ['Gene1'],
            'annotation': ['ann1']
        })
        mock_fetch_ids.return_value = mock_ids_df

        mock_fetch_interactions.return_value = pd.DataFrame()

        nb = NetworkBuilder(["Gene1"])

        with pytest.raises(ValueError, match="Please call .construct_network first to get the main network"):
            nb.get_enrichment_table()

    @patch('ppi_net_builder.src.graph.fetch_stringdb')
    @patch('ppi_net_builder.src.graph.fetch_string_ids')
    @patch('ppi_net_builder.src.graph.fetch_enrichment_figure')
    def test_save_enrichment_plot(self, mock_plot, mock_fetch_ids, mock_fetch_interactions):
        """Test enrichment plot saving."""
        mock_ids_df = pd.DataFrame({
            'gene_name': ['Gene1', 'Gene2', 'Gene3'],
            'STRING_ID': ['A', 'B', 'C'],
            'species': ['9606'] * 3,
            'species_name': ['Homo sapiens'] * 3,
            'preferred_name': ['Gene1', 'Gene2', 'Gene3'],
            'annotation': ['ann1', 'ann2', 'ann3']
        })
        mock_fetch_ids.return_value = mock_ids_df

        mock_interaction_df = pd.DataFrame({
            'stringId_A': ['A', 'B'],
            'stringId_B': ['B', 'C'],
            'preferredName_A': ['Gene1', 'Gene2'],
            'preferredName_B': ['Gene2', 'Gene3'],
            'score': [0.9, 0.8]
        })
        mock_fetch_interactions.return_value = mock_interaction_df

        nb = NetworkBuilder(["Gene1", "Gene2", "Gene3"])
        nb.construct_network()

        nb.save_enrichment_plot("test_plot.png", category="Process", n_terms=10)

        mock_plot.assert_called_once_with(
            genes=['A', 'B', 'C'],
            img_file_name="test_plot.png",
            version='12.0',
            species_id=9606,
            category="Process",
            n_terms=10,
            x_axis="FDR"
        )

    @patch('ppi_net_builder.src.graph.fetch_stringdb')
    @patch('ppi_net_builder.src.graph.fetch_string_ids')
    def test_get_interaction_table_calls_fetch(self, mock_fetch_ids, mock_fetch_interactions):
        """Test that get_interaction_table calls fetch_stringdb correctly."""
        mock_ids_df = pd.DataFrame({
            'gene_name': ['Gene1'],
            'STRING_ID': ['9606'],
            'species': ['9606'],
            'species_name': ['Homo sapiens'],
            'preferred_name': ['Gene1'],
            'annotation': ['ann1']
        })
        mock_fetch_ids.return_value = mock_ids_df

        mock_df = pd.DataFrame({
            'stringId_A': ['A'],
            'stringId_B': ['B'],
            'score': [0.9]
        })
        mock_fetch_interactions.return_value = mock_df

        nb = NetworkBuilder(["Gene1"])
        result = nb.get_interaction_table(required_score=700)

        assert result.equals(mock_df)
        mock_fetch_interactions.assert_called_with(
            genes=['9606'],  # Mock STRING ID
            method="network",
            version="12.0",
            species_id=9606,
            required_score=700
        )
