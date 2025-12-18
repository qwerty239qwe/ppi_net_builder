import pytest
import pandas as pd
from ppi_net_builder.src.data import Species, DataManager, common_species


class TestSpecies:
    """Test cases for the Species enum."""

    def test_species_values(self):
        """Test that species have correct NCBI taxonomy IDs."""
        assert Species.human.value == 9606
        assert Species.mouse.value == 10090
        assert Species.rat.value == 10116
        assert Species.yeast.value == 559292
        assert Species.fruit_fly.value == 7227
        assert Species.zebrafish.value == 7955
        assert Species.chicken.value == 9031
        assert Species.arabidopsis.value == 3702
        assert Species.e_coli.value == 562

    def test_species_names(self):
        """Test that species names are accessible."""
        assert Species.human.name == "human"
        assert Species.mouse.name == "mouse"
        assert Species.rat.name == "rat"

    def test_common_species_list(self):
        """Test that common_species list contains all species names."""
        expected_species = [
            "human", "mouse", "rat", "yeast", "fruit_fly",
            "zebrafish", "chicken", "arabidopsis", "e_coli"
        ]
        assert common_species == expected_species


class TestDataManager:
    """Test cases for the DataManager class."""

    def test_init_with_gene_list(self):
        """Test DataManager initialization with a list of genes."""
        genes = ["TP53", "BRCA1", "EGFR"]
        dm = DataManager(genes=genes, species="human")

        assert dm.genes == genes
        assert dm.df is None
        assert dm.species == 9606  # human NCBI ID
        assert dm.version == "12.0"
        assert dm._name_col == "gene_name"

    def test_init_with_dataframe(self):
        """Test DataManager initialization with a pandas DataFrame."""
        df = pd.DataFrame({"gene_name": ["TP53", "BRCA1", "EGFR"]})
        dm = DataManager(df=df, species="human")

        assert dm.df.equals(df)
        assert dm.genes is None
        assert dm.species == 9606

    def test_init_both_df_and_genes_raises_error(self):
        """Test that providing both df and genes raises ValueError."""
        df = pd.DataFrame({"gene_name": ["TP53"]})
        genes = ["TP53"]

        with pytest.raises(ValueError, match="Can only use either df or genes"):
            DataManager(df=df, genes=genes)

    def test_handle_species_with_enum(self):
        """Test species handling with Species enum."""
        result = DataManager.handle_species(Species.human)
        assert result == 9606

    def test_handle_species_with_int(self):
        """Test species handling with integer NCBI ID."""
        result = DataManager.handle_species(10090)  # mouse
        assert result == 10090

    def test_handle_species_with_string(self):
        """Test species handling with string species name."""
        result = DataManager.handle_species("mouse")
        assert result == 10090

    def test_handle_species_invalid_string(self):
        """Test that invalid species string raises AttributeError."""
        with pytest.raises(AttributeError):
            DataManager.handle_species("invalid_species")

    @pytest.mark.parametrize("species_input,expected", [
        ("human", 9606),
        ("mouse", 10090),
        ("rat", 10116),
        ("yeast", 559292),
        ("fruit_fly", 7227),
        ("zebrafish", 7955),
        ("chicken", 9031),
        ("arabidopsis", 3702),
        ("e_coli", 562),
    ])
    def test_handle_species_parametrized(self, species_input, expected):
        """Test species handling for all supported species."""
        result = DataManager.handle_species(species_input)
        assert result == expected

    def test_string_ids_property_with_df(self):
        """Test string_ids property when initialized with DataFrame."""
        df = pd.DataFrame({"gene_name": ["TP53", "BRCA1"]})
        dm = DataManager(df=df, species="human")

        # Mock the _gene_map since we can't actually call STRING API in tests
        dm._gene_map = {"TP53": "9606.ENSP00000269305", "BRCA1": "9606.ENSP00000471181"}

        expected_ids = ["9606.ENSP00000269305", "9606.ENSP00000471181"]
        assert dm.string_ids == expected_ids

    def test_string_ids_property_with_genes(self):
        """Test string_ids property when initialized with gene list."""
        genes = ["TP53", "BRCA1"]
        dm = DataManager(genes=genes, species="human")

        # Mock the _gene_map
        dm._gene_map = {"TP53": "9606.ENSP00000269305", "BRCA1": "9606.ENSP00000471181"}

        expected_ids = ["9606.ENSP00000269305", "9606.ENSP00000471181"]
        assert dm.string_ids == expected_ids

    def test_string_ids_with_none_values(self):
        """Test string_ids property handles None values in gene_map."""
        genes = ["TP53", "UNKNOWN_GENE", "BRCA1"]
        dm = DataManager(genes=genes, species="human")

        # Mock the _gene_map with some None values
        dm._gene_map = {
            "TP53": "9606.ENSP00000269305",
            "UNKNOWN_GENE": None,
            "BRCA1": "9606.ENSP00000471181"
        }

        expected_ids = ["9606.ENSP00000269305", "9606.ENSP00000471181"]
        assert dm.string_ids == expected_ids
