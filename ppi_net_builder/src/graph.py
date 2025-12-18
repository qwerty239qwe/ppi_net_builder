from ppi_net_builder.src.data import DataManager
from ppi_net_builder.src.fetch import fetch_stringdb, fetch_enrichment, fetch_enrichment_figure
import pandas as pd
import igraph as ig
from igraph import Graph

import typing as t


def require_attribute(attr_name):
    """
    Decorator to check if a specific attribute exists on the instance
    and is not None before calling the method.
    """
    def decorator(method):
        def wrapper(self, *args, **kwargs):
            if not hasattr(self, attr_name):
                raise AttributeError(f"The object is missing the required attribute '{attr_name}'.")
            if getattr(self, attr_name) is None:
                error_msg = "Please call .construct_network first to get the main network" if attr_name == "network" else f"The attribute '{attr_name}' must not be None."
                raise ValueError(error_msg)
            return method(self, *args, **kwargs)
        return wrapper
    return decorator


class NetworkBuilder:
    """A class for building and analyzing protein-protein interaction (PPI) networks.

    This class provides functionality to construct PPI networks from gene lists or DataFrames,
    extract subnetworks using community detection algorithms, and perform functional enrichment
    analysis. It integrates with STRING database for interaction data and enrichment analysis.

    Attributes:
        data (DataManager): Manages gene data, annotations, and species information.
        interaction_table (pd.DataFrame): DataFrame containing protein-protein interactions.
        network (Optional[Graph]): The main igraph Graph object representing the PPI network.
        network_vert_dic (Optional[Dict[int, str]]): Mapping from vertex indices to gene names.
        subnetworks (Dict[int, Graph]): Dictionary of extracted subnetworks indexed by ID.
    """

    def __init__(self,
                 genes_or_gene_df: t.Union[pd.DataFrame, t.List[str]],
                 gene_col: str = "gene_name",
                 annot_file_name: t.Optional[str] = None,
                 species: str = "human",
                 version: str = "12.0",
                 add_color_nodes: int = 10,
                 **kwargs) -> None:
        """Initialize the NetworkBuilder with gene data and configuration.

        Args:
            genes_or_gene_df: Either a pandas DataFrame containing gene information
                or a list of gene names. If DataFrame, must contain the gene column.
            gene_col: Name of the column containing gene names in the DataFrame.
                Defaults to "gene_name".
            annot_file_name: Optional path to a local annotation file. If None,
                annotations will be fetched from STRING database.
            species: Species name or NCBI taxonomy ID. Defaults to "human".
            version: STRING database version to use. Defaults to "12.0".
            add_color_nodes: Number of additional nodes to add for coloring in
                network visualization. Defaults to 10.
            **kwargs: Additional keyword arguments passed to STRING API calls.

        Raises:
            ValueError: If both DataFrame and genes are provided in conflicting ways.
        """
        self.data: DataManager = DataManager(df=genes_or_gene_df if isinstance(genes_or_gene_df, pd.DataFrame) else None,
                                             annot_file_name=annot_file_name,
                                             genes=genes_or_gene_df if not isinstance(genes_or_gene_df, pd.DataFrame) else None,
                                             name_col=gene_col,
                                             version=version,
                                             species=species)
        self.interaction_table: pd.DataFrame = self.get_interaction_table(add_color_nodes=add_color_nodes, **kwargs)
        self.network: t.Optional[Graph] = None
        self.network_vert_dic: t.Optional[t.Dict[int, str]] = None
        self.subnetworks: t.Dict[int, Graph] = {}

    def save_interaction_image(self,
                               img_file_name: str,
                               use_main_network: bool = True,
                               subnetwork_id: t.Optional[int] = None,
                               **kwargs) -> None:
        """Save a network visualization image to file.

        Args:
            img_file_name: Path and filename for the output image file.
            use_main_network: If True, visualize the main network. If False,
                visualize a specific subnetwork. Defaults to True.
            subnetwork_id: ID of the subnetwork to visualize when use_main_network
                is False. Required if use_main_network is False.
            **kwargs: Additional keyword arguments passed to the STRING API.

        Raises:
            ValueError: If use_main_network is False but subnetwork_id is None,
                or if the specified subnetwork doesn't exist.
        """
        if use_main_network:
            _ = fetch_stringdb(genes=self.data.string_ids,
                               method="network",
                               output="image",
                               img_file_name=img_file_name,
                               version=self.data.version,
                               species_id=self.data.species,
                               **kwargs)
        fetch_stringdb(genes=self.subnetworks[subnetwork_id].get_vertex_dataframe()["gene_name"].to_list(),
                       method="network",
                       output="image",
                       img_file_name=img_file_name,
                       version=self.data.version,
                       species_id=self.data.species,
                       **kwargs)
    
    def get_interaction_table(self, **kwargs) -> pd.DataFrame:
        """Fetch protein-protein interaction data from STRING database.

        Args:
            **kwargs: Additional keyword arguments passed to the STRING API.
                Common parameters include:
                - required_score: Minimum interaction score threshold (0-1000)
                - network_type: Type of interactions ("functional" or "physical")
                - add_nodes: Number of additional nodes to include

        Returns:
            pd.DataFrame: DataFrame containing interaction data with columns:
                - stringId_A: STRING ID of first protein
                - stringId_B: STRING ID of second protein
                - preferredName_A: Name of first protein
                - preferredName_B: Name of second protein
                - score: Interaction confidence score
        """
        return fetch_stringdb(genes=self.data.string_ids,
                              method="network",
                              version=self.data.version,
                              species_id=self.data.species,
                              **kwargs)

    def construct_network(self) -> None:
        """Construct an igraph Graph object from the interaction data.

        This method creates a network graph from the protein-protein interaction
        data fetched from STRING. Each unique protein becomes a node (vertex),
        and each interaction becomes an edge.

        The resulting graph is stored in self.network and includes vertex attributes
        for gene names. A mapping between vertex indices and gene names is also
        created and stored in self.network_vert_dic.
        """
        all_gene_names = sorted(list(set(self.interaction_table["stringId_A"]) | set(self.interaction_table["stringId_B"])))
        vertex_dic = {g: i for i, g in enumerate(all_gene_names)}
        edges = [[vertex_dic[a], vertex_dic[b]] for a, b in zip(self.interaction_table["stringId_A"],
                                                                self.interaction_table["stringId_B"])]

        self.network, self.network_vert_dic = Graph(edges=edges), {v: k for k, v in vertex_dic.items()}
        self.network.vs["gene_name"] = list(self.network_vert_dic.values())

    @require_attribute("network")
    def extract_subnets(self, method: str = "fastgreedy", **kwargs) -> None:
        """Extract subnetworks using community detection algorithms.

        This method applies community detection to identify densely connected
        subnetworks within the main PPI network. The detected communities are
        sorted by size (largest first) and stored in self.subnetworks.

        Args:
            method: Community detection algorithm to use. Supported methods
                include "fastgreedy", "multilevel", "infomap", "leading_eigenvector",
                "walktrap", etc. Defaults to "fastgreedy".
            **kwargs: Additional keyword arguments passed to the community
                detection algorithm.

        Raises:
            ValueError: If the network has not been constructed yet.

        Note:
            Requires that construct_network() has been called first.
        """
        vdg = getattr(self.network, f"community_{method}")(**kwargs)
        vclusters = vdg.as_clustering() if method == 'fastgreedy' else vdg
        subnetworks = list(reversed(sorted(vclusters.subgraphs(), key=lambda x: x.vcount())))
        self.subnetworks = {i: net for i, net in enumerate(subnetworks)}
    
    @require_attribute("network")
    def get_enrichment_table(self,
                             use_main_network: bool = True,
                             subnetwork_id: t.Optional[int] = None) -> pd.DataFrame:
        """Fetch functional enrichment analysis results from STRING database.

        Performs gene ontology and pathway enrichment analysis on the genes
        in the main network or a specific subnetwork.

        Args:
            use_main_network: If True, analyze the main network. If False,
                analyze a specific subnetwork. Defaults to True.
            subnetwork_id: ID of the subnetwork to analyze when use_main_network
                is False. Required if use_main_network is False.

        Returns:
            pd.DataFrame: DataFrame containing enrichment results with columns
                including term ID, description, FDR, and other enrichment metrics.

        Raises:
            ValueError: If the network has not been constructed yet, or if
                use_main_network is False but subnetwork_id is None/invalid.

        Note:
            Requires that construct_network() has been called first.
        """
        if use_main_network:
            genes_to_use = self.network.get_vertex_dataframe()["gene_name"].to_list()
        else:
            genes_to_use = self.subnetworks[subnetwork_id].get_vertex_dataframe()["gene_name"].to_list()
        df = fetch_enrichment(genes=genes_to_use,
                              version=self.data.version,
                              species_id=self.data.species)
        return df
    
    @require_attribute("network")
    def save_enrichment_plot(self,
                             img_file_name: str,
                             use_main_network: bool = True,
                             subnetwork_id: t.Optional[int] = None,
                             category: str = "Process",
                             n_terms: int = 20,
                             x_axis: t.Literal["FDR", "signal", "strength", "gene_count"] = "FDR",
                             **kwargs) -> None:
        """Save functional enrichment analysis plot to file.

        Creates a visualization plot of the enrichment analysis results and saves
        it to the specified file. The plot shows the most significant terms from
        the chosen category.

        Args:
            img_file_name: Path and filename for the output plot image file.
            use_main_network: If True, analyze the main network. If False,
                analyze a specific subnetwork. Defaults to True.
            subnetwork_id: ID of the subnetwork to analyze when use_main_network
                is False. Required if use_main_network is False.
            category: Enrichment category to plot. Options include "Process",
                "Component", "Function", "KEGG", "RCTM", etc. Defaults to "Process".
            n_terms: Number of top enrichment terms to include in the plot.
                Defaults to 20.
            x_axis: Metric to use for the x-axis. Options are "FDR", "signal",
                "strength", "gene_count". Defaults to "FDR".
            **kwargs: Additional keyword arguments passed to the plotting function.

        Raises:
            ValueError: If the network has not been constructed yet, or if
                use_main_network is False but subnetwork_id is None/invalid.

        Note:
            Requires that construct_network() has been called first.
        """
        if use_main_network:
            genes_to_use = self.network.get_vertex_dataframe()["gene_name"].to_list()
        else:
            genes_to_use = self.subnetworks[subnetwork_id].get_vertex_dataframe()["gene_name"].to_list()

        fetch_enrichment_figure(genes=genes_to_use,
                                img_file_name=img_file_name,
                                version=self.data.version,
                                species_id=self.data.species,
                                category=category,
                                n_terms=n_terms,
                                x_axis=x_axis,
                                **kwargs)
