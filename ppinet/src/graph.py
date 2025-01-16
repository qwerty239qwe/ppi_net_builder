from ppinet.src.data import DataManager
from ppinet.src.fetch import fetch_stringdb, fetch_enrichment, fetch_enrichment_figure
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
    def __init__(self,
                 genes_or_gene_df,
                 gene_col="gene_name",
                 annot_file_name=None,
                 species="human",
                 version="12.0",
                 add_color_nodes=10,
                 **kwargs):
        self.data = DataManager(df=genes_or_gene_df if isinstance(genes_or_gene_df, pd.DataFrame) else None,
                                annot_file_name=annot_file_name,
                                genes=genes_or_gene_df if not isinstance(genes_or_gene_df, pd.DataFrame) else None,
                                name_col=gene_col,
                                version=version,
                                species=species)
        self.interaction_table = self.get_interaction_table(add_color_nodes=add_color_nodes, **kwargs)
        self.network = None
        self.network_vert_dic = None
        self.subnetworks = {}

    def save_interaction_image(self, 
                               img_file_name,
                               use_main_network=True, 
                               subnetwork_id=None,
                               **kwargs):
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
    
    def get_interaction_table(self, **kwargs):
        return fetch_stringdb(genes=self.data.string_ids, 
                              method="network", 
                              version=self.data.version, 
                              species_id=self.data.species,
                              **kwargs)

    def construct_network(self):
        all_gene_names = sorted(list(set(self.interaction_table["stringId_A"]) | set(self.interaction_table["stringId_B"])))
        vertex_dic = {g: i for i, g in enumerate(all_gene_names)}
        edges = [[vertex_dic[a], vertex_dic[b]] for a, b in zip(self.interaction_table["stringId_A"], 
                                                                self.interaction_table["stringId_B"])]
        
        self.network, self.network_vert_dic = Graph(edges=edges), {v: k for k, v in vertex_dic.items()}
        self.network.vs["gene_name"] = list(self.network_vert_dic.values())

    @require_attribute("network")
    def extract_subnets(self, method="fastgreedy", **kwargs):
        vdg = getattr(self.network, f"community_{method}")(**kwargs)
        vclusters = vdg.as_clustering() if method == 'fastgreedy' else vdg
        subnetworks = list(reversed(sorted(vclusters.subgraphs(), key=lambda x: x.vcount())))
        self.subnetworks = {i: net for i, net in enumerate(subnetworks)}
    
    @require_attribute("network")
    def get_enrichment_table(self,
                             use_main_network=True, 
                             subnetwork_id=None,):
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
                             img_file_name,
                             use_main_network=True, 
                             subnetwork_id=None,
                             category="Process",
                             n_terms=20,
                             x_axis: t.Literal["FDR", "signal", "strength", "gene_count"] = "FDR",
                             **kwargs):
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