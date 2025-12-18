from ppi_net_builder.src.fetch import fetch_string_ids
from enum import Enum
import typing as t


class Species(Enum):
    human = 9606
    mouse = 10090
    rat = 10116


common_species = [sp.name for sp in Species]


class DataManager:
    def __init__(self,
                 df=None,
                 genes=None,
                 name_col="gene_name",
                 annot_file_name=None,
                 version="12.0",
                 species: t.Union[t.Literal["human", "mouse", "rat"], int, Species]=Species.human):
        self.df = df
        self.genes = genes

        if self.df is not None and self.genes is not None:
            raise ValueError("Can only use either df or genes, got both of them.")

        self._name_col = name_col
        self.annot_file_name = annot_file_name
        self.version = version
        self.species = self.handle_species(species)
        self._gene_map = self.translate_genes(self.df,
                                              name_col=self._name_col,
                                              genes=genes,
                                              annot_file_name=self.annot_file_name,
                                              version=self.version)
    
    @property
    def string_ids(self):
        if self.df is None:
            return [self._gene_map[gn] for gn in self.genes if gn in self._gene_map and self._gene_map[gn] is not None]
        return [self._gene_map[gn] for gn in self.df[self._name_col] if gn in self._gene_map and self._gene_map[gn] is not None]
    
    @staticmethod
    def handle_species(species) -> int:
        if isinstance(species, int):
            return species
        if isinstance(species, Species):
            return int(species.value)
        if isinstance(species, str):
            return int(getattr(Species, species).value)

    @staticmethod
    def translate_genes(df, name_col, genes, annot_file_name, version):
        genes = list(set(df[name_col].to_list())) if genes is None else genes
        annot_df = fetch_string_ids(genes, 
                                    version=version, 
                                    file_name=annot_file_name)
        g_map = {g: None for g in genes}
        g_map.update(dict(zip(annot_df["gene_name"], annot_df["STRING_ID"])))
        return g_map