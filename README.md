# PPI-net-builder v0.1.0
Protein-Protein Interaction (PPI) network construction and analysis using STRING-db (https://string-db.org/) in Python.


## Installation
```bash
pip install ppi-net-builder
```

## Usage
```python
from ppi_net_builder import NetworkBuilder

genes = ["p53", "BRCA1", "cdk2", "Q99835"]
annot_file_name = "./annotation.csv"  # this will be created later
nb = NetworkBuilder(genes,
                    annot_file_name=annot_file_name,
                    add_color_nodes=10)

# construct a ppi network, and find high-modality subnetworks
nb.construct_network()
nb.extract_subnets()
print(nb.subnets)

# do enrichment analysis on the main network or subnets
enrich_df = nb.get_enrichment_table()
print(enrich_df.head())


nb.save_enrichment_plot(img_file_name="enrichment.png")

```

