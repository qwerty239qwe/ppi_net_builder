
# PPI-net-builder
![PyPI version](https://img.shields.io/pypi/v/ppi-net-builder)
![Python version](https://img.shields.io/badge/python-3.9%2B-blue)
![License](https://img.shields.io/pypi/l/ppi-net-builder)
[![codecov](https://codecov.io/github/qwerty239qwe/ppi_net_builder/graph/badge.svg?token=1KG24F4GYU)](https://codecov.io/github/qwerty239qwe/ppi_net_builder)

**PPI-net-builder** is a Python package for constructing and analyzing protein-protein interaction (PPI) networks using the [STRING-db](https://string-db.org/) database.

## Features
- Construct PPI networks from a list of genes.
- Extract high-modality subnetworks.
- Perform enrichment analysis on the main network or its subnetworks.
- Save enrichment plots for visualization.

## Installation
Install the package using pip:

```bash
pip install ppi-net-builder
```

## Usage
Here is an example demonstrating how to use **PPI-net-builder**:

```python
from ppi_net_builder import NetworkBuilder

# List of genes to construct the network
genes = ["p53", "BRCA1", "cdk2", "Q99835"]

# Specify the annotation file name (this will be created later)
annot_file_name = "./annotation.csv"

# Initialize the NetworkBuilder
nb = NetworkBuilder(genes,
                    annot_file_name=annot_file_name,
                    add_color_nodes=10)

# Construct a PPI network and find high-modality subnetworks
nb.construct_network()
nb.extract_subnets()
print(nb.subnets)

# Perform enrichment analysis on the main network or subnetworks
enrich_df = nb.get_enrichment_table()
print(enrich_df.head())

# Save enrichment analysis plots
nb.save_enrichment_plot(img_file_name="enrichment.png")
```

## Requirements
- Python 3.9 or higher
- Dependencies are automatically installed with the package.

## Contributing
Contributions are welcome! Feel free to submit issues or pull requests to improve the package.

## License
This project is licensed under the MIT License. See the LICENSE file for details.

