import pandas as pd
import requests
from tqdm import tqdm
from pathlib import Path
import typing as t


def _format_req_url(version, method, output_format="tsv-no-header"):
    version = '-'.join(version.split('.')) if version is not None else version
    string_api_url = f"https://version-{version}.string-db.org/api" if version is not None else "https://string-db.org/api"
    output_format = output_format
    method = method
    return "/".join([string_api_url, output_format, method])
    


def fetch_string_ids(genes: t.Sequence,
                    version="12.0",
                    species_id=9606,
                    file_name=None,
                    batch_size=100
                    ) -> pd.DataFrame:
    prev_res = None
    if file_name is not None and Path(file_name).is_file():
        prev_res = pd.read_csv(file_name, index_col=0)
        n_genes = len(genes)
        genes = list(set(genes) - set(prev_res["gene_name"]))
        print(f"Found {n_genes - len(genes)} genes in the existing annotation file")

        if len(genes) == 0:
            return prev_res
    print(f"Trying to get {len(genes)} genes and add them into the existing annotation file")
    
    dfs = []
    request_url = _format_req_url(version=version, method="get_string_ids")
    for i in tqdm(range((len(genes) // batch_size) + 1)):
        params = {
            "identifiers" : "\r".join(genes[i * batch_size: (i+1) * batch_size]), # your protein list
            "species" : species_id, # species NCBI identifier 
            "limit" : 1, # only one (best) identifier per input protein
            "echo_query" : 1, # see your input identifiers in the output
            # "caller_identity" : "www.awesome_app.org" # your app name
        }
        results = requests.post(request_url, data=params)
        if not results.ok:
            continue

        dfs.append(pd.DataFrame([line.split("\t") for line in results.text.split("\n")]).dropna())

    col_id_map = {0: "gene_name", 
                2: "STRING_ID", 
                3: "species", 
                4: "species_name", 
                5: "preferred_name", 
                6: "annotation"}
    if len(dfs) != 0:
        stringdb_ids = pd.concat(dfs, axis=0).dropna(how='any')
        stringdb_ids = stringdb_ids.rename(columns=col_id_map).drop(columns=[1])
    else:
        stringdb_ids = pd.DataFrame(columns=list(col_id_map.values()))
        
    if prev_res is not None:
        stringdb_ids = pd.concat([prev_res, stringdb_ids], axis=0)
    if file_name is not None:
        stringdb_ids.to_csv(file_name)
    return stringdb_ids


def fetch_stringdb(genes, 
                   method = "network", 
                   output="tsv", 
                   img_file_name=None,
                   version="12.0",
                   species_id=9606,
                   **params):
    assert method in ["network", "ppi_enrichment", "interaction_partners"]
    request_url = _format_req_url(version=version, output_format=output, method=method)
    params = {
        "identifiers" : "%0d".join(genes), # your proteins
        "species" : species_id, # species NCBI identifier 
        "caller_identity" : "www.awesome_app.org", # your app name,
        **params
    }
    response = requests.post(request_url, data=params)
    if not response.ok:
        print(f"Request is not successful. Status Code: {response.status_code}")
        print(request_url, params)
        return pd.DataFrame(columns=[])
    
    if output == "image":
        with open(img_file_name, 'wb') as fh:
            fh.write(response.content)
        return
    
    data = [line.split("\t") for line in response.text.split("\n")[1:]]
    columns=response.text.split("\n")[0].split("\t")
    if len(data[0]) != len(columns):
        return pd.DataFrame(columns=columns)
    return pd.DataFrame(data, columns=columns).dropna()


def fetch_enrichment(genes,
                     version="12.0",
                     species_id=9606,):
    request_url = _format_req_url(version=version, output_format="tsv", method="enrichment")
    params = {
        "identifiers" : "%0d".join(genes), # your proteins
        "species" : species_id, # species NCBI identifier 
        "caller_identity" : "www.awesome_app.org" # your app name
    }
    response = requests.post(request_url, data=params)
    data = [line.split("\t") for line in response.text.split("\n")[1:]]
    columns=response.text.split("\n")[0].split("\t")
    if len(data[0]) != len(columns):
        return pd.DataFrame(columns=columns)
    return pd.DataFrame(data, columns=columns).dropna()


def fetch_enrichment_figure(genes, 
                            img_file_name=None,
                            version="12.0",
                            species_id=9606,
                            category="Process",
                            n_terms=20,
                            x_axis: t.Literal["FDR", "signal", "strength", "gene_count"] = "FDR",
                            **params):
    request_url = _format_req_url(version=version, output_format="image", method="enrichmentfigure")
    params = {
        "identifiers" : "%0d".join(genes), # your proteins
        "species" : species_id, # species NCBI identifier 
        "caller_identity" : "www.awesome_app.org", # your app name,
        "number_of_term_shown": n_terms,
        "category": category,
        "x_axis": x_axis,
        **params
    }
    response = requests.post(request_url, data=params)
    if not response.ok:
        print(f"Request is not successful. Status Code: {response.status_code}")
    
    with open(img_file_name, 'wb') as fh:
        fh.write(response.content)
    return