"""Functions called to set-up the datasets used in the scEVE paper.

    Run this script after download_data.sh

    2024/04/08 @yanisaspic"""

import os
import pandas as pd

downloads_dir = "./etc/setup_data/source"
data_dir = "./data"


# In the scEVE papers, all cells are associated to a unique id,
# with the following structure: {label}_{n}, where
# label corresponds to the ground truth of the authors.
# n is a value between 1 and the maximum number of cells of a label.

def get_cell_id(cell_label: str, label_counter: dict[str, int]) -> str:
    """Given the label of a cell, assign it a unique id."""
    cell_id = f"{cell_label}_{label_counter[cell_label]}"
    label_counter[cell_label] += 1
    return cell_id


def get_cell_ids(cell_labels: list[str]) -> list[str]:
    """Given a list of cell labels, assign a unique id to each cell w.r.t. its label."""
    label_counter = {label: 1 for label in set(cell_labels)}
    cell_ids = [get_cell_id(cell_label, label_counter) for cell_label in cell_labels]
    print(label_counter)
    return cell_ids


def setup_baron():
    """Set-up the Baron (2016) datasets.
    accession: GSE84133
    sequencing: inDrop
    doi: 10.1016/j.cels.2016.08.011
    
    Baron_HumPan:
    genes: 20,125
    clusters: 14
    	_1: 1,937 cells
    	_2: 1,724 cells
    	_3: 3,605 cells
    	_4: 1,303 cells
    
    Baron_MouPan:
    genes: 14,878
    clusters: 13
    	_1: 822 cells
    	_2: 1,064 cells
    """
    baron_dir = f"{downloads_dir}/Baron"
    paths = os.listdir(baron_dir)
    
    for p in paths:
    	data = pd.read_csv(f"{baron_dir}/{p}", index_col=0)
    	label = p.split("_")[1]
    	
    	species = label[:3].capitalize()
    	number = label[-1]
    	
    	cell_ids = get_cell_ids(data.assigned_cluster)
    	data.index = cell_ids
    	data = data.drop(["barcode", "assigned_cluster"], axis=1)
    	data = data.T
    	
    	out = f"{data_dir}/Baron_{species}Pan_{number}.csv"
    	print(out)
    	print(data)
    	data.to_csv(out)

    return data
    

# Li and Tasic datasets are set-up according to the Hemberg Lab pipeline.
# see:    https://hemberg-lab.github.io/scRNA.seq.datasets/
#         https://github.com/hemberg-lab/scRNA.seq.datasets

def setup_li():
    """Set-up the Li (2017) dataset.
    accession: GSE81861
    cells: 561
    genes: 55,186
    clusters: 9
    sequencing: SMARTer
    doi: 10.1038/ng.3818"""
    li_dir = f"{downloads_dir}/Li"
    data = pd.read_csv(f"{li_dir}/data.csv", index_col=0)

    trim_label = lambda label, sep: label.split(sep)[1]
    data.index = [trim_label(gene_label, "_") for gene_label in data.index]
    data.columns = [trim_label(cell_label, "__") for cell_label in data.columns]
    data.columns = get_cell_ids(data.columns)
    data = data[~data.index.duplicated(keep="first")]
    data.to_csv(f"{data_dir}/Li_HumCRC_a.csv")
    return data


def setup_tasic():
    """Set-up the Tasic (2016) dataset.
    accession: GSE71585
    cells: 1,679
    genes: 24,057
    clusters: 18
    sequencing: SMARTer
    doi: 10.1038/nn.4216
    """
    tasic_dir = f"{downloads_dir}/Tasic"
    data = pd.read_csv(f"{tasic_dir}/genes_counts.csv", index_col=0)

    # column -1: short_name
    cell_metadata = pd.read_csv(f"{tasic_dir}/cell_metadata.csv", index_col=-1)
    short_name_to_long_name = cell_metadata["long_name"].to_dict()
    # column 1: cluster_id
    cluster_metadata = pd.read_csv(f"{tasic_dir}/cluster_metadata.csv", index_col=1)
    cluster_id_to_group = cluster_metadata["group"].to_dict()
    # column 0: short_name
    cell_classification = pd.read_csv(
        f"{tasic_dir}/cell_classification.csv", index_col=0
    )
    short_name_to_cluster_id = cell_classification["primary"].to_dict()
    long_name_to_group = {
        short_name_to_long_name[short_name]: cluster_id_to_group[cluster_id]
        for short_name, cluster_id in short_name_to_cluster_id.items()
    }

    cell_labels = [long_name_to_group[long_name] for long_name in data.columns]
    cell_ids = get_cell_ids(cell_labels)
    data.columns = cell_ids
    data.index = [feature.replace("_", "-") for feature in data.index]
    data.to_csv(f"{data_dir}/Tasic_MouBra.csv")
    return data


# Run after downloading the required files:
setup_baron()
setup_li()
setup_tasic()
