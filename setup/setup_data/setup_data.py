"""Functions called to set-up the datasets used in the scEVE paper.

    Run this script after download_data.sh

    2024/04/08 @yanisaspic"""

import os
import pandas as pd


baron_dir = "./downloads/Baron"
li_dir = "./downloads/Li"
tasic_dir = "./downloads/Tasic"
camp_dir = "./downloads/Camp"
lake_dir = "./downloads/Lake"


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
    return cell_ids


# scEFSC datasets are set-up according to the Hemberg Lab pipeline.
# see:    https://hemberg-lab.github.io/scRNA.seq.datasets/
#         https://github.com/hemberg-lab/scRNA.seq.datasets

def setup_baron(baron_dir):
    """Set-up the Baron (2016) dataset.
    accession: GSE84133
    cells: 8,569
    genes: 20,125
    clusters: 14
    sequencing: inDrop
    doi: 10.1016/j.cels.2016.08.011
    """
    paths = os.listdir(baron_dir)
    datasets = [pd.read_csv(f"{baron_dir}/{p}", index_col=0) for p in paths]
    data = pd.concat(datasets)
    # paths must include only the 4 human datasets: GSM2230757_humanX_umifm_counts.csv,
    # with X ranging from 1 to 4.

    cell_ids = get_cell_ids(data.assigned_cluster)
    data.index = cell_ids
    data = data.drop(["barcode", "assigned_cluster"], axis=1)
    data = data.T
    data.to_csv("../../data/Baron_HumPan.csv")
    return data


def setup_li(li_dir):
    """Set-up the Li (2017) dataset.
    accession: GSE81861
    cells: 561
    genes: 55,186
    clusters: 9
    sequencing: SMARTer
    doi: 10.1038/ng.3818"""
    data = pd.read_csv(f"{li_dir}/data.csv", index_col=0)

    trim_label = lambda label, sep: label.split(sep)[1]
    data.index = [trim_label(gene_label, "_") for gene_label in data.index]
    data.columns = [trim_label(cell_label, "__") for cell_label in data.columns]
    data.columns = get_cell_ids(data.columns)
    data = data[~data.index.duplicated(keep="first")]
    data.to_csv("../../data/Li_HumCRC_a.csv")
    return data


def setup_tasic(tasic_dir):
    """Set-up the Tasic (2016) dataset.
    accession: GSE71585
    cells: 1,679
    genes: 24,057
    clusters: 18
    sequencing: SMARTer
    doi: 10.1038/nn.4216
    """
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
    data.to_csv("../../data/Tasic_MouBra.csv")
    return data


def setup_camp(camp_dir):
    """Set-up the Camp (2017) dataset.
    accession: GSE81252
    cells: 777
    genes: 19,020
    clusters: 7
    sequencing: SMARTer
    doi: 10.1038/nature22796
    """
    datasets = [pd.read_csv(f"{camp_dir}/data_1.csv", index_col=0),
                pd.read_csv(f"{camp_dir}/data_2.csv", index_col=0)]
    data = pd.concat(datasets)
    data = data[~data.index.duplicated()]

    cell_symbols_1 = data.experiment
    cell_symbols_2 = data.assignment_LB
    cell_symbols_1[~cell_symbols_2.isna()] = cell_symbols_2
    data = data.drop(["experiment", "assignment_LB"], axis=1)

    symbol_to_label = {"ih": "immature hepatoblast",
                       "mh": "mature hepatocyte",
                       "de": "definitive endoderm",
                       "ec": "endothelial",
                       "he": "hepatic endoderm",
                       "mc": "mesenchymal stem cell",
                       "ip": "ipsc"}
    cell_labels = cell_symbols_1.apply(lambda symbol: symbol_to_label[symbol[:2].lower()])
    cell_ids = get_cell_ids(cell_labels)

    data.index = cell_ids
    data = data.T
    data.to_csv("../../data/Camp_MouLiv.csv")
    return data


def setup_lake(lake_dir):
    """Set-up the Lake (2017) dataset.
    accession: phs000833.v3.p1
    cells: 3,042
    genes: 25,051
    clusters: 16
    sequencing: Fluidigm C1
    doi: 10.1126/science.aaf1204 
    """
    data = pd.read_csv(f"{lake_dir}/data.csv", index_col=0, sep="\t")
    data = data[~data.index.duplicated()]

    # column 1: Published_Sample_Name
    annotations = pd.read_csv(f"{lake_dir}/annotations.txt", index_col=1, sep="\t")
    cells_of_interest = annotations.index.intersection(data.columns)
    annotations = annotations.loc[cells_of_interest]
    data = data[cells_of_interest]

    cell_ids = get_cell_ids(annotations.SubGroup)
    data.columns = cell_ids
    data.to_csv("../../data/Lake_MouBra.csv")
    return data
    

def setup_scEFSC_datasets():
    """Set-up the datasets used in the scEFSC paper."""
    setup_baron(baron_dir)
    setup_li(li_dir)
    setup_tasic(tasic_dir)
    setup_camp(camp_dir)
    setup_lake(lake_dir)


# Run after downloading the required files:
setup_lake(lake_dir)