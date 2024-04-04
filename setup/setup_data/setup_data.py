"""Functions called to set-up the datasets used in the scEVE paper.

    Run this script after download_data.sh

    2024/04/04 @yanisaspic"""

import os
import pandas as pd
import xml.etree.ElementTree as ET

darmanis_dir = "./downloads/Darmanis"
baron_dir = "./downloads/Baron"
li_dir = "./downloads/Li"


# In the scEVE paper, all cells are associated to a unique id,
# with the following structure: {label}_{n}, where
# label corresponds to the ground truth of the authors.
# n is a value between 1 and the maximum number of cells of a label.
#__________________________________________________________________________
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


# Darmanis dataset is set-up with the cell types (e.g. Myeloid)
# assigned by the authors as labels. Note that for each cell,
# the t-SNE clusters assigned are also available in the metadata.
#__________________________________________________________________________
def get_smartseq2_metadata(element: ET.Element) -> tuple[str, str]:
    """Read an element tree of a Smart-Seq2 miniML file to get the cell id and the cell type corresponding to a sample."""
    for child in element:
        tag = child.tag.split("}")[1]  # ignore miniML prefix
        if tag == "Title":
            cell_id = str(child.text).strip()
        elif tag == "Channel":
            element_of_interest = child
            break
    characteristics = {}
    for child in element_of_interest:
        if "tag" in child.attrib:
            characteristics[child.attrib["tag"]] = str(child.text).strip()
    return cell_id, characteristics["cell type"]


def GET_SMARTSEQ2_METADATA(path: str) -> dict[str, str]:
    """Read a Smart-Seq2 miniML file to get a dict associating a molecular profile id to its manually assigned cell type."""
    tree = ET.parse(path)
    metadata = {}
    for element in tree.getroot():
        tag = element.tag.split("}")[1]  # ignore miniML prefix
        if tag == "Sample":
            cell_id, cell_type = get_smartseq2_metadata(element)
            metadata[cell_id] = cell_type
    return metadata


def setup_darmanis(darmanis_dir):
    """Set-up the Darmanis (2017) dataset.
    accession: GSE84465
    cells: 3,589
    genes: 23,460
    clusters: 7
    sequencing: Smart-Seq2
    doi: 10.1016/j.celrep.2017.10.030
    """
    data = pd.read_csv(
        f"{darmanis_dir}/data.csv", index_col=0, sep=" "
    )
    data = data.drop(data.tail(5).index)
    # drop the rows 'no_feature', 'ambiguous', 'too_low_aQual', 'not_aligned' and 'alignment_not_unique'

    metadata = GET_SMARTSEQ2_METADATA(f"{darmanis_dir}/GSE84465_family.xml")
    cell_ids = get_cell_ids(metadata.values())
    data.columns = cell_ids
    data.to_csv("../../data/Darmanis_HumGBM.csv")


# scEFSC datasets are set-up according to the Hemberg Lab pipeline.
# see:    https://hemberg-lab.github.io/scRNA.seq.datasets/
#         https://github.com/hemberg-lab/scRNA.seq.datasets
#__________________________________________________________________________
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
    data.to_csv("../../data/Li_HumCRC.csv")


# Run after downloading the required files.
#__________________________________________________________________________
setup_darmanis(darmanis_dir)
setup_baron(baron_dir)
setup_li(li_dir)