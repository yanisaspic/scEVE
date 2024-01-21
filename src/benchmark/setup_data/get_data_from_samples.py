"""Run this script to prepare the samples before the analysis.

The rows of a prepared sample are gene symbols, except for the last row: it corresponds to a Cell Ontology label.
The names of the columns of a prepared sample correspond to a Cell Ontology id and an integer X separated by a dot (e.g. CL:0000000.X)

    2023/05/04 @yanisaspic"""


import pandas as pd
import ground_truth as types
from collections import Counter
import xml.etree.ElementTree as ET
import pandas as pd


def get_smartseq2_sample_metadata(element: ET.Element) -> tuple[str, str]:
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


def get_smartseq2_metadata(path: str) -> dict[str, str]:
    """Read a Smart-Seq2 miniML file to get a dict associating a molecular profile id to its manually assigned cell type."""
    tree = ET.parse(path)
    metadata = {}
    for element in tree.getroot():
        tag = element.tag.split("}")[1]  # ignore miniML prefix
        if tag == "Sample":
            cell_id, cell_type = get_smartseq2_sample_metadata(element)
            metadata[cell_id] = cell_type
    return metadata


def get_unique_cl_id(
    id: str,
    metadata: dict[str, str],
    cell_types_ids: dict[str, str],
    counter: dict[str, int],
) -> pd.DataFrame:
    """Given the unique id of a molecular profile, return a valid unique CL id."""
    cell_type = metadata[id]
    unique_cl_id = f"{cell_types_ids[cell_type]}_{counter[cell_type]}"
    counter[cell_type] += 1
    return unique_cl_id


def apply_unique_cl_ids(
    data: pd.DataFrame, metadata: dict[str, str], cell_types_ids: dict[str, str]
) -> pd.DataFrame:
    """Replace the unique ids of each molecular profile in a scRNA-seq matrix by using the Cell Ontology (CL) ids."""
    counter = {cell_type: 1 for cell_type in metadata.values()}
    unique_cl_ids = [
        get_unique_cl_id(col, metadata, cell_types_ids, counter) for col in data.columns
    ]
    data.columns = unique_cl_ids
    return data


DARMANIS_DATA_PATH = "samples/Darmanis_GBM.csv"
DARMANIS_METADATA_PATH = "samples/Darmanis_GBM.xml"
JERBYARNON_DATA_PATH = "samples/Jerby-Arnon_MLM.csv"
JERBYARNON_METADATA_PATH = "samples/Jerby-Arnon_MLM.xml"


#############################
#   DARMANIS et al. (2017)  #
#   Cancer: Glioblastoma    #
#   Technology: Smart-Seq2  #
#############################
darmanis_data = pd.read_csv(DARMANIS_DATA_PATH, index_col=0, sep=" ")
darmanis_metadata = parser.get_smartseq2_metadata(DARMANIS_METADATA_PATH)
darmanis_cell_distribution = Counter(darmanis_metadata.values())

# rows
darmanis_data.drop(darmanis_data.tail(5).index, inplace=True)
# columns
darmanis_data.columns = darmanis_metadata.keys()
darmanis_data = parser.apply_unique_cl_ids(
    darmanis_data, darmanis_metadata, types.DARMANIS_CELL_TYPES
)
for axis in [0, 1]:
    darmanis_data.sort_index(axis=axis, inplace=True)
darmanis_data.to_csv("data/Darmanis_GBM.csv", index=True)


#################################
#   JERBY-ARNON et al. (2018)   #
#   Cancer: Melanoma            #
#   Technology: Smart-Seq2      #
#################################
jerbyarnon_data = pd.read_csv(JERBYARNON_DATA_PATH, index_col=0)
jerbyarnon_metadata = parser.get_smartseq2_metadata(JERBYARNON_METADATA_PATH)
jerbyarnon_cell_distribution = Counter(jerbyarnon_metadata.values())

# rows: prevent errors wrt gene APOBEC3A_B
jerbyarnon_data.index = jerbyarnon_data.index.str.replace("_", "-")
# columns
jerbyarnon_data.drop(
    [col for col in jerbyarnon_metadata if jerbyarnon_metadata[col] == "Undetermined"],
    axis=1,
    inplace=True,
)
jerbyarnon_data = parser.apply_unique_cl_ids(
    jerbyarnon_data, jerbyarnon_metadata, types.JERBYARNON_CELL_TYPES
)
for axis in [0, 1]:
    jerbyarnon_data.sort_index(axis=axis, inplace=True)
jerbyarnon_data.to_csv("data/Jerby-Arnon_MLM.csv", index=True)
