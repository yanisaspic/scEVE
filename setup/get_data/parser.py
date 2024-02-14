"""Library to get the ground truth for the Darmanis, Jerby-Arnon and Tirosh datasets (Smart-Seq2 cancer).

    2023/01/25 @yanisaspic"""

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
