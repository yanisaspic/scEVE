"""Run this script to prepare the samples before the analysis.

The rows of a prepared sample are gene symbols, except for the last row: it corresponds to a Cell Ontology label.
The names of the columns of a prepared sample correspond to a Cell Ontology id and an integer X separated by a dot (e.g. CL:0000000.X)

    2023/05/04 @yanisaspic"""


import pandas as pd
import vars.cell_types as types
import utils.parser as parser
from collections import Counter


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
