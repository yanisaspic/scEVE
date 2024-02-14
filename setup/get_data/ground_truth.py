"""Associate non-standardised cell labels in a dataset to its standard id in the Cell Ontology.

    2023/05/05 @yanisaspic"""

MOCK_CELL_TYPES = {"Oligodendrocyte": "CL:0000128",
                   "Neoplastic": "CL:0001063"}

DARMANIS_CELL_TYPES = {
    "Immune cell": "Mye",  # "CL:0000763",  # myeloid cells
    "Oligodendrocyte": "Oli",  # "CL:0000128",
    "Vascular": "Vas",  # "CL:0002139",  # endothelial cell of vascular tree
    "Neuron": "Neu",  # "CL:0000540",
    "Astocyte": "Ast",  # "CL:0000127",  # astrocyte
    "OPC": "OLP",  # "CL:0002453",  # oligodendrocyte precursor cells
    "Neoplastic": "Can",  # "CL:0001063",  # neoplastic cell
}  # see https://doi.org/10.1016/j.celrep.2017.10.030

JERBYARNON_CELL_TYPES = {
    "T.CD4": "CD4",  # "CL:0000624",  # CD4-positive, alpha-beta T cell
    "T.cell": "T",  # "CL:0000084",
    "B.cell": "B",  # "CL:0000236",
    "Mal": "Can",  # "CL:0001064",  # malignant cell
    "NK": "NK",  # "CL:0000623",  # natural killer cell
    "Macrophage": "Mac",  # "CL:0000235",
    "CAF": "CAF",  # "CL:0000057",  # fibroblast
    "T.CD8": "CD8",  # "CL:0000625",  # CD8-positive, alpha-beta T cell
    "Endo.": "End",  # "CL:0000115",  # endothelial cell
}  # see https://doi.org/10.1016/j.cell.2018.09.006

TIROSH_CELL_TYPES = {
    0: "CL:0001064",  # malignant cell
    1: "CL:0000084",  # T-cell
    2: "CL:0000236",  # B-cell
    3: "CL:0000235",  # macrophage
    4: "CL:0000115",  # endothelial cell
    5: "CL:0000057",  # fibroblast
    6: "CL:0000623",  # natural killer cell
}  # see https://doi.org/10.1126/science.aad0501