"""
queries the ensemble gene annotations

pathway-gene data uses hgnc annottaion while gene expression uses ensemble. This file
queries mygene.info and builds a dictionary of the equivalent annotations.

Adham Beyki
PRaDA, A@I@ - Deakin University
2018-11-16
"""

import mygene
import pandas as pd
from collections import defaultdict
from joblib import Parallel, delayed


PATHWAY_PATH = "../../data/pathways/hgncToGObp.csv"
df = pd.read_csv(PATHWAY_PATH)
genes = df['SYMBOL'].tolist()

mg = mygene.MyGeneInfo()

# single thread
# for gene in tqdm(genes):
#     result = mg.query(gene, scopes="symbol", fields=["ensembl"], species="human", verbose=False)
#     hgnc_name = gene
#     for hit in result["hits"]:
#         if "ensembl" in hit and "gene" in hit["ensembl"]:
#             hgnc_to_ensembl1[hgnc_name].append(hit["ensembl"]["gene"])


def get_hgnc_to_ensembl(gene):
    result = mg.query(gene, scopes="symbol", fields=["ensembl"], species="human", verbose=False)
    hgnc_name = gene
    ensembl_name = []
    for hit in result["hits"]:
        if "ensembl" in hit and "gene" in hit["ensembl"]:
            ensembl_name.append(hit["ensembl"]["gene"])
    print("{} -> {}".format(hgnc_name, ensembl_name))
    return hgnc_name, ensembl_name

print("working...")
annotation_tuples = Parallel(n_jobs=40, prefer='threads')(
    delayed(get_hgnc_to_ensembl)(gene) for gene in genes
)
hgnc_to_ensembl = dict(annotation_tuples)
pd.to_pickle(hgnc_to_ensembl, "hgnc_to_ensembl.pkl")