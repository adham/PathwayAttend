"""
Apply NMF on pathway data and measure the agreement wit PAM50

Adham Beyki
Deakin UNiversity
odinay@gmail.com
2018-11-05
"""

import numpy as np
import pandas as pd
from sklearn.decomposition import NMF
from sklearn.feature_extraction.text import TfidfTransformer
from sklearn.cluster import KMeans


PATHWAY_DATA ="../data/pathway_data.pkl"
PAM50_PATH = "../data/BRCA.547.PAM50.SigClust.Subtypes.txt"


def main():

    pathway_data = pd.read_pickle(PATHWAY_DATA)

    # prep X
    pathway_subjects = pathway_data['train_subjects'] + pathway_data['test_subjects']
    X = np.vstack([
        pathway_data['X_train'], pathway_data['X_test']
    ])
    tfidf_transformer = TfidfTransformer()
    X = tfidf_transformer.fit_transform(X)

    # fit NMF
    model = NMF(n_components=10, init='random', random_state=1)
    W = model.fit_transform(X)
    H = model.components_

    # cluster the NMF weights
    kmeans = KMeans(n_clusters=5)
    kmeans.fit(W)
    preds = kmeans.labels_

    # find the common subjects and the associated ground truth and prediction
    PAM50_df = pd.read_csv(PAM50_PATH, delimiter='\t')
    PAM50_subjects = PAM50_df['Sample'].tolist()
    subjects = [s for s in pathway_subjects if s in PAM50_subjects]

    idxs1 = [pathway_subjects.index(s) for s in subjects]
    idxs2 = [PAM50_subjects.index(s) for s in subjects]\



    1/0







    1/0



if __name__ == "__main__":
    main()