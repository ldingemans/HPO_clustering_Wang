# HPO_clustering_Wang

These are the scripts use to perform the HPO clustering method as used originally used in https://pubmed.ncbi.nlm.nih.gov/33513338/ and later in other studies.

The workflow is as follows:

1) First use the `get_data.py` script in the Scripts directory to process a phenotips style JSON and get a processable CSV file.
2) Since there are usually quite a lot of HPO terms, which sometimes are very alike, we need to perform some kind of dimensionality reduction technique before we cluster. Otherwise, terms that are very similar would be seen as seperate features and this could lead to bias. Therefore, we try to group the HPO terms in new features. Hereto, we first use `wang_sim_matrx.R` to calculate the Wang similarity score between all HPO terms in the study and plot these, so we can graphically see which terms are alike. Using these information, we can make new grouping variables, as seen in `hpo_groups.xlsx` in the Data directory.
3) The grouped features can then be used in `clustering_in_r.r` to run the actually clustering algorithm and perform the permutation tests.
