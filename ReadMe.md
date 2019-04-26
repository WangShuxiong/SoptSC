### This is the MATLAB software of SoptSC 
#### An optimization based method for the inference of clustering, cell lineage, pseudotime and cell-cell communication network from single cell RNA-seq data. 

* SoptSC is also avaiable as a R package at: https://mkarikom.github.io/RSoptSC/
* For citation, please refer to 
	> Wang, S., Karikomi, M., MacLean, A.L. and Nie, Q., Cell lineage and communication network inference via optimization for single-cell transcriptomics, Nucleic acids research, 2019.
	> https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkz204/5421812

* Key features of SoptSC
	> 1. Estimation of the number of clusters from the single-cell data.
	> 2. Identify cell subpopulations.
	> 3. Identify marker genes for each cell subpopulation.
	> 4. Infer pseudotime (temporal ordering of cells) in an unsupervised manner: i.e., initial cell is not required.
	> 5. Infer lineage in an unsupervised manner: i.e., initial cluster is not required. 
	> 6. Infer signaling network given a group of Ligand-Receptor pairs and their downstream target genes (to be updated).
	----------------------

- Requires: 
	> MATLAB R2017b or later (where function tSNE is included)
	> or Update the Statistics and Machine Learning Toolbox such that the function tSNE (t-Distributed Stochastic Neighbor Embedding) is included.
	> Please refer to each M file for more detailed descriptions of the corresponding function.

- Main functions:
	> 1.1. runexample.m --- an example run of SoptSC on a specific dataset (see here for QUICKSTART). \
	> 1.2. runexample_signaling.m --- an example run of SoptSC on a a specific dataset for signaling network inference. \
	> 2. SoptSC_cluster.m  --- the M file contains all components of the algorithm for clusters. \
	> 3. Lineage_Ptime.m ---the M file contains all components of the algorithm for pseudotime and lineage inference. \
	> 4. symnmf2 --- non-negative matrix factorization (NMF) tool from: [1] Da Kuang, Chris Ding, Haesun Park, Symmetric Nonnegative Matrix Factorization for Graph Clustering, The 12th SIAM International Conference on Data Mining. (SDM '12), pp. 106--117. \
	> 5. NNDSVD --- SVD-based initialization for NMF from: [2] C. Boutsidis and E. Gallopoulos, SVD-based initialization: A head start for nonnegative matrix factorization, Pattern Recognition, Elsevier. \
	> 6. Plot_cluster.m --- Plot clusters on 2-dimensional space. \
	> 7. Plot_eigengap.m --- Plot eigenvalues of the truncated graph Laplacian to illustrate the estimation of number of clusters. \
	> 8. plot_lineage.m --- plot cluster and pseudotime on the inferred lineage tree. \
	> 9. plot_lineage_marker.m --- plot marker genes on the inferred lineage tree. \
	> 10. plot_marker.m --- plot marker genes on the low dimensional space of cells. \
	> 11. newmap.mat --- colormap used for plot_genes. 
	-----
> **Please feel free to contact us if you have any question: shuxionw 'at' uci 'dot' edu**
