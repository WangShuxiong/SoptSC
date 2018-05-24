This is a MATLAB implementation of SoptSC, a computational method for single cell data analysis. More specifically, SoptSC can

1. Identify the number of clusters from the input data.
2. Identify cell subpopulations.
3. Identify marker genes for each cell subpopulation.
4. Infer pseudotime (temporal ordering of cells) in an unsupervised manner: i.e., initial cell is not required.
5. Infer lineage in an unsupervised manner: i.e., initial cluster is not required. 
6. Infer signaling network given a group of Ligand-Receptor pairs and their downstream target genes (to be updated).


=======
The current SoptSC manuscript is on bioRxiv: 
Wang, S., MacLean, A.L. and Nie, Q., SoptSC: Similarity matrix optimization for clustering, lineage, and signaling inference.
https://www.biorxiv.org/content/early/2018/05/12/168922

=======
Requires: 
MATLAB R2017b or later (where function tSNE is included)

or Update the Statistics and Machine Learning Toolbox such that the function
tSNE (t-Distributed Stochastic Neighbor Embedding) is included.

=======
Please refer to each M file for more detailed descriptions of the corresponding function.

This directory includes:

   1) Example.m -------- an example run of SOptSC on a specific dataset (see here for QUICKSTART).
   
   2) SoptSC_cluster.m  --------the M file contains all components of the algorithm for clusters. 
    Please refer to this file for further information.
   
   3) Lineage_Ptime.m --------the M file contains all components of the algorithm for pseudotime and lineage inference.
  
   4) pca1.m     -------- PCA algrotihm from the Matlab Toolbox for Dimensionality reduction.
 
   5) symnmf2   -------- non-negative matrix factorization (NMF) tool from:
	                      [1] Da Kuang, Chris Ding, Haesun Park, Symmetric Nonnegative Matrix
			      Factorization for Graph Clustering, The 12th SIAM International Conference
			      on Data Mining. (SDM '12), pp. 106--117.
   6) NNDSVD    -------- SVD-based initialization for NMF from:
	                      [2] C. Boutsidis and E. Gallopoulos, SVD-based initialization: A head
    	                      start for nonnegative matrix factorization, Pattern Recognition, Elsevier.
   7) Plot_cluster.m   --------- Plot clusters on 2-dimensional space.
   8) Plot_eigengap.m  ———————— Plot eigenvalues of the truncated graph Laplacian to illustrate the estimation of number of clusters.
   9) plot_lineage.m ———————— plot cluster and pseudotime on the inferred lineage tree.
   10) plot_lineage_marker.m ———————— plot marker genes on the inferred lineage tree.
   11) plot_marker.m ———————— plot marker genes on the low dimensional space of cells.
   12) newmap.mat    --------- colormap used for plot_genes.
   
   

Please refer to Example.m for instructions on how to use this code.


Please feel free to send us an email if you have any trouble in running our code. 
The email address for correspondence is shuxionw@uci.edu
