This is a MATLAB implementation of SoptSC, a computational method for single cell data analysis. More specifically, SoptSC can

	1. Identify the number of clusters from the input data.
	2. Identify cell subpopulations.
	3. Identify marker genes for each cell subpopulation.
	4. Infer pseudotime (temporal ordering of cells) in an unsupervised manner: 	i.e., initial cell is not required.
	5. Infer lineage in an unsupervised manner: i.e., initial cluster is not 	required. 
	6. Infer signaling network given a group of Ligand-Receptor pairs and
	   their downstream target genes (to be updated).


=======
Requires MATLAB R2015b or later.
=======
Please refer to each M file for more detailed descriptions of the corresponding function.

This directory includes:

   1) example.m -------- an example run of SOptSC on a specific dataset (see here for QUICKSTART).
   
   2) SOptSC_cluster.m  --------the M file contains all components of the algorithm for clusters. 
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
   7) Plotgenes.m   --------- Plot selected genes specified by users.
   8) newmap.mat    --------- colormap used for plot_genes.

   9) Example data can be found in Data
   10) The results are saved in Results


Please refer to example.m for instructions on how to use this code.


Please feel free to send us an email if you have any trouble in running our code. 
The email address for correspondence is shuxionw@uci.edu
