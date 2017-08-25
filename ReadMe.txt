This is a MATLAB implementation of SoptSC, a computational method for single cell data analysis.

<<<<<<< HEAD
The code file mainly includes:
1.1) SOptSC_cluster.m  --------the M file contains all components of the algorithm for clusters. 
    Please refer to this file for further information.
1.2) SOptSC_pseudotime.m --------the M file contains all components of the algorithm for pseudotime 
    estimation.

2) example.m -------- an example on how to run SOptSC on a specific data.
3) pca.m     -------- PCA algrotihm from the Matlab Toolbox for Dimensionality reduction: 
                      http://homepage.tudelft.nl/19j49
4) symnmf2   -------- Non-negative matrix factorization (NMF) tool from
	[1] Da Kuang, Chris Ding, Haesun Park, Symmetric Nonnegative Matrix Factorization 
            for Graph Clustering, The 12th SIAM International Conference on Data Mining 
            (SDM '12), pp. 106--117.
=======
Requires MATLAB R2015b or later.
>>>>>>> 7058b5954aade1b3820a96d9e2308e45f3ab95ef

This directory includes:

<<<<<<< HEAD
6) Example data is in file Data
7) Plotgenes.m   --------- Plot selected genes specified by users.
8) newmap.mat    --------- colormap used for plot_genes.


9) The computation results is saved in file Results.
=======
   1) example.m -------- an example run of SOptSC on a specific dataset (see here for QUICKSTART).
   2) SOptSC.m  --------  this file contains all the components of the algorithm. Please refer
	                      to its contents for further information.
   3) pca.m     -------- PCA algrotihm from the Matlab Toolbox for Dimensionality reduction: 
                              http://homepage.tudelft.nl/19j49
   3) symnmf2   -------- non-negative matrix factorization (NMF) tool from:
	                      [1] Da Kuang, Chris Ding, Haesun Park, Symmetric Nonnegative Matrix
			      Factorization for Graph Clustering, The 12th SIAM International Conference
			      on Data Mining. (SDM '12), pp. 106--117.
   4) NNDSVD    -------- SVD-based initialization for NMF from:
	                      [2] C. Boutsidis and E. Gallopoulos, SVD-based initialization: A head
    	                      start for nonnegative matrix factorization, Pattern Recognition, Elsevier.
   5) Example data can be found in Data/
   6) The results are saved in Results/
>>>>>>> 7058b5954aade1b3820a96d9e2308e45f3ab95ef


Please refer to example.m for instructions on how to use this code.


Please feel free to send us an email if you have any trouble in running our code. 
The email address for correspondence is shuxionw@uci.edu
