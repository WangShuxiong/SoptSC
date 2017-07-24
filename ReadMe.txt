This is a matlab implementation of a computational method, SoptSC, for single
cell data analysis.

The code file mainly includes:
1) SOptSC.m  --------the M file contains all components of the algorithm. Please refer
	     to this file for further information.
2) example.m -------- an example on how to run SOptSC on a specific data.
3) pca.m     -------- PCA algrotihm from the Matlab Toolbox for Dimensionality reduction: 
                      http://homepage.tudelft.nl/19j49
4) symnmf2   -------- Non-negative matrix factorization (NMF) tool from
	[1] Da Kuang, Chris Ding, Haesun Park, Symmetric Nonnegative Matrix Factorization 
            for Graph Clustering, The 12th SIAM International Conference on Data Mining 
            (SDM '12), pp. 106--117.

5) NNDSVD    -------- A SVD-based initialization for NMF from	 
	[2] C. Boutsidis and E. Gallopoulos, SVD-based initialization: A head
    	start for nonnegative matrix factorization, Pattern Recognition, Elsevier

6) Example data is in file Data
7) The computation results is saved in file Results.

 
Please refer to example.m for how to use the code.


Please feel free to send emails to us if you have any trouble in running our code. 
The correspondence email is shuxionw@uci.edu
