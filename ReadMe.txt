This is a MATLAB implementation of SoptSC, a computational method for single
cell data analysis.

Requires MATLAB R2015b or later.

This directory includes:
   1) SOptSC.m  -------- this file contains all the components of the algorithm. Please refer
	     to its contents for further information.
   2) example.m -------- an example on how to run SOptSC on a specific dataset.
   3) symnmf2   -------- non-negative matrix factorization (NMF) tool from
	[1] Da Kuang, Chris Ding, Haesun Park, Symmetric Nonnegative Matrix Factorization 
            for Graph Clustering, The 12th SIAM International Conference on Data Mining 
            (SDM '12), pp. 106--117.
   4) NNDSVD    -------- SVD-based initialization for NMF from	 
	[2] C. Boutsidis and E. Gallopoulos, SVD-based initialization: A head
    	start for nonnegative matrix factorization, Pattern Recognition, Elsevier
   5) Example data is in file Data
   6) The computation results is saved in file Results.
   
 
Please refer to example.m for instructions on how to use this code.


Please feel free to send us an email if you have any trouble in running our code. 
The email address for correspondence is shuxionw@uci.edu
