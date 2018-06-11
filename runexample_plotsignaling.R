
require(circlize)
source('plotP_functions.R')


## Input cluster labels
clusterlabels = read.table('Data/Joost_cluster_labels.txt')[,1]

## Set results folder
resfolder = 'Results_0612'

## Select ligand (or set of ligands)
Ligand = 'Bmp4_Bmpr2'
#Ligand = '5L_Wnt_all'
#Ligand = 'all_Tgfb1'


## Parameters
sigonly = FALSE   # If true => input is only signaler cells, and P will not be thresholded
plotclusternames = FALSE

subsample = TRUE
threshold = 0.05
nsample = 100

# Path to input data 'P'
P = read.table(paste0(resfolder, '/Signaling/P_',Ligand,'.txt'))
figname = paste0(resfolder, '/Signaling/P_', Ligand,'_',threshold,'.pdf')

## Run
df = plotSig(P, clusterlabels, threshold, sigonly, plotclusternames, subsample, nsample, outfile=figname)


