function [P,P_agg] = LR_Interaction(data,allgenes,L,R,TG_act,TG_inh)
% Compute cell-cell signaling network for given inputs of datasets
%   Input:
%   -- data: single cell data matrix with rows associated with genes and
%   columns associated with cells.
%   -- allgenes: cell array of all gene annotations.
%   -- L: cell array of a set of Ligands.
%   -- R: cell array of a set of Receptors corresponding to Ligands.
%   -- TG_act: Target genes of activator.
%   -- TG_inh: Target genes of inhibitors.
%
%   Output:
%   -- P: probability of individual ligand-receptor pair where P{i} 
%   is the probability of cell-cell interactions for the ith 
%   ligand-receptor pair.
%   -- P_agg: aggregate probability of cell-cell intereactions for all
%   given ligand-receptor pairs and their downstream target genes.

No_cells = size(data,2);
No_pairs = length(L);
P = cell(No_pairs,1);
P_agg = zeros(No_cells);
for i = 1:length(L)
    L_i = L{i};
    R_i = R{i};
    
    if nargin==4
        P1 = LR_Prob(data,allgenes,L_i,R_i);
        P1(P1<=1e-6) = 0;
        P{i} = P1;
        P_agg = P_agg + P1;
    elseif nargin==5
        TG_acti = TG_act{i};
        P1 = LRact_Prob(data,allgenes,L_i,R_i,TG_acti);
        P1(P1<=1e-6) = 0;
        P{i} = P1;
        P_agg = P_agg + P1;
    elseif nargin==6
        TG_acti = TG_act{i};
        TG_inhi = TG_inh{i};
        P1 = LRboth_Prob(data,allgenes,L_i,R_i,TG_acti,TG_inhi);
        P1(P1<=1e-6) = 0;
        P{i} = P1;
        P_agg = P_agg + P1;
    end
    
end
P_agg = P_agg./No_pairs;

