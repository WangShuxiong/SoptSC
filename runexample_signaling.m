% Example for using SoptSC for cell-cell and cluster-cluster
% signaling network inference from single cell data.
%
% After running example (data in example is from Joost paper), 
% run this script for signaling network inference for given pathways
% Tgfb, Wnt and Bmp.


%% Setup  
% save cluster labels to sigfolder
resfolder = 'Results';
sigfolder = [resfolder '/Signaling'];
dlmwrite('Data/Joost_cluster_labels.txt', cluster_label)


% Run each section separately or run them all 
%% Tgfb: ligand-receptor pairs and their target genes
Lig = {{'Tgfb1'},{'Tgfb1'},{'Tgfb2'},{'Tgfb2'}};
Rec = {{'Tgfbr1'},{'Tgfbr2'},{'Tgfbr1'},{'Tgfbr2'}};
Target = {{'Zeb2','Smad2','Wnt4','Wnt11','Bmp7','Sox9','Notch1'},...
          {'Zeb2','Smad2','Wnt4','Wnt11','Bmp7','Sox9','Notch1'},...
          {'Zeb2','Smad2','Wnt4','Wnt11','Bmp7','Sox9','Notch1'},...
          {'Zeb2','Smad2','Wnt4','Wnt11','Bmp7','Sox9','Notch1'}};

% Computing cell-cell interaction probability for given ligand-receptor pairs
[Pidv, Pall] = LR_Interaction(data, allgenes, Lig, Rec, Target);


%% Plot the (cell-cell and cluster-cluster) signaling network for Tgfb.
% Set the threshold such that the probability between cells or clusters
% less than the value of threshold is set to be zero. 
% Save the result figurs in folder Results\Signaling
threshold = 0.1;

plot_sig_network(Pidv,Pall,cluster_label,Lig,Rec,threshold,sigfolder)

%% Write P matrices to file
P_write = Pall;
dlmwrite([sigfolder '/P_all_' Lig{1}{1} '.txt'], P_write, 'delimiter','\t')

for i = 1:size(Pidv,1)
    dlmwrite([sigfolder '/P_' Lig{i}{1} '_' Rec{i}{1} '.txt'], Pidv{i}, 'delimiter','\t')
end

close all;




%% Wnt: ligand-receptor pairs and their target genes 
Lig = {{'Wnt3'},{'Wnt4'},{'Wnt5a'},{'Wnt6'},{'Wnt10a'}};
Rec = {{'Fzd1'},{'Fzd1'},{'Fzd1'},{'Fzd1'},{'Fzd1'}};
Target = {{'Ctnnb1','Lgr5','Runx2','Apc','Mmp7','Dkk1','Ccnd1'},...
        {'Ctnnb1','Lgr5','Runx2','Apc','Mmp7','Dkk1','Ccnd1'},...
        {'Ctnnb1','Lgr5','Runx2','Apc','Mmp7','Dkk1','Ccnd1'},...
        {'Ctnnb1','Lgr5','Runx2','Apc','Mmp7','Dkk1','Ccnd1'},...
        {'Ctnnb1','Lgr5','Runx2','Apc','Mmp7','Dkk1','Ccnd1'}};

% Computing cell-cell interaction probability for given ligand-receptor pairs
[Pidv, Pall] = LR_Interaction(data, allgenes, Lig, Rec, Target);


%
threshold = 0.1;

plot_sig_network(Pidv,Pall,cluster_label,Lig,Rec,threshold,sigfolder)

%% Write P matrices to file
P_write = Pall;
dlmwrite([sigfolder '/P_all_' Lig{1}{1} '.txt'], P_write, 'delimiter','\t')

for i = 1:size(Pidv,1)
    dlmwrite([sigfolder '/P_' Lig{i}{1} '_' Rec{i}{1} '.txt'], Pidv{i}, 'delimiter','\t')
end

close all;







%% Bmp: ligand-receptor pairs and their target genes
Lig = {{'Bmp1'},{'Bmp2'},{'Bmp4'},{'Bmp7'}};
Rec = {{'Bmpr2'},{'Bmpr2'},{'Bmpr2'},{'Bmpr2'}};
Target = {{'Crebbp','Fos','Id1','Jun','Runx1','Smad1','Smad5','Sox4','Cdh1'}, ...
          {'Crebbp','Fos','Id1','Jun','Runx1','Smad1','Smad5','Sox4','Cdh1'}, ...
          {'Crebbp','Fos','Id1','Jun','Runx1','Smad1','Smad5','Sox4','Cdh1'}, ...
          {'Crebbp','Fos','Id1','Jun','Runx1','Smad1','Smad5','Sox4','Cdh1'}};

% Computing cell-cell interaction probability for given ligand-receptor pairs
[Pidv, Pall] = LR_Interaction(data, allgenes, Lig, Rec, Target);


%
threshold = 0.03;

plot_sig_network(Pidv,Pall,cluster_label,Lig,Rec,threshold,sigfolder)

%% Write P matrices to file
P_write = Pall;
dlmwrite([sigfolder '/P_all_' Lig{1}{1} '.txt'], P_write, 'delimiter','\t')

for i = 1:size(Pidv,1)
    dlmwrite([sigfolder '/P_' Lig{i}{1} '_' Rec{i}{1} '.txt'], Pidv{i}, 'delimiter','\t')
end

%close all;
