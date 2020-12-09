function SoptSC_obj = creat_soptsc_object(Counts,gene_annotation,gene_ids,cell_annotation,opt)

% Filtering low-quality cells
counts_per_cell = sum(Counts);
counts_per_gene = sum(Counts,2);
aa = zeros(size(Counts));
aa(Counts >0 ) = 1;
genes_per_cell = sum(aa);   % number of genes expressed per cell
cells_per_gene = sum(aa,2); % number of cells expressed per gene

% plot count info
histogram(log(counts_per_cell+1));
xlabel('Counts Per Cell');
ylabel('log(Frequency)');

figure;
histogram(log(counts_per_gene+1));
xlabel('Counts Per Gene');
ylabel('log(Frequency)');

figure;
% plot(sort(genes_per_cell));
histogram(genes_per_cell);
xlabel('Cells');
ylabel('Frequency of Genes Per Cell');

figure;
histogram(cells_per_gene);
xlabel('Genes');
ylabel('Frequency of Cells Per Gene');


% filter out low-quality cells
alpha = 100;
if isfield(opt,'alpha')
    alpha = opt.alpha;
end
cell_idx1 = find(genes_per_cell > alpha);

% Remove low-quality cells based on expression of MT-genes
MT_idx = find(startsWith(gene_annotation, 'MT-'));
display(gene_annotation(MT_idx));
MT_score = sum(Counts(MT_idx,:))./sum(Counts);
 
beta = 0.1;
if isfield(opt,'beta')
    beta = opt.beta;
end
cell_idx2 = find(MT_score < beta);

%% Final selected cells
cell_selected = intersect(cell_idx1,cell_idx2);


SoptSC_obj.data = sparse(log(Counts(:,cell_selected)+1)); % log normalization
display('log-normalization');
SoptSC_obj.gene_annotation = gene_annotation;
SoptSC_obj.gene_ids = gene_ids;
SoptSC_obj.cell_annotation = cell_annotation(cell_selected);

