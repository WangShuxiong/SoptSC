function [reduction, umap, clusterIdentifiers]=run_umap(varargin)
%%RUN_UMAP reduces data matrices with 3+ parameters down to fewer
%   parameters using the algorithm UMAP (Uniform Manifold Approximation and
%   Projection).
%
%   [reduction,umap] = RUN_UMAP(csv_file_or_data,'NAME1',VALUE1,..., 
%   'NAMEN',VALUEN) 
%   
%
%   RETURN VALUES
%   Invoking run_umap produces 2 return values:
%   1) reduction, the actual data that UMAP reduces from the data 
%     specified by the input argument csv_file_or_data; 
%   2) umap, an instance of the UMAP class made ready for the invoker 
%     to save in a MATLAB file for further use as a template.
%
%
%   REQUIRED INPUT ARGUMENT
%   The argument csv_file_or_data is either 
%   A) a char array identifying a csv text file containing the data 
%       to be reduced. 
%   B) the actual data to be reduced; a numeric matrix.
%
%   If A) then the csv file needs the first line to be parameter names.
%
%   If run_umap is invoked with no arguments it then offers to download
%   example csv files and run one of them.  The programming examples in
%   this documentation use those files.
%
%
%   OPTIONAL INPUT ARGUMENTS
%   Some of these are identical to those in the original Python
%   implementation documented by the inventors in their document "Basic
%   UMAP parameters" which can be retrieved at
%   https://umap-learn.readthedocs.io/en/latest/parameters.html.
%   The optional argument name/value pairs are:
%
%   Name                    Value
%
%   'n_neighbors'           Controls local and global structure as 
%                           does the same parameter in the original 
%                           implementation. 
%                           Default is 30. 
%   
%   'min_dist'              Controls how tightly UMAP is allowed to 
%                           pack points together as does the same 
%                           parameter in the original implementation.
%                           Modifying this value requires the Curve Fitting
%                           Toolbox.
%                           Default is 0.3.
%
%   'metric'                Controls how distance is computed in the
%                           ambient space as does the same parameter in the
%                           original implementation. Accepted values for
%                           metric include 'euclidean', 'cosine',
%                           'cityblock', 'seuclidean', 'correlation',
%                           'jaccard', 'spearman', 'hamming'. These metrics
%                           are described in MATLAB's documentation for
%                           knnsearch.
%                           Default is 'euclidean'.
%
%   'randomize'             true/false.  If false run_umap invokes
%                           MATLAB's "rng default" command to ensure the
%                           same random sequence of numbers between
%                           invocations.
%                           Default is false.
%
%   'template_file'         This identifies a mat file with a saved
%                           instance of the UMAP class that run_umap
%                           previously produced. The instance must be be a
%                           suitable "training set" for the current "test
%                           set" of data supplied by the argument
%                           csv_file_or_data. Template processing
%                           accelerates the UMAP reduction and augments
%                           reproducibility. run_umap prechecks the
%                           suitability of the template's training set for
%                           the test set by checking the parameter names
%                           and the standard deviation distance of the
%                           means for each parameter.
%                           Default is empty ([]...no template).
%
%   'parameter_names'       Cell of char arrays to annotate each parameter
%                           in the data specified by csv_file_or_data. This
%                           is only needed if a template is being used or
%                           saved.
%                           Default is {}.
%                           
%   'verbose'               Accepted values are 'graphic', 'text', or
%                           'none'. If verbose=graphic then the data
%                           displays with probability coloring and contours
%                           as is conventional in flow cytometry analysis.
%                           If method=Java then the display refreshes as
%                           optimize_layout progresses and a progress bar
%                           is shown along with a handy cancel button. If
%                           verbose=text the progress is displayed in the
%                           MATLAB console as textual statements.
%                           Default is 'graphic'.
%                           
%   'method'                Picks 1 of our 6 implementations for UMAP's 
%                           time-consuming optimize_layout phase that does 
%                           stochastic gradient descent.  Accepted
%                           values are 'C++', 'Java', 'C', 'C vectorized', 
%                           'MATLAB' or 'MATLAB Vectorized'. We find our
%                           C++ code the fastest. We build it with clang++.
%                           See the shell scripts build.cmd and build and
%                           C++ source in the subfolder umap/sgdCpp_files.
%                           'C' and 'C vectorized' are produced by MATLAB's
%                           app "C coder".  For open source reasons these
%                           modules are only in the distribution found at
%                           http://cgworkspace.cytogenie.org/GetDown2/demo/umapDistribution.zip
%                           Only Java and C++ support the progress plots
%                           and cancellation options given by argument
%                           verbose=graphic; and the other methods are
%                           artifacts of the programming journey we took to
%                           find the fastest treatment for optimize_layout
%                           that was rapidly developable within the MATLAB
%                           environment.  The C and C vectorized methods
%                           are generated by the MATLAB C coder.
%                           Default is 'Java'.
%
%  'progress_callback'      A MATLAB function handle that run_umap
%                           invokes when method=Java and verbose=graphic. 
%                           The input/output contract for this function is
%                           keepComputing=progress_report(javaObjectOrString).
%                           The javaObjectOrString argument is a char array
%                           before optimize_layout phase starts and then
%                           when optimize_layout is running it is an
%                           instance of the java class
%                           StochasticGradientDescent.java. This instance's
%                           public methods getEpochsToDo, getEpochsDone and
%                           getEmbedding can be used to convey the state of
%                           progress as illustrated in the source code of
%                           run_umap for the function progress_report. If
%                           the function returns keepComputing=false then
%                           run_umap halts the processing.
%                           This argument only takes effect when 
%                           verbose=graphic AND (method=Java OR
%                           method=C++).
%                           Default is the function progress_report
%                           in run_umap.m.
%
%   'ask_to_save_template'  true/false instructs run_umap to ask/not ask
%                           to save a template PROVIDING method='Java',
%                           verbose='graphic', and template_file is empty.
%                           Default is false.
%
%   'label_column'          number identifying the column in the input data
%                           matrix which contains numeric identifiers to
%                           label the data for UMAP supervision mode.                    
%   `                       Default is 0, which indicates no label column.
%
%   'label_file'            the name of a properties file that contains
%                           the label names.  The property name/value
%                           format is identifier=false.
%                           Default is [].
%
%   'n_components'          The dimension of the space into which to embed
%                           the data.
%                           Default is 2.
%
%   'epsilon'               The epsilon parameter required by dbscan if
%                           output argument includes clusters and
%                           n_components > 2.
%                           Default is 0.6.
%
%   'save_template_file'    Fully qualified path of the file to save the
%                           resulting UMAP object as a template (alternate
%                           to output parameter).
%                           Default is [].
%
%   'match_supervisors'     A number indicating how to relabel data points 
%                           in the embedding data if the UMAP reduction 
%                           is guided by a template that in turn is guided 
%                           by supervisory labels.
%                           0 matches supervised and supervising data
%                             groupings by distance of medians. Supervising
%                             groupings are data points in the template's
%                             embedding that have the same supervisory
%                             label. Supervised groupings are DBM clusters
%                             in the final template-guided embedding. The
%                             publication that introduces DBM clustering is
%                             http://cgworkspace.cytogenie.org/GetDown2/demo/dbm.pdf.
%                           1 (default) matches groupings by quadratic 
%                             form dissimilarity.  The publication 
%                             that introduces QF dissimilarity is
%                             https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5818510/.
%                           2 matches supervised DBM clusters by assigning
%                             the label of the closest supervising data 
%                             point to each supervised data point and then 
%                             choosing the most frequent label in the 
%                             cluster.  Closeness is based on euclidean 
%                             distance in the supervised and supervising
%                             embedding data spaces.
%                           3 is similar to 2 except it only uses closeness
%                             to the supervising data point to relabel the
%                             supervised data points without the aid of DBM
%                             clustering.  Thus supervised groupings in the
%                             embedding space may have small fragments
%                             of data points that occur in distant 
%                             islands/clusters of the embedding.
%                           
%   'qf_dissimilarity'      Show QF dissimilarity scores between data
%                           groupings in the supervised and supervising 
%                           embeddings. The showing uses a sortable data 
%                           table as well as a histogram.
%                           Default is false.
%                           run_umap only consults this argument when it
%                           guides a reduction with a supervised
%                           template.
%                           
%   'qf_tree'               Show a dendrogram plot that represents the
%                           relatedness of data groupings in the
%                           supervising and supervised embeddings. The
%                           above documentation for the match_supervisors
%                           argument defines "data groupings". The
%                           publication that introduces the QF tree is
%                           https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6586874/.
%                           This uses phytree from MATLAB's Bioinformatics
%                           Toolbox, hence changing this value to true
%                           requires the Bioinformatics Toolbox.
%                           Default is false.
%                           run_umap only consults this argument when it
%                           guides a reduction with a supervised
%                           template.
%                           
%   joined_transform        true/false for a new transform method to avoid
%                           false positives when applying a template whose
%                           training data differs too much from test set
%                           data. This feature is not part of UMAP's
%                           original Python implementation. 
%                           Currently this not supported when
%                           method=C++.  Support for C++ is coming soon.
%                           Default is false.
%
%   python                  true/false to use UMAP's original
%                           implementation written in Python instead of
%                           this one written in MATLAB, C++ and Java.  The
%                           Python implementation is from Leland McInnes,
%                           John Healy, and James Melville.
%                           If true then certain arguments are ignored:
%                           joined_transform, method, verbose, 
%                           and progress_callback.
%                           Default is false.  
%                   
%
%   EXAMPLES 
%   Note these examples assume your current MATLAB folder is where
%   run_umap.m is stored.
%
%   1.  Download the example csv files and run sample10k.csv.
%
%       run_umap
%
%   2.  Reduce parameters for sample30k.csv and save as template.
%
%       [~, umap]=run_umap('sample30k.csv');
%       save('myTemplate30k.umap.mat', 'umap');
%
%   3.  Reduce parameters for sample130k.csv using prior template.
%
%       run_umap('sample130k.csv', 'template_file', 'myTemplate30k.umap.mat');
%
%   4.  Reduce parameters for sampleBalbcLabeled55k.csv supervised by
%       labels produced by Epp and save as a template. Epp is a more
%       conservative clustering technique described at
%       https://www.nature.com/articles/s42003-019-0467-6. Epp stands for
%       "Exhaustive Projection Pursuit".  By clustering exhaustively in 2
%       dimension pairs, this technique steers more carefully away from the
%       curse of dimensionality than does UMAP or t-SNE.
%
%       To use Epp you can download AutoGate from cytogenie.org which
%       contains tutorials on using Epp.
%
%       [~, umap]=run_umap('sampleBalbcLabeled55k.csv', 'label_column', 11, 'label_file', 'balbcLabels.properties');
%       save('myTemplateBalbcEpp55k.umap.mat', 'umap');
%
%   5.  Reduce parameters for sampleRag148k.csv using template that is
%       supervised by Epp.  This takes the clusters created by Epp on the
%       lymphocytes of a normal mouse strain (BALB/c) and applies them via
%       a template to a mouse strain (RAG) that has neither T cells nor B
%       cells.
%
%       run_umap('sampleRag148k.csv', 'template_file', 'myTemplateBalbcEpp55k.umap.mat');
%
%   6.  Reduce parameters for sample30k.csv and return cluster identifiers
%       using density-based merging described at 
%       http://cgworkspace.cytogenie.org/GetDown2/demo/dbm.pdf.
%
%       [~,~, clusterIds]=run_umap('sample30k.csv');
%
%   7.  Repeat sample 2 but for 3D output and return cluster identifiers
%       and save the result as 3D template.
%
%       [~, umap, clusterIds]=run_umap('sample30k.csv', 'n_components', 3);
%       save('myTemplate30k3D.umap.mat', 'umap');
%
%   8.  Repeat example 3 in 3D.
%
%       run_umap('sample130k.csv', 'template_file', 'myTemplate30k3D.umap.mat', 'n_components', 3);
%
%   9.  Reduce parameters and save template for sampleRagLabeled60k.csv
%       using labels produced by an expert biologist drawing manual gate 
%       sequences on lymphocyte data taken from a RAG mouse strain which 
%       has no T cells or B cells.
%
%       [~, umap]=run_umap('sampleRagLabeled60k.csv', 'label_column', 11, 'label_file', 'ragLabels.properties');
%       save('myTemplateRag60k.umap.mat', 'umap');
%
%   10. Reduce parameters for lymphocyte data taken from a BALB/c mouse
%       strain using template created in example 9.  This takes the
%       clusters created on the lymphocyte data of a knockout mouse strain
%       (RAG) with no B cells or T cells and applies them to a normal mouse
%       strain (BALB/c) which has both cell types.  This illustrates logic
%       to prevent false positives for data not seen when training/creating
%       supervised templates.  Choose to re-supervise to see effect.
%
%       run_umap('sample30k.csv', 'template_file', 'myTemplateRag60k.umap.mat');
%
%   11. Repeat example 10 but use joined_transform.
%
%       run_umap('sample30k.csv', 'template_file', 'myTemplateRag60k.umap.mat', 'joined_transform', true);
%
%   12. Run example 5 again showing QF tree and QF dissimilarity plots.
%
%       run_umap('sampleRag148k.csv', 'template_file', 'myTemplateBalbcEpp55k.umap.mat', 'qf_tree', true, 'qf_dissimilarity', true);
%
%   13. Compare our implementation with C++ method to the original Python 
%       implementation by repeating example 2 as follows.
%       
%       [~, umap]=run_umap('sample30k.csv', 'method', 'C++');
%       [~, umap]=run_umap('sample30k.csv', 'python', true);
%
%   14. Compare our implementation with C++ method to the original Python 
%       implementation by repeating example 4 as follows.
%
%       NOTE: The Python example below works for Python code retrieved in
%       September 2019. But as of February 10, 2020, there seems to be an
%       incompatibility between the Python packages "umap-learn" and
%       "numba". You may find that the numba package will throw an error
%       while doing supervised UMAP. One temporary workaround is to add a #
%       symbol at the start of line 500 of (Python
%       folder)\Lib\site-packages\umap\umap_.py, which will disable numba
%       for the function "fast_intersection" (though this will slow down
%       the Python code, making our time comparison less competitive).
%
%       run_umap('sampleBalbcLabeled55k.csv', 'label_column', 11, 'label_file', 'balbcLabels.properties', 'method', 'C++');
%       [~, umap]=run_umap('sampleBalbcLabeled55k.csv', 'label_column', 11, 'label_file', 'balbcLabels.properties', 'python', true);
%       save('pyTemplateBalbcEpp55k.umap.mat', 'umap');
%
%   15. Compare our implementation with C++ method to the original Python 
%       implementation by repeating example 5 as follows.
%
%       run_umap('sampleRag148k.csv', 'template_file', 'myTemplateBalbcEpp55k.umap.mat', 'method', 'C++');
%       run_umap('sampleRag148k.csv', 'template_file', 'pyTemplateBalbcEpp55k.umap.mat');
%
%
%   NOTE that you can do supervised UMAP and templates with n_components
%       ...but we have not had time to update the 3D GUI to show where the 
%       supervised regions fall.
%
%
%   REQUIRED PATHS
%   This distribution has 2 folders:  umap and util.  
%   You must set paths to these folders plus the java inside of umap.jar.
%   Assume you have put these 2 folders under /Users/Stephen.
%   The commands that MATLAB requires would be:
%
%   addpath /Users/Stephen/umap
%   addpath /Users/Stephen/util
%   javaaddpath('/Users/Stephen/umap/umap.jar');
%
%
%   ALGORITHMS
%   UMAP is the invention of Leland McInnes, John Healy and James Melville
%   at Canada's Tutte Institute for Mathematics and Computing.  See
%   https://umap-learn.readthedocs.io/en/latest/.
%
%   AUTHORSHIP
%   Primary Developer+math lead: Connor Meehan <connor.gw.meehan@gmail.com>
%   Secondary Developer:  Stephen Meehan <swmeehan@stanford.edu> 
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University. 
%   License: BSD 3 clause
%
%   IMPLEMENTATION NOTES
%   This is a total rewrite of the the original Python implementation
%   from Leland McInnes, John Healy and James Melville.  
%   This implementation is written in MATLAB, C++ and Java. 
%   The source is distributed openly on MathWorks File Exchange.
%   This implementation follows a very similar structure to the Python
%   implementation, and many of the function descriptions are nearly
%   identical. Leland McInnes has looked over it and considered it "a
%   fairly faithful direct translation of the original Python code (except
%   for the nearest neighbor search)". If you have UMAP's Python
%   implementation you can check how faithful and fast this
%   re-implementation is by using the argument python=true. When
%   python=false and method=C++ we observe superior performance on our Mac
%   and Windows laptops in most cases. For the cases of template-guided and
%   supervised parameter reduction the performance is significantly faster
%   than the Python implementation regardless of data size.  For the cases
%   of basic parameter reduction, however, performance is slower if the
%   data has less than 40,000 rows, but increasingly faster at greater
%   sizes. If you wish to have a simple user GUI to run all of these UMAP
%   features, download AutoGate at CytoGenie.org.
%   
%
clusterIdentifiers=[];
initJava;
pth=fileparts(mfilename('fullpath'));
pPth=fileparts(pth);
utilPath=fullfile(pPth, 'util');
addpath(utilPath);
globals=BasicMap.Global;
try
    props=fullfile(File.Home, 'run_umap.mat');
    globals.load(props);
catch ex
end
reduction=[];
umap=[];
p=parseArguments();
parse(p,varargin{:});
args=p.Results;   
plotting=strcmpi(args.verbose, 'graphic');
csv_file_or_data=args.csv_file_or_data;
save_template_file = args.save_template_file;
curAxes=[];
if plotting
    xLabel=[];
    yLabel=[];
    zLabel=[];
    fig=figure('name', 'Running UMAP ...');
    if args.qf_tree 
        movegui(fig, 'south')
    elseif args.qf_dissimilarity
        movegui(fig, 'center')
    else
        movegui(fig, 'onscreen');
    end
    curAxes=gca;
    if isempty(csv_file_or_data)
        if askYesOrNo(Html.WrapHr(...
                ['Should run_umap.m download example csv files<br>',...
                'from the Herzenberg Lab @ Stanford University<br><br>', ...
                '.. and then run one of them?']))
            csv_file_or_data=downloadCsv;
        end
        if isempty(csv_file_or_data)
            if plotting
                delete(fig);
            end
            globals.save;
            return;
        end
        if ~askYesOrNo(Html.Wrap([...
                'Test csv files have been downloaded:<ol>'...
                ' <li>sample10k<li>sample30k<li>sampleBalbcLabeled55k'...
                '<li>sample130k<li>sampleRag148k.csv<li>sampleRag55k.csv'...
                '</ol><br><center>Run UMAP on <b>sample10k</b> now?<hr></center>']))
            if plotting
                delete(fig);
            end
            globals.save;
            return;
        end
    end
end
if nargout>=3 && args.n_components>2
    % check for presence of DBSCAN
    if ~Density.HasDbScan(false)
        if plotting
            if ~askYesOrNo(Html.WrapHr(['DBSCAN for clustering in 3+D is '...
                    '<br>not downloaded ...Continue?']))
                delete(fig);
                globals.save;
                return;
            end
        end
        dispNoDbScan;
    end
end
if ischar(csv_file_or_data)
    if ~exist(csv_file_or_data, 'file')
        showMsg(Html.WrapHr(['The text file "<b>' csv_file_or_data ...
            '</b>"<br><font color="red"><i>can not be found !!</i></font>']));
        if plotting
            delete(fig);
        end
        globals.save;
        return;
    end
    
    t=readtable(csv_file_or_data, 'ReadVariableNames', true);
    inData=table2array(t);
    parameter_names=File.CsvNames(t);
else
    inData=csv_file_or_data;
    parameter_names=args.parameter_names;
end
newSubsetIdxs=[];
template_file=args.template_file;
if ~isempty(template_file)
    if ischar(template_file)
        if ~exist(template_file, 'file')
            showMsg(Html.WrapHr(['The template file "<b>' template_file ...
                '</b>"<br><font color="red"><i>can not be found !!</i></font>']));
            if plotting
                delete(fig);
            end
            globals.save;
            return;
        end
        if length(parameter_names)~=size(inData, 2)
            showMsg(Html.WrapHr(sprintf(['<b>Can not create '...
                'template</b> ...<br>'...
                '%d parameter_names... but data has %d parameters?'], ...
                length(parameter_names), size(inData,2))));
            if plotting
                delete(fig)
            end
            globals.save;
            return;
        end
        [umap, ~, canLoad, reOrgData, paramIdxs]=...
            Template.Get(inData, parameter_names, ...
            template_file, 3);
        if ~isempty(reOrgData)
            % column label order differed
            inData=reOrgData;
            if ~isempty(parameter_names) && ~isempty(paramIdxs)
                parameter_names=parameter_names(paramIdxs);
            end
        end
    elseif isa(template_file, 'UMAP')
        umap = template_file;
        canLoad = true;
    end
    if isempty(umap)
        if ~canLoad
        if plotting 
            showMsg(Html.WrapHr(['No template data found in <br>"<b>', ...
                template_file '</b>"']));
        else
            disp(['No template data found in ' template_file]);
        end
        end
        if plotting
            delete(fig);
        end
        globals.save;
        return;
    else
        args.n_components=umap.n_components;
        if ~isempty(umap.supervisors)
            %Connor's NEW joined_transform immunizes reduction from
            %false positives if items in the test set are TOO different
            %from the training set
            
            if ~args.joined_transform
                [percNewSubsets, unknownIdxs]=...
                    Template.CheckForUntrainedFalsePositives(umap, inData);
                if percNewSubsets>10
                    [choice, cancelled]=Template.Ask(percNewSubsets);
                    if cancelled
                        if plotting
                            delete(fig);
                        end
                        globals.save;
                        return;
                    end
                    if choice==2
                        umap.clearLimits;
                        newSubsetIdxs=unknownIdxs;
                        template_file=[];
                    end
                end
            else
                umap.supervisors.veryHighClusterDetail=false;
            end
        end
    end
else
    umap = UMAP;
    umap.dimNames=parameter_names;
end
umap.metric=args.metric;
umap.n_epochs=args.n_epochs;
umap.n_neighbors=args.n_neighbors;
umap.min_dist=args.min_dist;
umap.n_components=args.n_components;
if strcmpi('Java', args.method)
    if ~initJava
        args.method='C';
        showMsg(Html.WrapHr('Could not load umap.jar for Java method'), ...
            'Problem with JAVA...', 'south west', false, false);
    end
end
        
method=umap.setMethod(args.method);
umap.verbose=~strcmpi(args.verbose, 'none');
umap.random_state=~args.randomize;
%umap.negative_sample_rate=30;
tick=tic;
[R,C]=size(inData);

labelMap=[];
[nRows, nCols]=size(inData);
sCols=num2str(size(inData,2));
nParams=length(parameter_names);
good=nParams==0||nParams==nCols || (args.label_column>0 &&...
    (nParams==nCols-1 || nParams==nCols));
if ~good
    if args.label_column>0
        preAmble=sprintf(['# data columns=%d, # parameter_names=%d '...
            'since label_column=%d <br># parameter_names must be '...
            '%d or %d '],  nCols, nParams, args.label_column, ...
            nCols, nCols-1);
    else
        preAmble=sprintf(['# of data columns(%s) must equal '...
            '# of parameter_names(%d)'], sCols, nParams);
    end
    msg(Html.WrapHr(preAmble));
    assert(nParams==0||nParams==nCols || (args.label_column>0 &&...
        (nParams==nCols-1 || nParams==nCols)), preAmble);    
end
if ~isempty(newSubsetIdxs)
    hasLabels=true;
    labelCols=0;
    [labels, labelMap]=resupervise(umap, inData, newSubsetIdxs);
    nLabels=length(unique(labels));
elseif args.label_column>0
    if ~isempty(template_file)
        showMsg(Html.WrapHr(['Can not do supervised mode <br>'...
            'AND use prior template at<br>the same time!']));
        globals.save;
        return;
    end
    hasLabels=true;
    labelCols=1;
    good=args.label_column>0 && args.label_column<=nCols;
    if ~good
        msg(Html.WrapHr(['The input data has ' sCols 'columns ...<br>'...
            'THUS the label_column must be >=1 and <= ' sCols]));
        assert(args.label_column>0 && args.label_column<=nCols, [...
            'label_column must be >=1 and <= ' sCols]);
        
    end
    labels=inData(:, args.label_column);
    nLabels=length(unique(labels));
    if nLabels > .5*nRows
        preAmble=sprintf(...
            'WARNING:  %d is a LOT of distinct labels for a %dx%d matrix!',...
            nLabels, nRows, nCols);
        msg(preAmble);
        warning(preAmble);
    end
    inData(:, args.label_column)=[];
    if args.label_column<=nParams
        parameter_names(args.label_column)=[];
    end
    umap.dimNames=parameter_names;
    nLabels=length(unique(labels));
    if exist(args.label_file, 'file')
       map=java.util.Properties;
       try
           map.load(java.io.FileInputStream(args.label_file));
       catch ex
           showMsg(['Can not load ' args.label_file]);
           delete(fig);
           globals.save;
           return;
       end
       labelMap=map;
   end
else
    hasLabels=false;
    labelCols=0;
    nLabels=0;
end
if any(isnan(inData(:)))
    if plotting
        if isequal('Yes', questdlg({...
                'Data matrix has NAN values which',...
                'which cause odd effects on UMAP!','', ...
                'Try to remove nan values?'}))
            allNanColumns=all(isnan(inData));
            if any(allNanColumns)
                inData=inData(:, ~allNanColumns);
            end
            allNanRows=all(isnan(inData'));
            if any(allNanRows)
                inData=inData(~allNanRows,:);
            end
            [R,C]=size(inData);
        end
        if any(isnan(inData(:)))
            showMsg(Html.WrapHr(['Sorry...<br>can not proceed<br>'...
                '<br>NAN values exist... SIGH!']));
            globals.save;
            return;
        end
    else
        error('Can not proceed with NAN values');
    end
end

info=[String.encodeInteger(R) 'x' String.encodeInteger(C-labelCols)];
if ischar(csv_file_or_data)
    [~, fileName]=fileparts(csv_file_or_data);
    info=['UMAP on ' fileName ', ' info];
else
    info=['[UMAP on ' info];
end
if args.python
    info=[info ', Python'];
else
    info=[info ', ' method];
end
if plotting
    set(fig, 'NumberTitle', 'off', 'name', info );
    drawnow;
end
pause(.01);
info2=['(optimize\_layout method=' method ')'];
if strcmpi(method, 'C++')
    if ~StochasticGradientDescent.IsAvailable
        if ~askYesOrNo(Html.Wrap(...
                ['This C++ executable is missing or corrupt:'...
            '<br>"<b>' StochasticGradientDescent.GetCmd '</b>"'...
            '<br><br>Maybe try rebuilding by changing clang++ '...
            'to g++ in the build scripts in the same folder...<br>'...
            '<br><center>Try <b>method=Java</b> instead?</center><hr>']))
            return;
        end
        method='Java';
    end
end
if strcmpi(method, 'Java') || strcmpi(method, 'C++')
    if plotting
        umap.progress_callback=args.progress_callback;
        set(fig, 'NumberTitle', 'off', 'name', info);
        try
            nTh=edu.stanford.facs.swing.StochasticGradientDescent.EPOCH_REPORTS+3;
            figure(fig);
            if args.qf_tree
                puLocation='north++';
            else
                puLocation='south++';
            end
            path=BasicMap.Path;
            pu=PopUp(Html.WrapHr(sprintf(['Using UMAP to reduce '...
                ' <b>%d</b> parameters down to ' ...
                num2str(args.n_components) '...'], C-labelCols)), ...
                puLocation, 'Reducing parameters...', false, true, ...
                fullfile(path, 'genieSearch.png'));
            pu.initProgress(nTh);
            pu.pb.setStringPainted(true);
            pu.setTimeSpentTic;
            drawnow;
        catch ex
            args.method='C';
            method=umap.setMethod(args.method);
            showMsg(Html.WrapHr(['Could not load umap.jar for Java method'...
                '<br><br>Switching optimize_layout method to "C" ']), ...
                'Problem with JAVA...', 'south west', false, false);
        end
    end
end
tc=tic;
if plotting
    if ispc
        left=.21;
        width=.61;
        height=.145;
        lbl=annotation(fig, 'textbox','String', {['\color{blue}Running '...
            info], ['\fontsize{9}' info2]}, 'units', 'normalized', ...
            'position', [left .4 width height], 'fontSize', 11, ...
            'HorizontalAlignment', 'center');
    else
        left=.21;
        width=.61;
        height=.131;
        lbl=annotation(fig, 'textbox','String', {['\color{blue}Running '...
            info], ['\fontsize{11}' info2]}, 'units', 'normalized', ...
            'position', [left .4 width height], 'fontSize', 13, ...
            'HorizontalAlignment', 'center');
    end
    updatePlot;
end
paramAnnotation=[];
if umap.verbose
    txt=sprintf(['n\\_neighbors=\\color{blue}%d\\color{black}, '...
        'min\\_dist=\\color{blue}%s\\color{black}, '...
        'metric=\\color{blue}%s\\color{black},'...
        'randomize=\\color{blue}%d\\color{black}, '...
        'labels=\\color{blue}%d'], ...
        umap.n_neighbors, num2str(umap.min_dist), umap.metric,...
        ~umap.random_state, nLabels); 
    disp(txt);
    if plotting
        paramAnnotation=annotation(fig, 'textbox','String', txt,...
            'units', 'normalized', 'position', [.03 .94 .92 .05],...
            'fontSize', 9, 'HorizontalAlignment', 'center');
        drawnow;
    end
end
if ~isempty(template_file)
    if ~isempty(umap.pythonTemplate)
        args.python=true;
    end
    if ~args.python 
        if ~args.joined_transform
            reduction=umap.transform(inData);
        else
            reduction=umap.transform2(inData);
        end
    else
        if isempty(umap.pythonTemplate) || ~exist(umap.pythonTemplate, 'file')
            reduction=[];
            msg(Html.WrapHr('Python template file not found'),  8,...
                'south', 'Error...', 'error.png');
        else
            inFile=[tempname '.csv'];
            reduction=UmapPython.Go(inData,inFile, [], ...
                lower(umap.metric), umap.n_neighbors, ...
                umap.min_dist, umap.n_components, [], ...
                umap.pythonTemplate);
        end
    end
else
    if ~args.python
        if ~hasLabels
            reduction = umap.fit_transform(inData);
        else
            reduction = umap.fit_transform(inData, labels);
            if ~isempty(reduction)
                if ~isempty(labelMap)
                    umap.setSupervisors(labels, labelMap, curAxes);
                end
            end
        end
    else
        inFile=[tempname '.csv'];
        if ~hasLabels
            labels=[];
        end
        reduction=UmapPython.Go(inData,inFile, [], ...
            lower(umap.metric), umap.n_neighbors, ...
            umap.min_dist, umap.n_components, labels);
        pythonTemplate=fullfile(fileparts(inFile), ...
            [UmapPython.PYTHON_TEMPLATE '.umap']);
        if ~isempty(reduction)
            umap.embedding=reduction;
            umap.raw_data=inData;
            if ~isempty(labelMap)
                umap.setSupervisors(labels, labelMap, curAxes);
            end
        end
    end
end
if ~isempty(paramAnnotation)
    set(paramAnnotation, 'visible', 'on');
end
if ~isempty(reduction)
    if plotting
        figure(fig);
        delete(lbl);
        if strcmpi(method, 'Java') 
            pu.pb.setString('All done');
        end
        updatePlot(reduction, true)
        annotation(fig, 'textbox', 'String', ['Compute time=\color{blue}' ...
            String.MinutesSeconds(toc(tick))],'units', 'normalized', ...
            'position', [.65 .01 .33 .05], 'fontSize', 9)
        if isempty(template_file) && args.ask_to_save_template
            if isequal('Yes', questdlg({'Save this UMAP reduction', ...
                    'as template to accelerate reduction', ...
                    'for compatible other data sets?'}))
                if length(parameter_names)~=size(inData, 2)
                    showMsg(Html.WrapHr(sprintf(['<b>Can not create '...
                        'template</b> ...<br>'...
                        '%d parameter_names ...but data has %d parameters?'], ...
                        length(parameter_names), size(inData,2))));
                else
                    umap.prepareForTemplate(curAxes);
                    if ischar(csv_file_or_data)
                        Template.Save(umap, csv_file_or_data);
                    else
                        Template.Save(umap, fullfile(pwd, 'template.csv'));
                    end
                end
            end
        end
    end
    String.MinutesSeconds(toc(tick))
else
    msg('Parameter reduction was cancelled or not done');
end
if plotting
    if strcmpi(method, 'Java') || strcmpi(method, 'C++')
        pu.stop;
        pu.dlg.dispose;
    end
end
if (nargout>1 || ~isempty(save_template_file)) && ~isempty(reduction)
    if plotting
        umap.prepareForTemplate(curAxes);
    else
        umap.prepareForTemplate;
    end
    
    if  ~isempty(save_template_file)
        save(save_template_file, 'umap');
    end
    if args.python
        if  ~isempty(save_template_file)
            [f1, f2]=fileparts(save_template_file);
        elseif ischar(csv_file_or_data)
            [f1, f2]=fileparts(csv_file_or_data);
            if isempty(f1)
                f1=pwd;
            end
            f2=[f2 '.umap'];
        else
            f1=[];
        end
        if ~isempty(f1)
            umap.pythonTemplate=fullfile(f1, [f2 '.python']);
            movefile(pythonTemplate, umap.pythonTemplate, 'f');
        end
    end
    if nargout>2
        clusterIdentifiers=doClusters(reduction, plotting);
        if isempty(clusterIdentifiers)
            dispNoDbScan;
        end
    end
end     
globals.save;

    function clusterIds=doClusters(data, show)
        if isempty(data)
            clusterIds=[];
        else
            [numClusters, clusterIds]=findClusters(data);
            if show
                if ~isempty(clusterIds)
                    ClusterPlots.Go([], data, clusterIds, xLabel, ...
                        yLabel, zLabel);
                end
            end
        end
    end

    function updatePlot(data, lastCall)
        if nargin>0
            if isempty(xLabel)
                dimInfo=sprintf('  %dD\\rightarrow%dD', C-labelCols, ...
                    args.n_components);
                xLabel=['UMAP-X' dimInfo];
                yLabel=['UMAP-Y' dimInfo];
                if args.n_components>2
                    zLabel=['UMAP-Z' dimInfo];
                end
            end
            if args.n_components>2
                nD=size(data, 2);
                assert(nD==args.n_components);
                if args.frequencyDensity3D
                    Gui.PlotDensity3D(curAxes, data, 64, 'iso',...
                        xLabel, yLabel, zLabel);
                else
                    Gui.PlotNeighDist3D(curAxes, data, ...
                        args.n_neighbors);
                end
                if args.n_components>3
                    title(curAxes, ['NOTE:  Only 3 of \color{red}' ...
                        num2str(args.n_components) ...
                        ' dimensions being shown...']);
                end
            else
                if hasLabels
                    if nargin>1 && lastCall
                        ProbabilityDensity2.DrawLabeled(curAxes, data,...
                            labels, labelMap, true, true, [], [], -0.01,...
                            0.061, true,[], [], [], .05);
                        if ~isempty(umap.supervisors)
                            umap.supervisors.prepareForTemplate;
                            if args.qf_tree
                                umap.supervisors.inputData=umap.raw_data;
                                umap.supervisors.qfTreeSupervisors(true, ...
                                    [], 'UMAP training set');
                            end
                        end
                    else
                        ProbabilityDensity2.DrawLabeled(curAxes, data, ...
                            labels, labelMap);
                    end
                    if ~isempty(paramAnnotation)
                        set(paramAnnotation, 'visible', 'off');
                    end
                else
                    if isprop(umap, 'supervisors') ...
                            && ~isempty(umap.supervisors)
                        if nargin>1 && lastCall
                            if isempty(umap.supervisors.embedding)
                                umap.supervisors.embedding=umap.embedding;
                            end
                            umap.supervisors.computeAndMatchClusters(...
                                data, args.match_supervisors, pu);
                            [lbls, lblMap]=...
                                umap.supervisors.supervise(data, true);
                            
                            if args.qf_tree || args.qf_dissimilarity
                                umap.supervisors.inputData=umap.raw_data;
                                cascading={};
                                if args.qf_tree
                                    [~,qft]=umap.supervisors.qfTreeSupervisors(false, pu);
                                    if ~isempty(qft.fig)
                                        cascading{end+1}=qft.fig;
                                        Gui.CascadeFigs(cascading, false, true, 15, 2, true, false, fig);
                                        figure(fig);
                                        [~,qft]=umap.supervisors.qfTreeSupervisees(data, inData, false, pu);
                                        cascading{end+1}=qft.fig;
                                        Gui.CascadeFigs(cascading, false, true, 40, 2, true, false, fig);
                                    end
                                end
                                if args.qf_dissimilarity
                                    figure(fig);
                                    [~,qft]=umap.supervisors.qfDissimilarity(data, inData, false, pu);
                                    cascading{end+1}=qft.fig;
                                    Gui.CascadeFigs(cascading, false, true, 40, 2, true, false, fig);
                                    Gui.Locate(qft.qHistFig, qft.fig, 'south+'); 
                                    figure(qft.qHistFig);
                                    if ismac && 2==size(get(0, 'MonitorPositions'),1)
                                        drawnow;
                                        Gui.Locate(qft.qHistFig, qft.fig,...
                                            'south+');
                                    end
                                    figure(fig);
                                end
                            end
                            drawnow;
                            ProbabilityDensity2.DrawLabeled(curAxes, data,...
                                lbls, lblMap, true, true, [], [], -0.01,...
                                0.061, true,[]);
                        else
                            umap.supervisors.computeAndMatchClusters(data, 0, pu);
                            [lbls, lblMap]=...
                                umap.supervisors.supervise(data);
                            ProbabilityDensity2.DrawLabeled(curAxes,...
                                data, lbls, lblMap);
                        end
                        if ~isempty(paramAnnotation)
                            set(paramAnnotation, 'visible', 'off');
                        end
                    else
                        if nargin>1 && lastCall
                            ProbabilityDensity2.Draw(curAxes, data, ...
                                true, true, true, .05);
                        else
                            ProbabilityDensity2.Draw(curAxes, data);
                        end
                    end
                end
            end
        end
        if isprop(umap, 'xLimit')
            if ~isempty(umap.xLimit)
                xlim(curAxes, umap.xLimit);
                ylim(curAxes, umap.yLimit);
                if nargin>0
                    Gui.StretchLims(curAxes, data, .04)
                end
            end
        end
        xlabel(curAxes, xLabel);
        ylabel(curAxes, yLabel)
        if args.n_components>2
            zlabel(curAxes, zLabel);
        end
        grid(curAxes, 'on')
        set(curAxes, 'plotboxaspectratio', [1 1 1])
        if nargin>1 && lastCall
            if args.n_components==2
                if isprop(umap, 'supervisors') ...
                        && ~isempty(umap.supervisors)
                    if ~hasLabels
                        umap.supervisors.drawClusterBorders(curAxes);
                    end
                end
            end
        end
        drawnow;
    end

    function keepComputing=progress_report(embeddingOrStatusString)
        keepComputing=~pu.cancelled;
        if ischar(embeddingOrStatusString)
            if ~isequal(embeddingOrStatusString, ...
                    StochasticGradientDescent.FINDING_ISLANDS) 
                pu.pb.setValue(pu.pb.getValue+1);
            end
            pu.pb.setString(embeddingOrStatusString);
            pu.pack;
            pu.showTimeSpent;
            return;
        end 
        done=embeddingOrStatusString.getEpochsDone-1;
        toDo=embeddingOrStatusString.getEpochsToDo;
        pu.pb.setValue(3+(pu.pb.getMaximum*(done/toDo)));
        pu.pb.setMaximum(toDo);
        pu.pb.setString(sprintf('%d/%d epochs done', done, toDo));
        if isvalid(lbl)
            delete(lbl);
        end
        updatePlot(embeddingOrStatusString.getEmbedding);
        pu.showTimeSpent;
    end

    
    function file=downloadCsv
        if ispc
            prompt='Specify name & folder for saving zip file download';
        else
            prompt=Html.WrapHr(...
                ['Please specify the name and folder for the'...
                '<br>zip file being downloaded'...
                '<br>(which will be unzipped after)']);
        end
        [fldr, file]=FileBasics.UiPut(pwd, 'samplesFromHerzenbergLab.zip', ...
            prompt);
        if isnumeric(fldr)
            file=[];
            return;
        end
        pu=PopUp('Downloading & unzipping samples');        
        zipFile=fullfile(fldr, file);
        websave(zipFile, ...
            'http://cgworkspace.cytogenie.org/GetDown2/demo/samples.zip');
        unzip(zipFile);
        pu.close;
        file=fullfile(fldr, 'sample10k.csv');
    end

    function ok=validateCallback(x)
        ok=isequal('function_handle', class(x));
    end

    function ok=validateParameterNames(x)
        ok=false;
        if iscell(x)
            N=length(x);
            if N>0
                for i=1:N
                    if ~ischar(x{i})
                        ok=false;
                        return;
                    end
                end
                ok=true;
            end
        end
    end

    function [numClusters, clusters]=findClusters(data)
        clusters=[];
        epsilon=args.epsilon;
        if epsilon<1
            epsilon=1;
        end
        neighbors=args.n_neighbors;
        if args.n_components>2
            pu2=PopUp('Finding clusters with dbscan', 'north', ...
                'Clustering...', false);
            try
                clusters=dbscan(data, epsilon, neighbors);
            catch ex
                try
                    % see if DBSCAN is available
                    DBSCAN(randi(100, 50, 3), .5, 15);
                    %YES
                    pu2.setText(Html.WrapHr(['MATLAB''s dbscan not '...
                        'available ... <br>using DBSCAN from'...
                        'MathWorks File Exchange.<br><br><i>Note that'...
                        ' DBSCAN is quite slow!</i>']));
                    clusters=DBSCAN(data, epsilon, neighbors);
                catch ex
                    disp('No dbscan ... before r2019a');
                    if plotting
                        Density.DownloadDbScan;
                    end
                end
            end
            pu2.close;
            numClusters=sum(unique(clusters)>0);
        else
            [numClusters,clusters]=Density.ClusterVeryHigh(data);
        end
    end


    function p=parseArguments(varargin)
        p = inputParser;
        defaultMetric = 'euclidean';
        expectedMetric = {'precomputed', 'euclidean', 'l2', 'manhattan', 'l1',...
        'taxicab', 'cityblock', 'seuclidean', 'standardised_euclidean',...
        'chebychev', 'linfinity', 'linfty', 'linf', 'minkowski',...
        'mahalanobis', 'cosine', 'correlation', 'hamming', 'jaccard',...
        'spearman'};
        defaultVerbose= 'graphic';
        expectedVerbose = {'graphic','text','none'};
        defaultMethod='Java';
        expectedMethod={'Java', 'C', 'C vectorized', 'MATLAB', 'MATLAB vectorized',...
            'MATLAB experimental', 'MATLAB experimental 2', 'C++'};
        addOptional(p,'csv_file_or_data',[],@(x) ischar(x) || isnumeric(x));
        addParameter(p,'save_template_file',[], @ischar);
        addParameter(p,'ask_to_save_template', false, @islogical);
        addParameter(p,'randomize', false, @islogical);
        addParameter(p,'template_file',[], @(x) ischar(x) || isa(x, 'UMAP'));
        addParameter(p,'n_neighbors', 30, @(x) isnumeric(x) && x>2 && x<200);
        addParameter(p,'min_dist', .3, @(x) isnumeric(x) && x>.05 && x<.8);
        addParameter(p,'metric', defaultMetric, ...
            @(x) any(validatestring(x,expectedMetric)));
        addParameter(p,'n_epochs',[], @(x) isnumeric(x) && x>4);
        addParameter(p,'verbose',defaultVerbose,...
            @(x) any(validatestring(x,expectedVerbose)));
        addParameter(p,'method',defaultMethod,...
            @(x) any(validatestring(x,expectedMethod)));
        addParameter(p, 'parameter_names', {}, @validateParameterNames);
        addParameter(p, 'progress_callback', ...
            @(javaObject)progress_report(javaObject), @validateCallback);
        addParameter(p,'label_column',0,@(x) isnumeric(x) && x>0);
        addParameter(p,'label_file',[], @ischar);
        addParameter(p,'n_components', 2, @(x) isnumeric(x) && x>=2 && x<101);
        addParameter(p,'epsilon', .6, @(x) isnumeric(x) && x>1 && x<100);
        addParameter(p,'frequencyDensity3D', true, @islogical);
        addParameter(p,'match_supervisors', 1, @(x) isnumeric(x) && x>=0 && x<=3);
        addParameter(p,'qf_dissimilarity', false, @islogical);
        addParameter(p,'qf_tree', false, @islogical);
        addParameter(p,'joined_transform', false, @islogical);
        addParameter(p,'python', false, @islogical);
    end
end

function dispNoDbScan
        warning(['dbscan for clustering in 3+D is not available ... '...
            '\nDownload from MathWorks File Exchange: '...
            'https://www.mathworks.com/matlabcentral/fileexchange/52905-dbscan-clustering-algorithm']);
end
