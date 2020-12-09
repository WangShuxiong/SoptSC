%   AUTHORSHIP
%   Primary Developer: Stephen Meehan <swmeehan@stanford.edu> 
%   Math Lead & Secondary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Provided by the Herzenberg Lab at Stanford University 
%   License: BSD 3 clause
%
classdef Supervisors < handle
    properties(Constant)
        DEBUG_QF=false;
        DEV_UNIT_LIMIT=4;
        BORDER_DEV_UNIT_LIMIT=.2;
        VERBOSE=false;
        COLOR_NEW_SUBSET=[220 220 223];
        CLUSTER_LABELS_IF_EQUAL=false;
    end
    
    properties
        density;
        embedding;
        inputData;
        veryHighClusterDetail=true;
    end
    
    properties(SetAccess=private)
        uuid;
        sourceUuid; % universal id of source of supervising labels.  This is a
                    % immutable universally unique identifier
        sourceSubId;
        sourceDataId;
        sourceDescription;
        ids;
        labelMap;
        mdns;
        mads;
        means;
        stds;
        cnts;
        N;
        xLimit;
        yLimit;
        labels;
        nnLabels;% nearest neighbor match
        matchType=1;
    end
    
    methods
        function N=size(this)
            N=length(this.ids);
        end
        function this=Supervisors(labels, labelMap, embedding, ax)
            this.labelMap=labelMap;
            this.ids=unique(labels);
            this.labels=labels;
            N=length(this.ids);
            this.cnts=zeros(N,1);
            this.mdns=zeros(N,2);
            this.means=zeros(N,2);
            this.mads=zeros(N,2);
            this.stds=zeros(N,2);
            for i=1:N
                id=this.ids(i);
                l=labels==id;
                this.cnts(i)=sum(l);
                this.mdns(i,:)=median(embedding(l,:));
                this.mads(i,:)=mad(embedding(l,:),1);
                this.means(i,:)=mean(embedding(l,:));
                this.stds(i,:)=std(embedding(l,:),1);
            end
            this.N=N;
            if nargin<4 || isempty(ax) || ~ishandle(ax)
                mx=max(embedding);
                mn=min(embedding);
                this.xLimit=[mn(1) mx(1)];
                this.yLimit=[mn(2) mx(2)];
            else
                this.xLimit=xlim(ax);
                this.yLimit=ylim(ax);
            end
            this.embedding=embedding;
        end
        
        function setUniversalIdentifier(this, uuid, subId, dataId, description)
            assert(isempty(this.sourceUuid), ...
                'universal id is not empty (immutable means set once)');
            this.sourceUuid=uuid;
            if nargin>1 
                this.sourceSubId=subId;
                if nargin>2
                    this.sourceDataId=dataId;
                    if nargin>3
                        this.sourceDescription=description;
                    end
                end
            end
            this.uuid=java.util.UUID.randomUUID;
        end
        

        function [name, color, label]=getNameByMedian(this, mdn)
            [~, iMeHearty]=pdist2(this.density.clusterMdns, mdn, ...
                'euclidean', 'smallest', 1);
            name=this.density.clusterNames{iMeHearty};
            color=this.density.clusterColors{iMeHearty};
            label=this.density.clusterLabels{iMeHearty};
        end        
    end
    
    methods(Static)
        
        function clr=NewColor(id)
            clr=Supervisors.COLOR_NEW_SUBSET+id;            
            clr(clr<0)=0;
            clr(clr>252)=252;
        end
        
        function [mins, maxs]=GetMinsMaxs(data)
            [mins, maxs]=MatBasics.GetMinsMaxs(data, .15);
        end
        
        function [density, numClusts, clustIds]=Cluster(data, veryHigh)
            [mins, maxs]=Supervisors.GetMinsMaxs(data);
            if nargin<2 || veryHigh
                [numClusts, clustIds, density]=Density.ClusterVeryHigh(...
                    data, mins, maxs);
            else
                [numClusts, clustIds, density]=Density.ClusterMedium(...
                    data, mins, maxs);
            end
        end
    end
    
    methods
        function [names, lbls, clrs]=getQfTraining(this)
            names={};
            clrs=[];
            lbls=this.labels;
            N_=length(this.ids);
            for i=1:N_
                id=this.ids(i); 
                if id>0      
                    if this.cnts(i)>=11
                        key=num2str(id);
                        names{end+1}=this.labelMap.get(key);
                        clrs(end+1,:)=str2num(...
                            this.labelMap.get([key '.color']))/256;
                    else
                        lbls(this.labels==id)=0;
                    end
                end
            end
        end
        
        function [names, lbls, clrs]=getQfTrained(this, data)
            names={};
            clrs=[];
            [lbls, lblMap]=this.supervise(data);
            lbls(lbls<0)=0;
            ids_=unique(lbls);
            cnts_=histc(lbls, ids_);
            N_=length(ids_);
            for i=1:N_
                id=ids_(i);
                if id>0
                    if cnts_(i)>=20
                        key=num2str(id);
                        names{end+1}=lblMap.get(key);
                        clrs(end+1,:)=str2num(...
                            lblMap.get([key '.color']))/256;
                    else
                        lbls(lbls==id)=0;
                    end
                else
                    lbls(lbls==id)=0;
                end
            end
        end
        
        function qfMatchWithClusters(this, data, density, ...
                numClusters, clusterIds, pu)
            if nargin<5
                pu=[];
                if nargin<3
                    [density, numClusters, clusterIds]=Supervisors.Cluster(...
                        data, this.veryHighClusterDetail);
                end
            end
            this.density=density;
            [R,C]=size(clusterIds);
            if C>1 && R==1
                clusterIds=clusterIds';
            end
            [tNames, lbls, clrs]=this.getQfTraining;
            if isequal(data, this.embedding)
                matchStrategy=2;
            else
                matchStrategy=1;
            end
            qf=run_HiD_match(this.embedding, ...
                lbls, data, clusterIds, 'trainingNames', tNames, ...
                'matchStrategy', matchStrategy, 'log10', false, 'pu', pu);
            [tCnt, sCnt, supr1stIdx4Clue, clusterMatch]=qf.getMatches;
            if Supervisors.DEBUG_QF
                [~, ~, tMi, sMi, supr1stIdx4Clue2, supr1stId4Clue, supervisorsUnmatched, ...
                    clustersUnmatched, tNoMatchIds, sNoMatchIds]=qf.getMatches2;
                qf.tNames(supervisorsUnmatched)
                assert(isequal(tNames(ismember(qf.tIds, tNoMatchIds)),...
                    qf.tNames(supervisorsUnmatched)));
                tNames(supr1stIdx4Clue(supr1stIdx4Clue>0))
                %assert(isequal(this.ids(supr1stIdx4Clue+1), supr1stId4Clue'))
                assert(isequal(supr1stIdx4Clue,supr1stIdx4Clue2))
                assert(isequal(clustersUnmatched, sCnt==0))
                assert(isequal(supervisorsUnmatched, tCnt==0))
                [tQ, sQ, tF, sF]=qf.getScores;
                [d2, unmatched]=qf.getTableData(clrs)
                tN=length(tQ);
                sN=length(sQ);
                NN=min([tN sN]);
                for i=1:NN
                    assert(isequal(num2str(d2{i,4}), num2str(tQ(i))))
                    assert(isequal(num2str(d2{tN+i,4}), num2str(sQ(i))))
                    assert(isequal(num2str(d2{i,5}), num2str(tF(i))))
                    assert(isequal(num2str(d2{tN+i,5}), num2str(sF(i))))
                    assert(d2{i, 8}==qf.tSizes(i))
                    assert(d2{tN+i, 8}==qf.sSizes(i))
                end
            end
            cluMdns=zeros(numClusters, 2);
            clusterLabels=cell(1,numClusters);
            clusterNames=cell(1,numClusters);
            clusterColors=cell(1,numClusters);
            newSubsets=0;
            labels_=zeros(size(data, 1), 1);
            for i=1:numClusters
                l=clusterIds==i;
                if sCnt(i)==0
                    newSubsets=newSubsets+1;
                    clusterLabel=0-i;
                    clusterNames{i}=['New subset #' num2str(newSubsets) ];
                    clr=num2str(Supervisors.NewColor(newSubsets));
                else
                    clusterLabel=clusterMatch(i);
                    clusterNames{i}=tNames{supr1stIdx4Clue(i)};
                    clr=this.labelMap.get([num2str(clusterLabel) '.color']);
                    if Supervisors.DEBUG_QF
                        assert(isequal(clusterNames{i}, ...
                            this.labelMap.get(num2str(clusterLabel))))
                    end
                end
                clusterColors{i}=clr;
                clusterLabels{i}=clusterLabel;
                labels_(l)=clusterLabel;
                if Supervisors.VERBOSE
                    sum(l)
                end
                cluMdns(i,:)=median(data(l,:));
            end
            this.density.setLabels(labels_, clusterNames, ...
                clusterColors, cluMdns, clusterLabels);
            

            function d=normalize(d)
                mn=min(d);
                N_=length(mn);
                for j=1:N_
                    if mn(j)<=0
                        add_=1-mn(j);
                        d(:,j)=d(:,j)+add_;
                    end
                end
                for j=1:N_
                    mx_=max(d(:,j));
                    mn_=min(d(:,j));
                    r=mx_-mn_;
                    d(:,j)=(d(:,j)-mn_)/r;
                end
            end
        end
        
        function nnMatchWithClusters(this, data, density, ...
                numClusters, clusterIds)
            this.density=density;
            this.computeNearestNeighbors(data);
            clusterLabels=cell(1,numClusters);
            clusterNames=cell(1,numClusters);
            clusterColors=cell(1,numClusters);
            labels_=zeros(size(data, 1), 1);
            cluMdns=zeros(numClusters, 2);
            newSubsets=0;
            for i=1:numClusters
                l=clusterIds==i;
                cluMdns(i,:)=median(data(l,:));
                clueLabels=this.nnLabels(l);
                u=unique(clueLabels);
                clusterLabelCnts=histc(clueLabels, u);
                [mxCnt, mxI]=max(clusterLabelCnts);
                lbl=u(mxI);
                if lbl==0
                    newSubsets=newSubsets+1;
                    lbl=0-i;
                    clusterNames{i}=['New subset #' num2str(newSubsets) ];
                    clusterColors{i}=num2str(...
                        Supervisors.NewColor(newSubsets));
                else
                    key=num2str(lbl);
                    clusterNames{i}=this.labelMap.get(key);
                    clusterColors{i}=this.labelMap.get([key '.color']);
                end
                labels_(l)=lbl;
                clusterLabels{i}=lbl;
                fprintf(['%s (id=%d) has %d/%d events in cluster '...
                    '#%d''s %d events\n'], clusterNames{i}, ...
                    lbl, mxCnt, ...
                    sum(this.nnLabels==lbl), i, sum(l));
            end
            this.density.setLabels(labels_, clusterNames, ...
                clusterColors, cluMdns, clusterLabels);
        end
                
        function matchWithClusters(this, data, density, ...
                numClusters, clusterIds)
            this.density=density;
            clusterLabels=cell(1,numClusters);
            clusterNames=cell(1,numClusters);
            clusterColors=cell(1,numClusters);
            %[~,I]=pdist2(this.mdns, embedding, 'euclidean', 'smallest', 1);
            %labels2=this.ids(I);
            avgs=this.mdns;
            if isprop(this, 'mads')
                devs=this.mads;
            else
                devs=[];
            end
            if any(this.cnts<20)
                avgs(this.cnts<20, 1)=this.xLimit(2)*5;
                avgs(this.cnts<20, 2)=this.yLimit(2)*5;
            end
            cluMdns=zeros(numClusters, 2);
            cluDevs=zeros(numClusters, 2);
            for i=1:numClusters
                l=clusterIds==i;
                if Supervisors.VERBOSE
                    sum(l)
                end
                cluMdns(i,:)=median(data(l,:));
                cluDevs(i,:)=mad(data(l,:), 1);
            end
            hasDevs=~isempty(devs);
            [D, I]=pdist2(avgs, cluMdns, 'euclidean', 'smallest', 1);
            labels_=zeros(size(data, 1), 1);
            reChecks={};
            for i=1:numClusters
                labelIdx=I(i);
                label=this.ids(labelIdx);
                if label==0
                    label=0-i;
                else
                    key=num2str(label);
                    clusterLabels{i}=label;
                    clusterNames{i}=this.labelMap.get(key);
                    clusterColors{i}=this.labelMap.get([key '.color']);
                end
                l=clusterIds==i;
                if Supervisors.VERBOSE
                    sum(l)
                    this.labelMap.get(num2str(label))
                end
                if hasDevs
                    devDist=MatBasics.DevDist(cluMdns(i,:), cluDevs(i,:));
                    if any(D(i)>=devDist*Supervisors.DEV_UNIT_LIMIT)
                        reChecks{end+1}=struct('clustId', i, 'count',...
                            sum(l), 'label', label, 'labelIdx', labelIdx);
                        label=0-i;
                    end
                end
                labels_(l)=label;
            end
            if hasDevs
                N_=length(reChecks);
                while N_>0 
                    changes=[];
                    for i=1:N_
                        clustId=reChecks{i}.clustId;
                        label=reChecks{i}.label;
                        labelIdx=reChecks{i}.labelIdx;
                        closestLabelIdxs=labels_==label;
                        if Supervisors.VERBOSE
                            sum(closestLabelIdxs)
                            this.labelMap.get(num2str(label))
                            reChecks{i}
                        end
                        if any(closestLabelIdxs)
                            %Does this cluster with no label match
                            %   sit on the border of one with the 
                            %   closest label from the supervisor?
                            unlabeledClusterIdxs=clusterIds==clustId;
                            borderDistance=min(pdist2(...
                                data(unlabeledClusterIdxs, :), ...
                                data(closestLabelIdxs,:), 'euclidean', ...
                                'smallest', 1));
                            supervisorDevDistance=MatBasics.DevDist(...
                                avgs(labelIdx,:), devs(labelIdx,:));
                            limit=supervisorDevDistance*...
                                Supervisors.BORDER_DEV_UNIT_LIMIT;
                            if borderDistance<=limit
                                changes(end+1)=i;
                                labels_(unlabeledClusterIdxs)=label;
                            end
                        end
                    end
                    if isempty(changes)
                        break;
                    else
                        reChecks(changes)=[];
                    end
                    N_=length(reChecks);                    
                end
                newSubsetIds=unique(labels_(labels_<0));
                newSubsets=length(newSubsetIds);
                for i=1:newSubsets
                    clustId=0-newSubsetIds(i);
                    clusterLabels{clustId}=newSubsetIds(i);
                    clusterNames{clustId}=['New subset #' num2str(i)];
                    color_=Supervisors.NewColor(clustId);
                    if any(color_<0)
                        color_(color_<0)=0;
                    end
                    clusterColors{clustId}=num2str(color_);
                end
                N_=length(this.ids);
                remainder=data(labels_==0,:); 
                labels2=zeros(size(remainder,1),1);
                for i=1:N_
                    if ~any(find(labels_==this.ids(i),1))
                        label=this.ids(i);
                        if Supervisors.VERBOSE
                            this.labelMap.get(num2str(label))
                        end
                        [D2, I2]=pdist2(avgs(i,:), remainder, ...
                            'euclidean', 'smallest', 1);
                        pt=avgs(i,:)+devs(i,:);
                        devDist=pdist2(avgs(i,:), pt);
                        if any(D2<devDist*Supervisors.DEV_UNIT_LIMIT)
                            labels2(D2<devDist)=label;
                        end
                    end
                end
                if any(labels2)
                    labels_(labels_==0)=labels2;
                end
            end
            this.density.setLabels(labels_, clusterNames, ...
                clusterColors, cluMdns, clusterLabels);
        end
        
        function computeAndMatchClusters(this, data, matchType, pu)
            if nargin<4
                pu=[];
                if nargin<3
                    matchType=this.matchType;
                end
            end
            [dns, numClusters, clusterIds]=Supervisors.Cluster(...
                data, this.veryHighClusterDetail);
            this.matchClusters(data, dns, numClusters, clusterIds, ...
                matchType, pu);
        end
        
        function setMatchType(this, matchType)
            this.matchType=matchType;
        end
        
        function matchClusters(this, data, dns, numClusters, clusterIds, ...
                matchType, pu)
            if nargin<7
                pu=[];
                if nargin<6
                    matchType=this.matchType;
                else
                    this.matchType=matchType;
                end
            end
            if matchType==0
                this.matchWithClusters(data, dns, numClusters, clusterIds);
            elseif matchType==1
                this.qfMatchWithClusters(data, dns, numClusters, clusterIds, pu);
            elseif matchType>=2
                this.nnMatchWithClusters(data, dns, numClusters, clusterIds);
            end
        end
        
        function nnLbls=prepareForTemplate(this)
            nnLbls=this.nnLabels;
            this.nnLabels=[];
        end
        
        function setNearestNeighborLabels(this, nnLabels)
            this.nnLabels=nnLabels;
        end
        
        function resetNearestNeighbors(this, data)
            this.nnLabels=[];
            this.computeNearestNeighbors(data);
        end
        
        function computeNearestNeighbors(this, data)
            if length(this.nnLabels)~=size(data,1)
                this.nnLabels=[];
            end
            if isempty(this.nnLabels)
                if ~isequal(this.embedding, data)
                    pu=PopUp('Comupting nearest neighbors');
                    [~,II]=pdist2(this.embedding, data, 'euclidean', 'Smallest', 1);
                    this.nnLabels=this.labels(II);
                    pu.close;
                else
                    this.nnLabels=this.labels;
                end
            end
        end
        
        function [labels, labelMap, reSupervised]=supervise(...
                this, data, doHtml, matchType)
            reSupervised=false;
            if nargin<4
                matchType=this.matchType;
                if nargin<3
                    doHtml=false;
                end
            else
                this.matchType=matchType;
            end
            if length(this.density.labels) ~= size(data, 1)
                reSupervised=true;
                this.computeAndMatchClusters(data, matchType)
            end
            if matchType==3% nearest neighbor match alone (no cluster matching)
                this.computeNearestNeighbors(data);
                labels=this.nnLabels;
                labelMap=this.labelMap;
            else %0=median OR 2=nearest neighbor match to clusters
                labelMap=java.util.Properties;
                labels=this.density.labels;
                doMap(this.density.clusterColors, this.density.clusterNames);
            end
            
            function doMap(clusterColors, clusterNames)                
                ids_=unique(labels);
                N_=length(ids_);
                for i=1:N_
                    putInMap(ids_(i), clusterColors, clusterNames);
                end
            end
            
            function putInMap(id, clusterColors, clusterNames) 
                key=num2str(id);
                keyColor=[key '.color'];
                if id==0
                    if doHtml
                        name='<font color="62A162">unsupervised</font>';
                    else
                        name='\color[rgb]{0.4 0.65 0.4}\bf\itunsupervised';
                    end
                    color='92 92 128';
                elseif id<0
                    clustId=0-id;
                    nm=clusterNames{clustId};
                    if doHtml
                        name=['<font color="#4242BA"><i>' nm ' ?</i></font>'];
                    else
                        name=['\color[rgb]{0. 0.4 0.65}\bf\it' nm  ' ?'];
                    end
                    color=clusterColors(clustId);
                else
                    name=this.labelMap.get(key);
                    if doHtml
                        if String.Contains(name, '^{')
                            name=strrep(name, '^{', '<sup>');
                            name=strrep(name, '}', '</sup>');
                        end
                    end
                    color=this.labelMap.get(keyColor);
                end
                labelMap.put(java.lang.String(key), name);
                labelMap.put(keyColor, color);
            end
            
        end
        
        function [qf, qft]=qfTreeSupervisors(this, visible, pu, ttl)
            if nargin<4
                ttl='UMAP template''s training set';
                if nargin<3
                    pu=[];
                    if nargin<2
                        visible=true;
                    end
                end
            end
            if isempty(this.inputData)
                msg(Html.WrapHr(['Parameter reduction is out of date!<br>'...
                    '...redo with <u>this version</u> of the software!']));
                qf=[];
                qft=[];
                return;
            end
            [tNames, lbls, clrs]=this.getQfTraining;
            [qf, qft]=run_QfTree(this.inputData, lbls, {ttl}, ...
                'trainingNames', tNames, 'log10', true, 'colors', clrs, ...
                'pu', pu);
            if visible
                set(qft.fig, 'visible', 'on');
            end
        end
        
        function [qf, qft]=qfTreeSupervisees(this, embedding, rawData, ...
                visible, pu, ttl)
            if nargin<6
                ttl='"test set" found by UMAP supervised template';
                if nargin<5
                    pu=[];
                    if nargin<4
                        visible=true;
                    end
                end
            end
            [tNames, lbls, clrs]=this.getQfTrained(embedding);
            [qf, qft]=run_QfTree(rawData, lbls, {ttl}, ...
                'pu', pu, 'trainingNames', tNames, ...
                'log10', true, 'colors', clrs);
            if visible
                set(qft.fig, 'visible', 'on');
            end
        end

        function [qf, qft]=qfDissimilarity(this, reducedData, ...
                unreducedData, visible, pu, file)
            if nargin<6
                file=[];
                if nargin<5
                    pu=[];
                    if nargin<4
                        visible=true;
                    end
                end
            end
            if isempty(this.inputData)
                msg(Html.WrapHr(['Parameter reduction is out of date!<br>'...
                    '...redo with <u>this version</u> of the software!']));
                return;
            end
            [tNames, tLbls, clrs]=this.getQfTraining;
            matchStrategy=1;
            [sNames, sLbls]=this.getQfTrained(reducedData);
            if isempty(file) || ~exist(file, 'file')
                qf=run_HiD_match(this.inputData, ...
                    tLbls, unreducedData, sLbls, 'trainingNames', tNames, ...
                    'matchStrategy', matchStrategy, 'log10', true, ...
                    'testNames', sNames, 'pu', pu);
            else
                qf=QfTable.Load(file);
            end
            qft=QfTable(qf, clrs, [], gcf, visible);
            qft.doHistQF(visible);
            if ~isempty(file) && ~exist(file, 'file')
                qft.save(qf, file);
            end
        end
        
        
        function drawClusterBorders(this, ax)
            wasHeld=ishold(ax);
            if ~wasHeld
                hold(ax, 'on');
            end
            N_=length(this.density.clusterColors);
            for i=1:N_
                clr=(str2num(this.density.clusterColors{i})/256)*.85;
                gridEdge(this.density, true, i, clr, ax, .8, '.', '-', .5);
                if Supervisors.VERBOSE
                    str2num(this.density.clusterColors{i})
                    clr
                    disp('ok');
                end
            end
            if ~wasHeld
                hold(ax, 'off');
            end
        end

    end
end