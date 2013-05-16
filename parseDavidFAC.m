function [cellLabelPerCluster,cellGenesPerCluster,cellGenesPerClassPerCluster,matEnrichmentScorePerCluster,matOverlapMap] = parseDavidFAC(strFilePath,strMappingOption)
% help parseDavidFAC
%
% usage:
%
% [cellLabelPerCluster,cellGenesPerCluster,cellGenesPerClassPerCluster,matEnrichmentScorePerCluster,matOverlapMap] = parseDavidFAC(strFilePath,strMappingOption)
%
% where the strFilePath input is a string describing the path to an excel
% or tab-separated text file output of DAVID Functional Annotation
% Clustering. 
%
% In this file, optionally, the third column per header of each cluster can
% contain a manual cluster annotations. If this field is empty (default),
% the first annotation class will be the cluster name.
%
% strMappingOption is optional and can describe how to assign genes to
% clusters. Default = 'intersect', alternatively 'union' can be passed.
% 'Union' assigns all genes observed at least once to a cluster, whereas
% 'intersect' assigns only those genes present in all functional annotation
% classes per cluster. 
%
%
% Output:
%
% cellLabelPerCluster = n x 1 cellarray with labels per cluster.
%
% cellGenesPerCluster = n x 1 cellarray with genes per cluster.
%
% cellGenesPerClassPerCluster = n x m cellarray with genes per class (m is
%    highest class number per cluster observed) of each cluster (n). 
%
% matEnrichmentScorePerCluster = n x 1 matrix with DAVID EASE score per
%    cluster 
%
% matOverlapMap = n x n matrix, where n is the number of clusters, with for
%    each cluster pair the fraction of genes overlapping between the two
%    clusters. The fraction is calculated as the fraction of genes from the
%    smaller cluster
%
% [copyright Berend Snijder, 2011]

% init output
cellLabelPerCluster = {};
cellGenesPerCluster = {};
cellGenesPerClassPerCluster = {};
matEnrichmentScorePerCluster = [];
matOverlapMap = [];

% if no input or incorrect input is given, ask user for DAVID file.
if nargin==0 || isempty(strFilePath) || ~fileattrib(strFilePath)
    [FileName,PathName] = uigetfile( {'*.txt;*.xls;*.xlsx','DAVID Functional Annotation Clustering file (*.txt,*.xls,*.xlsx)';
   '*.txt',  'Text files (*.txt)'; ...
   '*.xls','Excel files (*.xls)'; ...
   '*.xlsx','Excel files (*.xlsx)'; ...
   '*.*',  'All Files (*.*)'}, ...
   'Please select a DAVID Functional Annotation Clustering file');
    if isequal(FileName,0) || ~fileattrib(fullfile(PathName,FileName))
        warning('bs:Bla','Next time please select a valid DAVID Functional Annotation Clustering file \n ')
        return
    else
        strFilePath = fullfile(PathName,FileName);
    end
end

% optional input on how to parse genes per cluster, default to 'union'
if nargin<2
    strMappingOption = 'intersect';
end

% report inputs
switch lower(strMappingOption)
    case 'union'
        fprintf('%s: parsing ''%s'' with ''%s'' (default).\n',mfilename,strFilePath,strMappingOption)
    case 'intersect'
        fprintf('%s: parsing ''%s'' with ''%s''.\n',mfilename,strFilePath,strMappingOption)
    otherwise
        warning('bs:Bla','unkown mapping option ''%s'', using default ''union''.\n',strMappingOption)
        strMappingOption = 'union';
        fprintf('%s: parsing ''%s'' with ''%s''.\n',mfilename,strFilePath,strMappingOption)
end


% read the DAVID functional annotation file (either excel or text)
if ~isempty(strfind(lower(strFilePath),'.xls'))
    fprintf('%s: \treading file as excel.\n',mfilename)
    [~,~,foo] = xlsread(strFilePath);
else
    fprintf('%s: \treading file as tab-delimited.\n',mfilename)
    % try to read as text tab-separated file
    fid = fopen(strFilePath);
    
    % the following code should give equal output as xlsread
    foo = {};
    iLine = 0;
    while true

        % read in the next line
        tline = fgetl(fid);

        % stop if this line isnt text, at eof it's -1
        if ~ischar(tline); break; end

        % skip empty lines
        if isempty(tline)
            iLine = iLine + 1;
            continue
        end
        
        C = textscan(tline, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s','delimiter','\t','BufSize',10000);
        matEmptyIX = cellfun(@isempty,C);
        C(~matEmptyIX) = cellfun(@(x) x{1}, C(~matEmptyIX), 'UniformOutput',false);
        C(matEmptyIX) = {[]};
        % convert strings describing numbers to numbers
        matPotentialNumeric = cellfun(@(x) isnumeric(str2double(x)) && ~isnan(str2double(x)), C);
        C(matPotentialNumeric) = cellfun(@(x) str2double(x), C(matPotentialNumeric),'UniformOutput',false);
        
        iLine = iLine + 1;
        foo(iLine,:) = C;
    end
    fclose(fid);    
end




cellLabelColumn = foo(:,2);
cellManualLabelColumn = foo(:,3);
cellGenesColumn = foo(:,6);

cellLabelPerCluster = cell(400,1);
cellGenesPerCluster = cell(400,1);
cellGenesPerClassPerCluster = cell(400,10);

iCluster = 0;

% get enrichment score per DAVID cluster
matIX = cellfun(@(x) ~isempty(strfind(x,'Annotation Cluster')),foo(:,1));
matEnrichmentScorePerCluster = cellfun(@(x) cell2mat(textscan(x,'Enrichment Score: %f')),foo(matIX,2));
 
for i = 1:length(cellGenesColumn)

    if strcmpi(cellGenesColumn{i},'Genes')
        iCluster = iCluster + 1;
        cellGenesPerCluster{iCluster,1} = [];
        
        if ~isempty(cellManualLabelColumn) & (numel(cellManualLabelColumn) >= i)
            % note that the manual label is placed one row above the
            % 'Genes' label containing row...
            if ~isempty(cellManualLabelColumn{i-1}) & ~isnumeric(cellManualLabelColumn{i-1}) & ~strcmpi(cellManualLabelColumn{i-1},'Count')
                cellLabelPerCluster{iCluster,1} = cellManualLabelColumn{i-1};
            else
                cellLabelPerCluster{iCluster,1} = {};
            end
        end
    else
        if ~isnan(cellGenesColumn{i}) & ~isempty(cellGenesColumn{i}) & ~strcmpi(cellGenesColumn{i},'Genes')

            % if we didn't find a 'gene' header, we didn't find any clusters...
            % throw error.
            if iCluster==0
                errordlg(sprintf('Failed to parse file ''%s''\n\n Are you sure it''s a DAVID Functional Annotation Cluster file?',strFilePath))
                error('%s: \t Failed to parse DAVID file ''%s'' (No clusters found)\n',mfilename,strFilePath)
            end            
            
            % DAVID nowadays separates identifiers by comma. Try numeric
            % first.
            matCurrentGenes = cell2mat(textscan(cellGenesColumn{i},'%d','delimiter',','));
            
            % If the file does not contain numeric gene identifiers, switch
            % to this.
            if numel(matCurrentGenes)<2
                matCurrentGenes = textscan(cellGenesColumn{i},'%s','delimiter',',');
                matCurrentGenes = matCurrentGenes{1};
            end
            
            switch lower(strMappingOption)
                case 'union'
                    cellGenesPerCluster{iCluster,1} = union(cellGenesPerCluster{iCluster,1},matCurrentGenes);
                case 'intersect'
                    if isempty(cellGenesPerCluster{iCluster,1})
                        cellGenesPerCluster{iCluster,1} = matCurrentGenes;
                    else
                        cellGenesPerCluster{iCluster,1} = intersect(cellGenesPerCluster{iCluster,1},matCurrentGenes);
                    end
                otherwise
                    warning('bs:Bla','Unknown mapping option given...')
                    cellGenesPerCluster{iCluster,1} = union(cellGenesPerCluster{iCluster,1},matCurrentGenes);
            end
                    
            
            % also store the genes for all the classes in a cluster
            matClassIX = find(cellfun(@isempty,cellGenesPerClassPerCluster(iCluster,:)),1,'first');
            if isempty(matClassIX); matClassIX = size(cellGenesPerClassPerCluster,2) + 1;end
            cellGenesPerClassPerCluster{iCluster,matClassIX} = matCurrentGenes;

            % we can annotate each cluster with the first - most enriched -
            % annotation (alternatively the shortest also works nice :)
            if isempty(cellLabelPerCluster{iCluster,1})% | length(cellLabelColumn{i}) < cellLabelPerCluster{iCluster,1}
                strLabel = cellLabelColumn{i};
                strLabel(1:find(ismember(strLabel,':~'),1,'last')) = [];
                cellLabelPerCluster{iCluster,1} = sprintf('[%d] %s (E=%.1f)',iCluster,strLabel,matEnrichmentScorePerCluster(iCluster));
            end
        end
        
    end
    
end

% some reformatting and removing of trailing / empty clusters
matEmptyRowsIX = cellfun(@isempty,cellGenesPerCluster) & cellfun(@isempty,cellLabelPerCluster);
cellGenesPerCluster(matEmptyRowsIX)=[];
% make sure each gene list is unique...
cellGenesPerCluster = cellfun(@unique,cellGenesPerCluster,'UniformOutput',false);
cellGenesPerClassPerCluster(matEmptyRowsIX,:)=[];
cellLabelPerCluster(matEmptyRowsIX)=[];

% add the number of genes in each DAVID cluster to the cluster labels
cellLabelPerCluster = cellfun(@(x,y) sprintf('%s (n=%d)',x,numel(y)),cellLabelPerCluster,cellGenesPerCluster,'UniformOutput',false);


%%%% calculate overlap map
% Add between class links if their genes overlap a lot (bi-clustering in a
% network). This also draws/adds otherwise ignored DAVID clusters.
matOverlapMap = zeros(length(cellGenesPerCluster));
for i = 1:length(cellGenesPerCluster)
    for ii = i:length(cellGenesPerCluster)
        
        % skip self links and process each link only once
        if i==ii;continue;end

        % calculate fraction of overlap as function of smallest cluster
        % size
        matFractionOfSmallestIdentical = sum(ismember(cellGenesPerCluster{i},cellGenesPerCluster{ii})) / ...
            min(length(cellGenesPerCluster{i}),length(cellGenesPerCluster{ii}));

        % make graph for cytoscape visualization
        matOverlapMap(i,ii) = matFractionOfSmallestIdentical;
        matOverlapMap(ii,i) = matFractionOfSmallestIdentical;

    end
end

