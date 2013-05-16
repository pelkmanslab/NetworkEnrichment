
% parse the output from the DAVID Functional Annotation Clustering and draw
% it in a Cytoscape network.
%
% see for DAVID: http://david.abcc.ncifcrf.gov/summary.jsp
% see for CYTOSCAPE: http://www.cytoscape.org

% This is an example david file (passing gene symbols). If the file is not
% found, parseDavidFAC will ask the user to locate a valid file.

% Berend Snijder, 2011.

%%%%%%%%%%%%%%%%%%%%%%%%%%
% a couple of parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%

% draw all overlap edges above the first threshold
intPrimaryThreshold = 0.75;
% and add the single highest edge for all unconnected nodes
intSecondaryThreshold = 0.1;

% wether to run demo mode, or ask the user to locate the files...
boolRunDemo = true;
% boolRunDemo = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%


if boolRunDemo
    strDAVIDFile = which('new_david_gene_symbols_high.txt');
else
    strDAVIDFile = '';
end

[cellLabelPerCluster,cellGenesPerCluster,cellGenesPerClassPerCluster,matEnrichmentScorePerCluster,matOverlapMap]= parseDavidFAC(strDAVIDFile);


% backup the original matOverlapMap
matOverlapMapBackup = matOverlapMap;

% kick out weakest links
matOverlapMap(matOverlapMap<intPrimaryThreshold) = 0;

% but add back for the unconnected links their strongest link
matUnconnectedIX = max(matOverlapMap')==0 & max(matOverlapMap)==0; %#ok<UDIM>
[matUnconnectedMaxVal,matMaxIX] = max(matOverlapMapBackup'); %#ok<UDIM>

% add back the single-strongest link (above the secondary threshold) for
% otherwise unconnected functional annotation clusters.
matBelowThresholdLinks = zeros(size(matOverlapMap));
for i = find(matUnconnectedIX & matUnconnectedMaxVal > intSecondaryThreshold)
    matOverlapMap(matMaxIX(i),i) = matUnconnectedMaxVal(i);
    matBelowThresholdLinks(matMaxIX(i),i) = true;
end

% make sure we only have single links in the network.
matOverlapMap = triu(matOverlapMap);

%%%%%%%%%
% open in cytoscape
matClusterSize = cellfun(@numel,cellGenesPerCluster);
cellClusterNumber = arrayfun(@(x) num2str(x),1:numel(matClusterSize),'UniformOutput',false)';
strVizmapFile = which('vizmap.props');
cytoscape(matOverlapMap,cellLabelPerCluster,matEnrichmentScorePerCluster,matClusterSize,cellClusterNumber,strVizmapFile,matBelowThresholdLinks)
%%%%%%%%%


