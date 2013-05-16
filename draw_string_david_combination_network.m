%%%
% Run this file to see a demo of how to combine STRING network interactions
% with DAVID functional annotation clustering.


% STRING: go to the string website (http://string-db.org), click multiple
% names, enter your list of gene symbols, click save (disk icon below
% network), and save the "Text Summary  (TXT - simple tab delimited
% flatfile)" on your local drive. 

% DAVID: go to the david website (http://david.abcc.ncifcrf.gov/tools.jsp),
% select "upload" on the left panel, (step 1) paste your list of gene symbols into
% the field, (step 2) select "official_gene_symbol" as gene identifier,
% (step 3) clik on "gene list", and (step 4) submit list. On the left, be
% sure to select the correct species, and clik the "Select Species"
% buttion. Then click on the  "functional annotation clustering" link, and
% click on the "functional annotation clustering" button. Set the
% classification stringency to "high" and press "rerun using options". Now
% click the "download file" link

% This exmaply only works if gene symbols are fed into DAVID, as STRING
% returns gene symbols.

% Berend Snijder, 2011.

%%%%%%%%%%%%%%%%%%%%%%%%%%
% a couple of parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%

% minimum global string interaction score
intMinString = 0.9;
% minimum within-david-cluster string interaction score
intWithinDavidMinString = 0;

% minimum david cluster overlap score
intMinDavidOverlap = 0.9;
% minimum david cluster enrichment score
intMinDavidEnrichment = 0.3;

% wether to run demo mode, or ask the user to locate the files...
boolRunDemo = true;
% boolRunDemo = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%


% Get the string network and the corresponding gene id and symbol seed
% lists.
if boolRunDemo
    strSTRINGFile = which('new_string9_genesymbols.txt');
else
    strSTRINGFile = '';
end
[matStringNetwork,matStringIdentifiers,cellStringGeneSymbols,cellENSPs] = parseStringNetwork(strSTRINGFile);

% Load the david functional annotation clustering data. 
if boolRunDemo
    strDAVIDFile = which('new_david_gene_symbols_high.txt');
else
    strDAVIDFile = '';
end
[cellLabelPerCluster,cellGenesPerCluster,cellGenesPerClassPerCluster,matEnrichmentScorePerCluster,matOverlapMap]= parseDavidFAC(strDAVIDFile);

% we need to combine the genes listed in STRING file and those in the
% DAVID file to get to a minimal list of all present gene identifiers.
% Note that genes regularly disappear due to gene identifier mapping
% mismatches. 
cellSeedGeneSymbols = unique(cat(1,cellGenesPerCluster{:},cellStringGeneSymbols));

% reformat string network according to full seed gene list
[~,b]=ismember(cellStringGeneSymbols,cellSeedGeneSymbols);
matStringNetwork2 = zeros(numel(cellSeedGeneSymbols));
matStringNetwork2(b,b) = matStringNetwork;
matStringNetwork = matStringNetwork2;


% identify kinase (or other 'special') cluster ids, that you want to
% disband.
matKinaseClusterIDs = find(cellfun(@(x) ~isempty(regexpi(x,'Nucleotide.*binding')),cellLabelPerCluster) | ...
                        cellfun(@(x) ~isempty(regexpi(x,'phosphate metabolic process')),cellLabelPerCluster) | ...
                        cellfun(@(x) ~isempty(regexpi(x,'phosphorylation')),cellLabelPerCluster));


% initialize combined DAVID and STRING graph
intDAVIDNodes = size(matOverlapMap,1);
intSTRINGNodes = size(matStringNetwork,1);

matCombinedGraph = zeros(intDAVIDNodes+intSTRINGNodes);

% add in string network above total threshold
matGRAPH1 = triu(matStringNetwork);
matGRAPH1(matGRAPH1<intMinString) = 0;
matCombinedGraph(1:intSTRINGNodes,1:intSTRINGNodes) = matGRAPH1;

% now we need to add the links between the genes and the clusters, but of
% course only to the strongest enriched networks...
matStringDavidLinks = zeros(intSTRINGNodes,intDAVIDNodes);

% add gene to david links, but skip those for kinase-like classes
cellAdjustedGenesPerCluster = cellGenesPerCluster;
cellAdjustedGenesPerCluster(matKinaseClusterIDs) = {[]};
cellGenesPerClassPerCluster(matKinaseClusterIDs,:) = {[]};

matGeneClusterID = nan(numel(cellSeedGeneSymbols),1);
cellGeneShape = cell(numel(cellSeedGeneSymbols),1);
matEmptyClassPerClusterIX = cellfun(@isempty,cellGenesPerClassPerCluster);
matClassPerClusterCount = sum(~matEmptyClassPerClusterIX,2);
for i = 1:numel(cellSeedGeneSymbols)
    
    % this is the alternative function, look for specificity index of gene
    % in cluster (i.e. how many 
    matPerClassPerClusterGeneCount = NaN(size(cellGenesPerClassPerCluster));
    matPerClassPerClusterGeneCount(~matEmptyClassPerClusterIX) = single(cellfun(@(x) ismember(cellSeedGeneSymbols{i},x),cellGenesPerClassPerCluster(~matEmptyClassPerClusterIX)));

    % different metrics according to which we could score... the first
    % cluster with the highest mean presence of all classes is preferred.
    matMeans = sum(matPerClassPerClusterGeneCount==1,2)./matClassPerClusterCount;
    matMeans(isnan(matMeans))=0;
    
    %  "max" finds the indices of the maximum values of A, and returns them
    %  in output vector I. If there are several identical maximum values,
    %  the index of the first one (i.e. highest enriched) found is
    %  returned.
    [intMaxValue,intCrossIX]  = max(matMeans);
    
    % if we don't find a maximum, set to [], i.e. flag it's unclustered
    if intMaxValue == 0; intCrossIX = []; end
    intCrossIX2 = intCrossIX;
    matStringDavidLinks(i,intCrossIX) = 2;
    
    % annotate david cluster to gene node mapping
    if isempty(intCrossIX)
        matGeneClusterID(i) = -1;
    else
        matGeneClusterID(i) = intCrossIX;
    end

    % create gene node shapes
    if ismember(intCrossIX2,matKinaseClusterIDs)
        cellGeneShape{i} = 'kinase';
    elseif isempty(intCrossIX)
        cellGeneShape{i} = 'homeless';
    else
        cellGeneShape{i} = 'regular';
    end    
    
end


%%%
% and add string interactions that are above a lower threshold for all
% genes present in the same david cluster.
for i = unique(matGeneClusterID')
    a = find(matGeneClusterID==i);
    matGRAPH1 = triu(matStringNetwork(a,a));
    matGRAPH1(matGRAPH1<intWithinDavidMinString) = 0;
    matCombinedGraph(a,a) = matGRAPH1;
end
%%%

%%%
% add links between genes and david classes
matCombinedGraph(1:intSTRINGNodes,intSTRINGNodes+1:end) = matStringDavidLinks;
% % keep a edge-feature to annotate gene-david links
matGeneToDavidEdges = zeros(size(matCombinedGraph));
matGeneToDavidEdges(1:intSTRINGNodes,intSTRINGNodes+1:end) = matStringDavidLinks;

% merge all labels for final combined graph
cellCombinedLabels = [cellSeedGeneSymbols;cellLabelPerCluster];
matCombinedClusterIDs = [matGeneClusterID;(1:intDAVIDNodes)'];
cellCombinedNodeShapes = [cellGeneShape;repmat({'cluster'},[intDAVIDNodes,1])];

% generate an attribute per gene that lists the class per gene
cellClusterPerGene = cell(size(cellCombinedLabels));
cellClusterPerGene(matCombinedClusterIDs>0) = cellLabelPerCluster(matCombinedClusterIDs(matCombinedClusterIDs>0));
cellClusterPerGene(matCombinedClusterIDs<=0) = {'Unclustered'};

% call cytoscape with predefined example vizmap file
strVizmapFile = which('vizmap_combined.props');
cytoscape(matCombinedGraph,cellCombinedLabels,matCombinedClusterIDs,cellCombinedNodeShapes,matGeneToDavidEdges,cellClusterPerGene,strVizmapFile)


