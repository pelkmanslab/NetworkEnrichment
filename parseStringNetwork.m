function [matStringNetwork,matStringIdentifiers,cellStringGeneSymbols,cellENSPs] = parseStringNetwork(strStringNetworkFile)
% help parseStringNetwork
%
% usage:
%
% [matStringNetwork,matStringIdentifiers,cellStringGeneSymbols,cellENSPs] = parseStringNetwork(strStringNetworkFile)
%
% where the strStringNetworkFile input is a string describing the path to
% a tab-separated text file stored from the STRING website string-db.org.
%
% Output:
%
% matStringNetwork = n x n matrix with network where edge-scores are STRING
% scores. 0 means no edge. n is the number of genes in your network.
%
% matStringIdentifiers = n x 1 matrix with STRING identifiers
%
% cellStringGeneSymbols = n x 1 cellarray with gene symbols per gene
%
% cellENSPs = n x 1 cellarray with ensemble protein ids per gene
%
% [copyright Berend Snijder, 2011]

% init output
matStringNetwork = [];
matStringIdentifiers = [];
cellStringGeneSymbols = {};
cellENSPs = {};

if nargin==0 || isempty(strStringNetworkFile) || ~fileattrib(strStringNetworkFile)
    [FileName,PathName] = uigetfile( {'*.txt','STRING Network file (*.txt)'}, 'Please select a STRING network file');
    if isequal(FileName,0) || ~fileattrib(fullfile(PathName,FileName))
        warning('bs:Bla','Next time please select a valid STRING network file \n ')
        return
    else
        strStringNetworkFile = fullfile(PathName,FileName);
    end
end

fprintf('%s: parsing file ''%s''\n',mfilename,strStringNetworkFile)

% open the string network file for reading
fid = fopen(strStringNetworkFile,'r');


% read in the header of file
tline = fgetl(fid);

i = 0;
cellData = cell(200,15);
while true 
    
    % read in next line (and skip header)
    tline = fgetl(fid);
    
    if ~ischar(tline); break; end
    
    i = i + 1;
    % this is the correct formatting
    try
        cellData(i,:) = textscan(tline, '%s %s %d %d %s %s %f %f %f %f %f %f %f %f %f','delimiter','\t','BufSize',10000);
    catch objFoo
        errordlg(sprintf('Failed to parse STRING network file ''%s''\n\n Are you sure it''s a STRING network file?',strStringNetworkFile))
        error('%s: \t Failed to parse STRING network file ''%s''\n',mfilename,strStringNetworkFile)
        objFoo
    end
    
    if isempty(cellData{i,end})
        continue
    end
    
    
end

% remove optional empty rows from pre-initialization of cell array
matEmptyRows = all(cellfun(@isempty,cellData),2);
cellData(matEmptyRows,:) = [];

fclose('all');

% get rid of nested cells
cellData(cellfun(@iscell,cellData)) = cellfun(@(x) x{1}, cellData(cellfun(@iscell,cellData)),'uniformoutput',false);

% let's make one network, with three identifiers, only storing the combined
% string9 scores
matAllStringIdentifiers = cell2mat(cellData(:,[3,4]));
cellGeneSymbols = cellData(:,[1,2]);
cellENSPs = cellData(:,[5,6]);
matCombinedScores = cell2mat(cellData(:,15));

% look for unique string identifiers.
[a,b] = unique(matAllStringIdentifiers(:));
% also create corresponding gene symbol and ensp identifier lists
matStringIdentifiers = a;
cellStringGeneSymbols = cellGeneSymbols(b);
cellENSPs = cellENSPs(b);

matStringNetwork = zeros(numel(matStringIdentifiers));
% fill in network
for i = 1:size(cellData,1)
    iX1 = find(matStringIdentifiers==matAllStringIdentifiers(i,1));
    iX2 = find(matStringIdentifiers==matAllStringIdentifiers(i,2));
    matStringNetwork(iX1,iX2) = matCombinedScores(i);
    matStringNetwork(iX2,iX1) = matCombinedScores(i);
end


end