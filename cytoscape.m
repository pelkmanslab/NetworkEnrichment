function cytoscape(matGraph,cellLabels,varargin)
% help cytoscape
%
% usage:
%
% cytoscape(matGraph, cellLabels, ...)
%
%   cytoscape(matGraph) loads the network matGraph in cytoscape. matGraph
%   must be a square matrix. Values will be imported as edge scores. 0 or
%   NaN will be considered the absence of an edge. Values will be imported
%   into cytoscape as edge interactionStrength attribute.
%
%   cytoscape(matGraph,cellLabels) uses the labels defined by cellLabels.
%   The length of cellLabels must be equal to both the number of rows and
%   columns in matGraph.
%
%   cytoscape(matGraph, cellLabels, ...) loads additional attributes.
%   If the variable input arguments has the shape of matGraph it will be
%   imported as edge attribute. If the variable input argument has the
%   shape of cellLabels, it will be imported as node attribute. Attribute
%   names will follow the corresponding variable names.
%
%   You can also pass the full path to a ".props" file, which
%   specifies vizmap properties to load.
%
%   run cytoscape() to see a demo.
%
%   NOTE: cytoscape() contains a hardcoded reference to the executable of
%   cytoscape on your machine. Edit this reference to match the setup of
%   your machine.
%
%   Cytoscape will only return control back-over to MATLAB if the
%   executable is quit, which is the default behaviour of a "system" call.
%   Edit the 'cytoscape.m' code to see an alternative solution.
%
%   example:
%
%
%      matGraph = [ 0, 1, 1; ...
%                   0, 0, 1; ...
%                   0, 0, 1; ]
%      cellLabels = {'Gene A', 'Gene B', 'Gene C'}
%      cellENSPs = {'ENSP0001', 'ENSP0002', 'ENSP0003'}
%      cytoscape(matGraph, cellLabels, cellENSPs)
%
%
%   Copyright of Berend Snijder, 2011.


    if nargin==0
        % demo mode
        fprintf('%s: demo mode with almost random graph of 100 nodes\n',mfilename)
        matGraph = sign(rand(100)-0.5);
        matGraph = matGraph + (magic(100) > 900);
        matGraph(rand(100)<0.98)=0;
    else
        fprintf('\n%s: starting...\n',mfilename)
    end

    if nargin<=1 || isempty(cellLabels)
        % default node labels, if none are passed
        cellLabels = strcat('node_',arrayfun(@(x) num2str(x),1:size(matGraph,1),'UniformOutput',false));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Figure out the path to the cytoscape executable:

    % Default path to the cytoscape executable called by this function.
    % adjust this to your needs. The value stored here will be tested here
    % and only if invalid will the cfg-file-method be used.
    if ispc
        % example PC executable location
        strCytoscapeExecutable = 'C:\Program Files\Cytoscape_v2.8.2\Cytoscape.exe';
    elseif ismac
        % example mac executable location
        strCytoscapeExecutable = '/Applications/Cytoscape_v2.8.2/cytoscape.sh';
    else
        strCytoscapeExecutable = 'null';
    end
    
    % By default we will store the location of the cytoscape executable in
    % a txt file. If this txt file is not present, ask the user to locate
    % the cytoscape executable and go from there.
    strCfgPath = fullfile(fileparts(which(mfilename)),'cytoscape.txt');
    
    % if the cfg file does not exist, ask the user for the cytoscape
    % executable location and write that path to the cfg file.
    if ~fileattrib(strCytoscapeExecutable)
        if ~fileattrib(strCfgPath)
            strCurrentPath=cd;
            if ispc
                if fileattrib('C:\Program Files\'); cd('C:\Program Files\'); end
                [filename, pathname] = uigetfile({'Cytoscape.exe','Cytoscape executable (Cytoscape.exe)'}, 'Please locate your Cytoscape.exe file');
            else
                % cytoscape.sh is the executable for both mac and linux.
                if fileattrib('/Applications/'); cd('/Applications/'); end
                [filename, pathname] = uigetfile({'cytoscape.sh','Cytoscape executable (cytoscape.sh)'}, 'Please locate your cytoscape.sh file');
            end
            cd(strCurrentPath);

            % quit if user did not supply the path to a cytoscape executable
            if isequal(filename,0)
               errordlg(sprintf('This function cannot run unless you locate your Cytoscape executable.\n\nIf you do not have Cytoscape installed, you can download it for free from http://www.cytoscape.org.'),'Cytoscape is required')
               return
            end

            % store user selected path to the cytoscape executable
            fid = fopen(strCfgPath,'w');
            if fid==0
                warndlg(sprintf('Could not store the path to your cytoscape executable in ''%s''. This means next time we will have to ask you again. Alternatively edit the cytoscape.m file directly.',strCfgPath),'Warning')
            else
                strCytoscapeExecutable = fullfile(pathname,filename);
                fprintf(fid,'%s',strCytoscapeExecutable);
                fclose(fid);
            end

        elseif fileattrib(strCfgPath) && ~fileattrib(strCytoscapeExecutable)
            % read the cytoscape.cfg file
            fprintf('%s: reading cytoscape path from cfg file ''%s''\n',mfilename,strCfgPath)
            fid = fopen(strCfgPath);
            strCytoscapeExecutable = fgetl(fid);
            fclose(fid);
            strCytoscapeExecutable = strtrim(strCytoscapeExecutable);
            
            % if the stored file contains a wrong path, delete cfg file...
            % at least next time user can correct the mistake.
            if ~fileattrib(strCytoscapeExecutable)
                delete(strCfgPath);
                errordlg(sprintf('The stored path to the cytoscape executable (''%s'') is incorrect. Please restart the cytoscape() matlab function.',strCytoscapeExecutable),'Aborting')
                return
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % store original labels and pass as canonicalName
    cellOrigLabels = cellLabels;

    % replace all spaces and equal signs with underscores in node labels
    cellLabels = strrep(cellLabels,' ','_');
    cellLabels = strrep(cellLabels,'=','_');
    cellLabels = strrep(cellLabels,'''','_');    
    cellLabels = strrep(cellLabels,'+','_');
    cellLabels = strrep(cellLabels,',','_');

    % let's also say NaNs are no interactions.
    matGraph(isnan(matGraph)) = 0;

    % find interactions
    [iX,jX] = find(matGraph);

    % create temporary network file (SIF format)
    strFileName = sprintf('%s.sif',tempname);
    % open file for writing
    fid = fopen(strFileName,'w');
    % for each interaction in biograph, write new line in SIF format
    for i = 1:length(iX)
        % transform into SIF format
        strLine = sprintf('%s interacts %s\n',cellLabels{iX(i)},cellLabels{jX(i)});
        % write line
        fwrite(fid, strLine, 'char'); % write line
    end
    % add each individual node, to make sure unconnected nodes are present
    % in graph
    for i = 1:length(cellLabels)
        % transform into SIF format
        strLine = sprintf('%s\n',cellLabels{i});
        % write line
        fwrite(fid, strLine, 'char'); % write line
    end    
    % close file
    fclose(fid);

    
    %%%
    % we should add the matGraph values as interaction strength edge
    % attribute
    cellAttributeFiles = {};% use for deleting temp attribute files.
    strAttributeCommand = '';%#ok<NASGU> % append to system call
    
    % edge attribute detected
    strAttrFileName = sprintf('%s.eda',tempname);
    % add as edge attribute file to command line
    strAttributeCommand = sprintf('-e %s',strAttrFileName);
    % open file for writing
    fidAttr = fopen(strAttrFileName,'w');
    % format header with attribute description
    strHeaderLine = sprintf('%s (class=%s)\n','interactionStrength','Double');
    % write header line
    fwrite(fidAttr, strHeaderLine, 'char'); % write header line
    % for each interaction in biograph, write new line in SIF format
    for ii = 1:length(iX)
        % transform into SIF format
        strLine = sprintf('%s (interacts) %s = %.2g\n',cellLabels{iX(ii)},cellLabels{jX(ii)},matGraph(iX(ii),jX(ii)));
        % write line
        fwrite(fidAttr, strLine, 'char'); % write line
    end
    % close file
    fclose(fidAttr); 
    % add to list for cleanup from file system
    cellAttributeFiles(end+1) = {strAttrFileName}; %#ok<AGROW>
    %%%
    
    
    %%%
    % Add original labels as canonicalName by default
    strAttrFileName = sprintf('%s.noa',tempname);
    % add as node attribute file to command line
    strAttributeCommand = sprintf('%s -n %s',strAttributeCommand,strAttrFileName);
    % open file for writing
    fidAttr = fopen(strAttrFileName,'w');
    % write header line
    fwrite(fidAttr, sprintf('canonicalName (class=String)\n'), 'char'); % write header line
    % for each interaction in biograph, write new line in SIF format
    for ii = 1:length(cellOrigLabels)
        % transform into SIF format
        strLine = sprintf('%s = %s\n',cellLabels{ii},cellOrigLabels{ii});
        % write line
        fwrite(fidAttr, strLine, 'char'); % write line
    end
    % close file
    fclose(fidAttr);
%     go(strAttrFileName)
    % add to list for cleanup from file system
    cellAttributeFiles(end+1) = {strAttrFileName}; %#ok<AGROW>        
    %%%
    
    
    
    
    
    
    %%%
    % other attributes that were input via varargin:
    % "Node attribute file names often use the suffix ".noa", while edge
    % attribute file names use the suffix ".eda"."
    for i = 1:length(varargin)
        % if varargin contains stuff with equal size to GRAPH or LABELS,
        % create node or edge attribute file
        
        attrData = varargin{i};
        
        % name attribute according to input variable name, or default to
        % "Attribute#", where # is input number.
        strAttributeName = inputname(i+2);
        if isempty(strAttributeName)
            strAttributeName = sprintf('Attribute%d',i);
        end
        strClassName = class(varargin{i});
        % format header with attribute description
        strHeaderLine = sprintf('%s (class=%s)\n',strAttributeName,strClassName);        
        
        fprintf('%s: processing attribute: %s',mfilename,strHeaderLine)
        
        if isequal(size(attrData),size(matGraph))
            fprintf('%s: -- edge attribute\n',mfilename)
            % edge attribute detected
            strAttrFileName = sprintf('%s.eda',tempname);
            % add as edge attribute file to command line
            strAttributeCommand = sprintf('%s -e "%s"',strAttributeCommand,strAttrFileName);
            % open file for writing
            fidAttr = fopen(strAttrFileName,'w');
            % write header line
            fwrite(fidAttr, strHeaderLine, 'char'); % write header line
            % for each interaction in biograph, write new line in SIF format
            for ii = 1:length(iX)
                % transform into SIF format
                if isnumeric(attrData)
                    if isequal(round(attrData),attrData)
                        strLine = sprintf('%s (interacts) %s = %d\n',cellLabels{iX(ii)},cellLabels{jX(ii)},attrData(iX(ii),jX(ii)));
                    else
                        strLine = sprintf('%s (interacts) %s = %.2g\n',cellLabels{iX(ii)},cellLabels{jX(ii)},attrData(iX(ii),jX(ii)));
                    end
                else
                    strLine = sprintf('%s (interacts) %s = %s\n',cellLabels{iX(ii)},cellLabels{jX(ii)},attrData(iX(ii),jX(ii)));
                end
                % write line
                fwrite(fidAttr, strLine, 'char'); % write line
            end
            % close file
            fclose(fidAttr);            
            % add to list for cleanup from file system
            cellAttributeFiles(end+1) = {strAttrFileName}; %#ok<AGROW>
            
        elseif (isequal(size(attrData),size(cellLabels)) || isequal(size(attrData'),size(cellLabels))) && ~ischar(attrData)
            fprintf('%s: -- node attribute\n',mfilename)
            % node attribute detected
            strAttrFileName = sprintf('%s.noa',tempname);
            % add as node attribute file to command line
            strAttributeCommand = sprintf('%s -n "%s"',strAttributeCommand,strAttrFileName);
            % open file for writing
            fidAttr = fopen(strAttrFileName,'w');
            % write header line
            fwrite(fidAttr, strHeaderLine, 'char'); % write header line
            % for each interaction in biograph, write new line in SIF format
            for ii = 1:length(attrData)
                if iscell(attrData)
                    tempAttrData = attrData{ii};
                else
                    tempAttrData = attrData(ii);
                end
                % transform into SIF format
                if isnumeric(tempAttrData)
                    if isequal(round(tempAttrData),tempAttrData)
                        strLine = sprintf('%s = %d\n',cellLabels{ii},tempAttrData);
                    else
                        strLine = sprintf('%s = %.2f\n',cellLabels{ii},tempAttrData);
                    end
                else
                    strLine = sprintf('%s = %s\n',cellLabels{ii},tempAttrData);
                end
                % write line
                fwrite(fidAttr, strLine, 'char'); % write line
            end
            % close file
            fclose(fidAttr);            
            % add to list for cleanup from file system
            cellAttributeFiles(end+1) = {strAttrFileName}; %#ok<AGROW>     
        elseif ischar(attrData) && ~isempty(findstr(lower(attrData),'.props')) && fileattrib(attrData)
            fprintf('%s: -- vizmap file "%s"\n',mfilename,attrData)
            strAttributeCommand = sprintf('%s -V "%s"',strAttributeCommand,attrData);
        else
            warning('BS:Bla','unknown attribute type for input ''%s'' (%s). size must equal graph or label input size',strAttributeName,strClassName)
            continue
        end
        %%%
    
    end
    
    
    
    
    % check if it exists
    if ~fileattrib(strFileName)
        error('temporary network file not found')
    end
    
    % call cytoscape with network file (machine specific, change on your
    % setup!)
    fprintf('%s: running %s',mfilename,sprintf('"%s" -N "%s" %s\n',strCytoscapeExecutable,strFileName,strAttributeCommand))

    %%%%%%%%%%%%%%
    %%% First strategy: don't pass control back to matlab until cytoscape
    %%% is quit
    % note, the ampersand (&) forces control back to MATLAB immediately
    fprintf('\n\tNOTE: MATLAB WILL NOT RETURN CONTROL UNTIL YOU CLOSE CYTOSCAPE. \n\t\t  (You can change this behaviour in the cytoscape.m file)\n')    
    boolFailed = system(sprintf('"%s" -N "%s" %s',strCytoscapeExecutable,strFileName,strAttributeCommand));
    % if application call failed, remove stored application path..
    if boolFailed~=0
        disp('system call failed, cleaning up reference to cytoscape executable.')
        delete(strCfgPath)
    end
    %%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%
    %%% Alternative strategy: pass control back after 30 seconds
%     % note, the ampersand (&) forces control back to MATLAB immediately
%     system(sprintf('"%s" -N "%s" %s &',strCytoscapeExecutable,strFileName,strAttributeCommand));
%     
%     % Give cytoscape 30 seconds to load all the files, then clean up the
%     % temporary files and continue processing in matlab
%     pause(30)
    %%%%%%%%%%%%%%

    % when cytoscape is closed, delete temporary network file
    delete(strFileName);
    
    % delete all created temporary node and edge attribute files
    if ~isempty(cellAttributeFiles)
        cellfun(@delete,cellAttributeFiles);
    end
end