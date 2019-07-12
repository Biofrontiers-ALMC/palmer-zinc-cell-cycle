classdef Plate < handle
    %PLATE  Class to hold experimental data from an imaging plate
    %
    %  P = PLATE(dataDir) will import all suitable MAT-files in the data
    %  directory. Each MAT-file will be imported as a new Well object.
    
    properties (Access = private)

        Wells
        
    end
    
    properties (Dependent)
        
        NumWells
        FirstFrame
        LastFrame
        NumFrames
        TotalTracks
        
        DataCols
        
        NumTracks
    end
    
    properties (SetAccess = private)
        
        CellType
        RestingFRETinterval = [-8, -3];
        MaxFRETinterval = [1. 5];
        
        MitosisInterval = 5;
        
        CDK2PostMitosisWindow = [3, 8];    %In frames
        CDK2ratioThreshold = 0.5;
                
        CDK2PreMitosisWindow = [-5 -1];
        
        CDK2Window = [20 30];
        
    end
    
    methods
        
        function obj = Plate(varargin)
            %Constructor function
            
            if nargin > 0
                obj.importData(varargin{:});
            end
            
        end
        
        function importData(obj, varargin)
            %IMPORTDATA  Import tracked data from MAT-file
            %
            %  P.IMPORTDATA(dataDir) will import all suitable MAT-files
            %  from the directory supplied. Each MAT-file is treated as
            %  belonging to a different movie/well.
            %
            %  P.IMPORTDATA(S) where S is either a string or a character
            %  array will import the file(s) specified.
            %
            %  P.IMPORTDATA without any arguments will open a dialog box
            %  allowing you to specify which files to import.
            
            ip = inputParser;
            ip.addOptional('ImportLocation','',@(x) ischar(x) || iscellstr(x) );
            ip.parse(varargin{:});            
            
            if isempty(ip.Results.ImportLocation)
                %Open dialog box for user to select files to import
                
                [importFiles, importDir] = uigetfile({'*.mat','MAT-Files (*.mat)';...
                    '*.*','All Files (*.*)'}...
                    ,'Select files to import','MultiSelect','on');
                
                if ~iscell(importFiles)
                    if importFiles == 0
                        %Cancelled
                        return;
                    end
                end
                
            elseif iscell(ip.Results.ImportLocation)
                
                importFiles = ip.Results.ImportLocation;
                importDir = '';
                
            elseif isdir(ip.Results.ImportLocation)
                importDir = ip.Results.ImportLocation;
                
                %Get list of MAT files in the directory
                importFiles = dir(fullfile(importDir,'*.mat'));
                importFiles = {importFiles.name};
                
            elseif ischar(ip.Results.ImportLocation)
                
                importFiles = ip.Results.ImportLocation;
                importDir = '';
                
            end
            
            %Convert the filename to a cell if its not
            if ~iscell(importFiles)
                importFiles = {importFiles};
            end
            
            for iF = 1:numel(importFiles)
                fprintf('Importing file %s...\n',importFiles{iF});
               
                %Validate the file
                fileVarsInfo = who('-file',fullfile(importDir, importFiles{iF}));
                
                if ismember('trackData',fileVarsInfo)
                    
                    %Get well location from filename
                    wellLoc = regexp(importFiles{iF},'Well(?<loc>[a-zA-Z][0-9][0-9])','names');
                    
                    %Create a Well object from the file
                    obj.addWell(wellLoc.loc, fullfile(importDir, importFiles{iF}));
                    
                    %Load metadata (if not already populated)
                    if isempty(obj.CellType)
                        if ismember('metadata',fileVarsInfo)
                            
                            %Load the file into memory
                            load(fullfile(importDir, importFiles{iF}),'metadata')
                            
                            %Get cell type from the metadata
                            obj.CellType = metadata.CellType;
                            
                        end
                    end
                else
                    fprintf('\bNo track data\n');
                end
                
                fprintf('\bDONE\n');
            end
            
        end
        
        function szOut = size(obj, varargin)
            %SIZE  Get the number of rows x number of cols of wells
            
            szOut = size(obj.Wells, varargin{:});
            
        end
        
        function numWells = numel(obj)
            %NUMEL  Get number of elements in the plate
            %
            %  N = NUMEL(P) returns the number of well elements in the
            %  plate. The number of elements includes empty/missing wells.
            
            numWells = numel(obj.Wells);
            
        end
        
        function numPopulatedWells = get.NumWells(obj)
            %GET.NUMWELLS  Get number of wells with tracks
            
            numPopulatedWells = nnz(~cellfun(@isempty,obj.Wells));
            
        end
        
        function totalNumTracks = get.TotalTracks(obj)
            %GET.TOTALTRACKS  Get total number of tracks in the plate
            
            totalNumTracks = 0;
            for iW = 1:numel(obj.Wells)
                if ~isempty(obj.Wells{iW})
                    totalNumTracks = totalNumTracks + obj.Wells{iW}.NumTracks;
                end
            end
                    
        end
        
        function numFrames = get.NumFrames(obj)
            
            %Get a populated well
            wellIdx = find(~cellfun(@isempty,obj.Wells),1,'first');
            
            numFrames = obj.Wells{wellIdx}.NumFrames;          
            
        end
        
        function firstFrame = get.FirstFrame(obj)
            
            %Get a populated well
            wellIdx = find(~cellfun(@isempty,obj.Wells),1,'first');
            
            firstFrame = obj.Wells{wellIdx}.FirstFrame;      
            
        end
                
        function lastFrame = get.LastFrame(obj)
            
            %Get a populated well
            wellIdx = find(~cellfun(@isempty,obj.Wells),1,'first');
            
            lastFrame = obj.Wells{wellIdx}.LastFrame;      
            
        end
        
        function isPlateEmpty = isempty(obj)
            
            isPlateEmpty = numel(obj) == 0;
            
        end
        
        function dataCols = get.DataCols(obj)
            
            for iW = 1:numel(obj.Wells)
                if ~isempty(obj.Wells{iW})
                    dataCols = fieldnames(obj.Wells{iW}.getDataTable.Data);
                    return;
                end
            end
            
        end
        
        function numTracks = get.NumTracks(obj)
            
            numTracks = zeros(size(obj));
            %Get populated wells
            wellIdx = find(~cellfun(@isempty,obj.Wells));
            
            for iW = wellIdx'
                numTracks(iW) = obj.Wells{iW}.NumTracks;          
            end
            
        end
        
        %--- Well and track functions ---%
        
        function addWell(obj, wellLoc, input)
            %ADDWELL  Add a well to the plate
            %
            %  P.ADDWELL(L, A) inserts the TrackArray object to
            %  the well location L specified. L can either be a string
            %  specifiying the well location (e.g. 'A01') or a 1x2 vector
            %  specifiying the [row, col] subindices.
            
            %Parse well location
            if ischar(wellLoc)
                [rowNum, colNum] = obj.loc2ind(wellLoc);
                
            elseif isnumeric(wellLoc) && numel(wellLoc) == 2
                rowNum = wellLoc(1);
                colNum = wellLoc(2);            
                
            else
                error('Plate:addWell:InvalidWellLocation',...
                    'Expected well location to be a string (e.g. ''A01'') or [row, col] indices.')
            end
            
            obj.Wells{rowNum,colNum} = Well;
            
            %Set resting FRET interval and mitosis interval
            obj.Wells{rowNum,colNum}.RestingFRETinterval = obj.RestingFRETinterval;
            obj.Wells{rowNum,colNum}.MitosisInterval = obj.MitosisInterval;
            
            obj.Wells{rowNum,colNum} = obj.Wells{rowNum,colNum}.importData(input);
            
        end
        
        function trackArray = getWell(obj, wellLoc)
            %GETWELL  Get specified well
            %
            %  A = P.GETWELL(L) gets the well specified by L. L can either
            %  be a string with the well location (e.g. ''A01''), a
            %  subindex [row, col] vector (e.g. [2, 1]) or an index.
            
            if ischar(wellLoc)
                %Convert char to row/col
                trackArray = obj.Wells{obj.loc2ind(wellLoc)};
            elseif isnumeric(wellLoc)
                if numel(wellLoc) == 1
                    trackArray = obj.Wells{wellLoc};
                elseif numel(wellLoc) == 2
                    trackArray = obj.Wells{wellLoc(1), wellLoc(2)};
                end                
            end
            
        end
        
        function trackOut = getTrack(obj, wellLoc, trackInd)
            %GETTRACK  Get a track from a well
            %
            %  T = P.GETTRACK(L, I) gets track with index I from well L. L
            %  can either be a string containing the well location (e.g.
            %  ''A01'') or an index.
                        
            trackOut = obj.getWell(wellLoc).getTrack(trackInd);
                        
        end
        
        function dataTableOut = getDataTable(obj, wellLoc)
            %GETDATATABLE  Get data table from well
            %
            %  DT = P.GETDATATABLE(WellLocation)
            
            dataTableOut = obj.getWell(wellLoc).getDataTable;
            
        end
        
         %--- Plate-wide analysis functions---%
        
        function dataOut = queryData(obj, rows, cols, reqData, varargin)
             %QUERYDATA  Query selected Wells for data
             %
             %  D = P.QUERYDATA(R, C, P) will query the well(s) in the
             %  row(s) R and column(s) C specified for the data P.
             %
             %  D = P.QUERYDATA(R, C, P, F) will query the well(s) in the
             %  row(s) R and column(s) C specified for the data P,
             %  satisfying the optional filtering statements F
             
             dataOut = [];
             
             for iRow = rows
                 for iCol = cols
                     if ~isempty(obj.Wells{iRow, iCol})
                         dataOut = [dataOut; obj.Wells{iRow, iCol}.queryData(reqData, varargin{:})];
                     else
                        error('Plate:queryData:NonexistantWell',...
                            'Well in row %d, col %d does not exist.',...
                            iRow, iCol);                     
                     end
                 end
             end             
             
        end
        
        function plotTrack(obj, row, col, trackIdxs, reqData, varargin)
            %PLOTTRACK  Plots track(s) from a well
            %
            %  P.PLOTTRACK(R, C, I, P, options) will plot track(s) I from
            %  the well in row R and column C, with the specified options.
            
            obj.Wells{row, col}.plotTrack(trackIdxs, reqData, varargin{:});
            
        end
        
        function plotTrackByLoc(obj, wellLoc, trackIdxs, reqData, varargin)
            %PLOTTRACK  Plots track(s) from a well
            %
            %  P.PLOTTRACK(R, L, P, options) will plot track(s) I from
            %  the well in row R and column C, with the specified options.
            
            [row, col] = loc2ind(obj, wellLoc);
            
            obj.Wells{row, col}.plotTrack(trackIdxs, reqData, varargin{:});
            
        end
        
        function cellCnts = getCellCounts(obj,varargin)
            %GETCELLCOUNTS  Get cell counts over time
            %
            %  M = P.GETCELLCOUNTS will return a matrix M, with each
            %  row containing the cell counts for each well over time.
            %
            %  To only get cell counts from a specific number of wells,
            %  specify their location, M = P.GETCELLCOUNTS(wellLocation).
            %
            %  Example:
            %
            %   cellCnts = P.GETCELLCOUNTS('A06','A07'...).
            
            if isempty(varargin)
                %Initialize output struct
                cellCnts(size(obj, 1), size(obj, 2)).CellCounts = [];
                
                for iRow = 1:size(obj, 1)
                    for iCol = 1:size(obj, 2)
                        if ~isempty(obj.Wells{iRow, iCol})
                            cellCnts(iRow, iCol).CellCounts = obj.Wells{iRow, iCol}.getCellCounts;
                        end
                    end
                end
                
            else
                
                normalizeCounts = false;
                
                cellCnts = nan(numel(varargin), obj.NumFrames);
                for ii = 1:numel(varargin)
                    
                    if strcmpi(varargin{ii},'normalize')
                        normalizeCounts = true;
                        cellCnts(ii,:) = [];
                    else
                        if ~isempty(obj.getWell(varargin{ii}))
                            cellCnts(ii,:) = obj.getWell(varargin{ii}).getCellCounts;
                        end
                    end
                end
                
                if normalizeCounts
                    
                    for ii = 1:size(cellCnts,1)
                        cellCnts(ii,:) = cellCnts(ii,:) ./ cellCnts(ii,1);
                    end
                    
                end
                
            end
                        
        end
        
        function setRestingFRETinterval(obj, restingInterval)
            %SETRESTINGFRETINTERVAL  Set the resting FRET interval
            %
            %  P.SETRESTINGFRETINTERVAL([M, N]) will set the window for
            %  calculating the resting FRET to the interval [M, N]. Both M
            %  and N must be < 0 and M must be < N.
            
            for iW = 1:numel(obj.Wells)
                if ~isempty(obj.Wells{iW})
                    obj.Wells{iW}.RestingFRETinterval = restingInterval;                    
                end                
            end
            
            %If no errors, update the current resting FRET interval
            obj.RestingFRETinterval = restingInterval;
        end
        
        function setMaxFRETinterval(obj, maxFRETinterval)
            %SETRESTINGFRETINTERVAL  Set the resting FRET interval
            %
            %  P.SETRESTINGFRETINTERVAL([M, N]) will set the window for
            %  calculating the resting FRET to the interval [M, N]. Both M
            %  and N must be < 0 and M must be < N.
            
            for iW = 1:numel(obj.Wells)
                if ~isempty(obj.Wells{iW})
                    obj.Wells{iW}.MaxFRETinterval = maxFRETinterval;
                end
            end
            
            %If no errors, update the current resting FRET interval
            obj.MaxFRETinterval = maxFRETinterval;
        end
        
        function setMitosisinterval(obj, mitosisInterval)
            %SETRESTINGFRETINTERVAL  Set the resting FRET interval
            %
            %  P.SETRESTINGFRETINTERVAL([M, N]) will set the window for
            %  calculating the resting FRET to the interval [M, N]. Both M
            %  and N must be < 0 and M must be < N.
            
            for iW = 1:numel(obj.Wells)
                if ~isempty(obj.Wells{iW})
                    obj.Wells{iW}.MitosisInterval = mitosisInterval;                    
                end                
            end
            
            %If no errors, update the current resting FRET interval
            obj.MitosisInterval = mitosisInterval;
        end
        
        
        function setCDK2PostMitosis(obj, newInterval, newThreshold)
            %SETCDK2POSTMITOSIS  Set post-mitosis CDK2 levels
            %
            %  SETCDK2POSTMITOSIS(OBJ, INTERVAL, THRESHOLD) will set the
            %  post-mitosis CDK2 interval (window) and threshold to the new
            %  values specified. Changing these values will affect the
            %  MaxCDK2AfterMitosis and CDK2Classification values.
            %
            %  THRESHOLD is the CDK2 threshold level for the track to be
            %  classified as 'active' or 'emerging' (if the CDK2 level
            %  rises above the THRESHOLD value) or 'quiescent' (if the
            %  CDK2 level stays below the THRESHOLD value).
            %
            %  Example:
            %
            %  %Show current CDK2 parameters
            %  P.CDK2PostMitosisWindow
            %  P.CDK2ratioThreshold
            %
            %  %Change the window to be between 4 and 10 frames and ratio
            %  %to 0.5
            %  P.setCDK2PostMitosis([4 10], 0.5);
                        
            for iW = 1:numel(obj.Wells)
                if ~isempty(obj.Wells{iW})
                    obj.Wells{iW} = obj.Wells{iW}.setCDK2PostMitosis(newInterval, newThreshold);
                end
            end
            
            %If no errors, update the current CDK2 data
            obj.CDK2ratioThreshold = newThreshold;
            obj.CDK2PostMitosisWindow = newInterval;
            
        end
        
        function setCDK2PreMitosis(obj, newInterval)
            %SETCDK2PREMITOSIS  Set pre-mitosis CDK2 levels
            %
            %  SETCDK2PREMITOSIS(OBJ, INTERVAL) will set the pre-mitosis
            %  CDK2 interval (window) to the new value specified. Changing
            %  this window will affect the MaxCDK2BeforeMitosis.
            %
            %  Example:
            %
            %  %Show current CDK2 parameters
            %  P.CDK2PreMitosisWindow
            %
            %  %Change the window to be between -10 and -1 frames (note:
            %  %these values are negative, indicating the number of frames
            %  %before mitosis.)
            %  P.setCDK2PreMitosis([-10 -1]);
            
            for iW = 1:numel(obj.Wells)
                if ~isempty(obj.Wells{iW})
                   obj.Wells{iW} = obj.Wells{iW}.setCDK2PreMitosis(newInterval);
                end
            end
            
            %If no errors, update the current CDK2 data
            obj.CDK2PreMitosisWindow = newInterval;
            
        end
        
        function setCDK2Window(obj, newInterval)
            
            for iW = 1:numel(obj.Wells)
                if ~isempty(obj.Wells{iW})
                    obj.Wells{iW} = obj.Wells{iW}.setCDK2Window(newInterval);
                end
            end
            
            %If no errors, update the current CDK2 data
            obj.CDK2Window = newInterval;
        end
        
        
        function setDeltaT(obj, deltaT)
            %SETDELTAT  Set delta T for timestamp information
            %
            %  P.SETDELTAT(dT) will set the timestamp information of each
            %  Well in the Plate object to dT:(N*dT) where N is the number
            %  of frames in the Well.
            
            for iW = 1:numel(obj.Wells)
                if ~isempty(obj.Wells{iW})
                    obj.Wells{iW} = obj.Wells{iW}.setDeltaT(deltaT);
                end
            end
            
        end
        
        function tsOut = getTimestamps(obj, wellID)
            %GETTIMESTAMPS  Get timestamps from a specified well
            %
            %  TS = P.GETTIMESTAMPS(WellID) will get the timestamps vector
            %  from the well specified.
            %
            %  Example:
            %
            %   TS = P.GETTIMESTAMPS('A06') will return the timestamps from
            %   well A06 into the vector TS.           
            
            wellInd = obj.loc2ind(wellID);
            
            tsOut = obj.getWell(wellInd).getTimestamps;
            
        end
        
        function show(obj, varargin)
            %SHOW  Visual representation of plate data
            %
            %  SHOW(P) will plot a map showing existing wells as filled
            %  circles and empty wells as empty circles. The map will
            %  extend only as far as the size of the plate object (e.g. if
            %  the plate only has one row, only one row will be drawn).
            %
            %  SHOW(P, type) will show different maps depending on the
            %  value of 'type'.
            %
            %  For example:
            %
            %  %Show a map representing number of cell counts in each well
            %  show(obj, 'cell counts')
            %
            %  %Show a map representing number of existing wells (default)
            %  show(obj, 'map')
            
            %Parse the variable input to determine the type of plot to make
            if ~isempty(varargin)
                
                switch lower(varargin{1})
                    
                    case {'cell count','counts','proliferation','cellcounts'}
                        type = 'cellCount';

                        numTracks = zeros(size(obj));
                        for ii = 1:numel(obj)
                            
                            if ~isempty(obj.getWell(ii))
                                numTracks(ii) = obj.getWell(ii).NumTracks;
                            end
                            
                        end
                        maxNumTracks = max(numTracks(:));
                    
                    case {'mitosis events','events','divisions'}
                        type = 'mitosisEvents';
                        
                    case {'map','filled','existing', ''}
                        type = 'map';
                        
                    otherwise
                        error('Plate:show:InvalidType',...
                            '''%s'' is not a recognized map type. See help Plate.show for valid types.',varargin{1});                        
                        
                end
                
            else
                
                type = 'map';
                
            end
            
            
            %%%% Well map design:
            %   Each circle is a 64x64. Circle outer radius is 25 px.
            %   Inner radius is 23 px.
            
            headerSizePx = 30;  %Height and width of the column & row headers
           
            %Initialize the well map
            wellMap = zeros(64 * size(obj,1) + headerSizePx, 64 * size(obj,2) + headerSizePx, 'uint8');
                        
            %Draw the wells
            for iR = 1:size(obj,1)
                for iC = 1:size(obj,2)
                    
                    switch type
                        
                        case 'cell count'
                            
                            if isempty(obj.getWell([iR,iC]))
                                circOut = makeCircle(0,'heat');
                                
                            else
                                
                                circOut = makeCircle( (numel(obj.getWell([iR,iC])) ./ maxNumTracks)*255,'heat');
                                
                            end
                            
                        
%                         case 'mitosisEvents'
%                             
%                             if isempty(obj.getWell([iR,iC]))
%                                 circOut = makeCircle(0,'heat');
%                             else
%                                 circOut = makeCircle(numel(obj.MitosisEvents{iR,iC})./obj.TotalMitosisEvents * 255,'heat');
%                             end
                            
                        case 'map'
                            
                            %Determine the type of circle to draw based on whether
                            %the location exists
                            if ~isempty(obj.getWell([iR,iC]))
                                circOut = makeCircle(true,'flat');
                            else
                                circOut = makeCircle(false,'flat');
                            end
                            
                            circOut = circOut .* 255;
                        
                    end
                    
                    %Add the circle to the map
                    wellMap((iR - 1) *64 + headerSizePx + 1:iR *64 + headerSizePx,...
                        (iC - 1) * 64 + headerSizePx + 1: iC * 64 + headerSizePx) =  circOut;
                end
            end

            %Draw the column and row numbers
            for iC = 1:size(obj,2)
                wellMap = insertText(wellMap,[headerSizePx + 32 + (iC - 1)* 64, headerSizePx/2],...
                    iC,...
                    'AnchorPoint','center',...
                    'BoxOpacity',0,'TextColor','white','FontSize',18);
            end

            for iR = 1:size(obj,1)
                wellMap = insertText(wellMap,[headerSizePx/2,headerSizePx + 32 + (iR - 1)* 64],...
                    char(iR + 64),...
                    'AnchorPoint','center',...
                    'BoxOpacity',0,'TextColor','white','FontSize',18);
            end            
            
            %Plot the map
            imshow(wellMap,[])
            
            %--- Start subfunction ---%
            function circOut = makeCircle(fillValue, fillType)
                %Subfunction to draw filled/unfilled 64x64 circles
                
                if islogical(fillValue)
                    circOut = false(64);
                    fillType = 'flat';
                else
                    %Make a color image
                    circOut = zeros(64,'uint8');
                end
                
                xx = (1:64) - 32;
                [xx,yy] = meshgrid(xx,xx);
                
                switch fillType
                    
                    case 'flat'
                        
                        circOut(xx.^2 + yy.^2 <= 25.^2) = 1;
                        circOut(xx.^2 + yy.^2 <= 23.^2) = fillValue;
                        
                    case 'heat'
                        
                        circOut = fillValue .* exp(-(xx.^2 + yy.^2)/15.^2);
                                                
                        %circOut(xx.^2 + yy.^2 <= 23.^2) = fillValue;
                        
                        %Draw the circle around each well
                        %circOut(xx.^2 + yy.^2 <= 25.^2) = 255;
                                            
                end
                
                
                
            end
            %--- End subfunction ---%
            
        end
        
        
        function [trackData, timeVec] = getTrackData(obj, row, col, trackIdxs, trackProp)
            
            [trackData, timeVec] = obj.Wells{row,col}.getTrackData(trackIdxs, trackProp);
            trackData = cat(1,trackData{:});
            
        end
        
        %TODO- split track tidy up
        %Integrate the missed mitosis event detection
        
        function splitTrack(obj, row, col, trackIdx, splitFrame)
            %SPLITTRACK  Split a track
            %
            %  SPLITTRACK(OBJ, ROW, COL, TRACKID, FRAME) will split the
            %  track in the well specified. TRACKID can be a vector of
            %  track IDs. FRAME should be the same size as TRACKID and
            %  contain the frame number to split on. New tracks will be
            %  added to the end of the TrackArray.
            
            currWell = obj.Wells{row, col};
            
            %Split the specified tracks
            for ii = 1:numel(trackIdx)
            
                currWell.TrackArray = splitMitosis(currWell.TrackArray, trackIdx(ii), splitFrame(ii));
               
            end
            
            obj.Wells{row, col} = Well(currWell.TrackArray);
                        
        end
        
        function fixMitosis(obj, varargin)
            %FIXMITOSIS  Attempts to identify and split missing mitosis
            %events
            
            ip = inputParser;
            ip.addParameter('MinPeakProminence', 1000);
            ip.addParameter('MaxPeakWidth', 7);
            ip.addParameter('MinPeakDistance', 10);
            ip.addParameter('TestMode', false);
            parse(ip, varargin{:});
            
            
            
            for iWell = 1:numel(obj.Wells)
                
                %fprintf('Looking for missed mitosis events in well %s \n', obj.Wells{iWell}.Location);
                
                if ~isempty(obj.Wells{iWell})
                    
                    currWell = obj.Wells{iWell};
                    TA = getTrackArray(currWell);
                    
                    tracksToSplit = [];
                    frameToSplit = [];                    
                    
                    for iTrack = 1:TA.NumTracks
                        
                        %Search the nuclear intensity channel to find peaks
                        ct = getTrack(TA, iTrack);
                        data = getData(ct, 'NuclInt');
                        
                        dataSm = data;
                        %dataSm = smooth(data, 3);
                        
                        if numel(data) > 20
                            [pks, loc] = findpeaks(dataSm,...
                                'MinPeakProminence', ip.Results.MinPeakProminence,...
                                'SortStr', 'ascend',...
                                'MaxPeakWidth', ip.Results.MaxPeakWidth,...
                                'MinPeakDistance', ip.Results.MinPeakDistance,...
                                'Npeaks', 1);
                            
                            %Remove the peak location that is close to the supposed mitosis
                            %event so that we don't split the track where it already is
                            remIdx = abs(loc - numel(data)) <= 4;
                            loc(remIdx) = [];
                            pks(remIdx) = [];
                            
                            if numel(pks) > 0
                                
                                if ip.Results.TestMode
                                    %Plot remaining peaks to check
                                    tt = 1:numel(data);
                                    
                                    cdk2Data = getData(ct, 'CDK2ratioCorr');
                                    
                                    subplot(1,2,1)
                                    plot(tt, data, tt(loc), pks, 'x')
                                    ylabel('Nuclear intensity')
                                    
                                    subplot(1,2,2)
                                    plot(tt, cdk2Data, tt(loc), cdk2Data(loc), 'x')
                                    ylabel('CDK2 ratio')
                                    title('Press key for next plot')
                                    
                                    pause
                                    
                                else
                                    
                                    if ~isempty(loc)
                                        
                                        %Let's try splitting the track and adding it to a new array
                                        tracksToSplit = [tracksToSplit ones(1, numel(loc)) .* iTrack];
                                        
                                        %Split the track at the second mitosis point
                                        frameToSplit = [frameToSplit, loc' + ct.FirstFrame - 1];
                                        
                                    end
                                    
                                end
                                
                            end
                        end
                        
                    end
                    
                    %Split the specified tracks
                    for ii = 1:numel(tracksToSplit)
                        currWell.TrackArray = splitMitosis(currWell.TrackArray, tracksToSplit(ii), frameToSplit(ii));
                    end
                    
                    fprintf('Fixed %d tracks \n', numel(tracksToSplit));
                    
                    if numel(tracksToSplit) > 0
                        obj.Wells{iWell} = Well(currWell.TrackArray);
                    end
                end
                
            end
        end
        
        
    end
    
    methods (Access = private)
        
        function varargout = loc2ind(obj,wellLoc)
            %LOC2IND  Converts well location string to index
            %
            %  I = LOC2IND(S) converts the well location in string S into
            %  an index I.
            %
            %  [I, J] = LOC2IND(S) converts the well location into the row
            %  index I and col index J.
            
            if ~ischar(wellLoc)
                error('Plate:loc2ind:LocationNotChar',...
                    'Expected the well location to be a char.')                
            end
            
            if nargout == 1
                varargout{1} = sub2ind(size(obj),...
                    (wellLoc(1)) - 64,...
                    str2double(wellLoc(2:3)));
                
            elseif nargout == 2
                varargout{1} = double(wellLoc(1)) - 64;
                varargout{2} = str2double(wellLoc(2:3));
            end
        end
        
    end
    
end