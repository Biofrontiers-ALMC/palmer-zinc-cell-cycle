classdef Well
    %WELL  Class representing data from a single well
    %
    %  This class will import and analyze data from a well.
    
    properties %(Access = private)
        
        TrackArray
        DataTable = ZNSDataTable;
        
    end
    
    properties (Dependent)
       
        NumTracks
        NumFrames
        FirstFrame
        LastFrame
        NumRecords
        RestingFRETinterval
        MitosisInterval
        MaxFRETinterval
        
    end
    
    methods
        
        function obj = Well(varargin)
            %Constructor
            
            if ~isempty(varargin)
                
                obj = obj.importData(varargin{1});
                
            end
                        
        end
        
        function obj = importData(obj, varargin)
            %IMPORTDATA  Import data from MAT-file
            %
            %  W = W.IMPORTDATA(S) imports data from the MAT-file
            %  specified. If the filename is left empty, a dialog box will
            %  pop-up asking the user to select a file.
            
            if isempty(varargin)
                
                %Get a filename if none was specified
                [filename, pathname] = uigetfile({'*.mat','MAT-file (*.mat)'},...
                    'Select tracked data file');
                
                if ~ischar(filename)
                    %Cancelled
                    return;
                end
                
                filename = fullfile(pathname, filename);
            
            elseif ischar(varargin{1})
                
                filename = varargin{1};
                
            elseif isa(varargin{1},'ZNSTrackArray')
                
                trackData = varargin{1};
                   
            else
                 error('Well:importData:UnknownInput',...
                        'Expected input to be either a filename or a ZNSTrackArray object.')                    
            end
            
            %If a filename was specified, load the track data
            if exist('filename','var')
                %Check that this is a valid file
                load(filename, 'trackData');
                
                if ~exist('trackData','var')
                    error('Well:importData:TrackDataNotFound',...
                        'The MAT-file does not appear to contain track data.')
                end
            end
            
            %Populate the track array and data table properties
            obj.TrackArray = trackData;
            obj.DataTable = obj.DataTable.analyzeTrackArray(trackData);            
                        
        end
        
        function numTracks = get.NumTracks(obj)
            %GET.NUMTRACKS  Get number of tracks
            
            if isempty(obj.TrackArray)
                numTracks = NaN;
                return;
            end
            
            numTracks = numel(obj.TrackArray);
            
        end
        
        function numRecords = get.NumRecords(obj)
            %GET.NUMRECORDS  Get number of records in data table
            
            if isempty(obj.DataTable)
                numRecords = NaN;    
                return;
            end
         
            numRecords = numel(obj.DataTable);
            
        end
        
        function numFrames = get.NumFrames(obj)
            
            numFrames = obj.TrackArray.NumFrames;
            
        end
        
        function firstFrame = get.FirstFrame(obj)
            
            firstFrame = obj.TrackArray.FirstFrame;
            
        end
        
        function lastFrame = get.LastFrame(obj)
            
            lastFrame = obj.TrackArray.LastFrame;
            
        end
        
        function trackOut = getTrack(obj, trackIdx)
            %GETTRACK  Get specified track
            %
            %  T = W.GETTRACK(I) returns the ZNSTrack object with index i
            %  to T.
            
            trackOut = obj.TrackArray.getTrack(trackIdx);
            
        end
        
        function dataTableOut = getDataTable(obj)
            %GETDATATABLE  Returns the data table
            
            dataTableOut = obj.DataTable;
            
        end
        
        function trackArrayOut = getTrackArray(obj)
            %GETTRACKARRAY  Returns the TrackArray
            
            trackArrayOut = obj.TrackArray;
            
        end
        
        function [cellCnts, tVec] = getCellCounts(obj, varargin)
            %GETCELLCOUNTS  Gets cell counts
            %
            %  [C, T] = W.GETCELLCOUNTS will return the cell counts in C,
            %  and an (optional) vector of frames in T.
            
            cellCnts = obj.TrackArray.getCellCounts;
            
            if nargout > 1
                tVec = obj.TrackArray.FirstFrame:obj.TrackArray.LastFrame;
            end
            
            if ~isempty(varargin) && strcmpi(varargin{1}, 'normalize')
                cellCnts = cellCnts ./ (cellCnts(1));
                
            end
            
        end
        
        function [dataOut, recIdx] = queryData(obj, reqData, varargin)
            %QUERY  Send a query to the data table
            %
            %  D = W.QUERY(P)
            
            [dataOut, isValid] = obj.DataTable.get(reqData, varargin{:});
            
            if nargout == 2
                recIdx = find(isValid);
            end
        end
       
        function obj = setDeltaT(obj, deltaT)
            %SETDELTAT  Manually set the time between frames
            %
            %  W = W.SETDELTAT(dT) will set the time between frames to dT,
            %  where dT is in seconds.
            %
            %  
            
            %Regenerate a timestamp vector for the TrackArray
            tsVec = (1:obj.TrackArray.NumFrames) * deltaT;
            
            obj.TrackArray = obj.TrackArray.setTimestamp(tsVec,'s');
            
        end
        
        function tsOut = getTimestamps(obj)
            
            tsOut = obj.TrackArray.getTimestamp;
            
        end
        
        function restingFRETint = get.RestingFRETinterval(obj)
            
            restingFRETint = obj.DataTable.RestingFRETinterval;
            
        end
        
        function mitosisInt = get.MitosisInterval(obj)
            
            mitosisInt = obj.DataTable.MitosisInterval;
            
        end
        
        function maxFRETint = get.MaxFRETinterval(obj)
            
            maxFRETint = obj.DataTable.MaxFRETinterval;
        
        
        end
        
        function isWellEmpty = isempty(obj)

            isWellEmpty = isempty(obj.TrackArray);
            
        end
        
        function obj = set.RestingFRETinterval(obj, newInterval)
            
            if isempty(obj.DataTable)
                obj.DataTable = obj.DataTable.setRestingFRETinterval(newInterval);            
            else
                obj.DataTable = obj.DataTable.setRestingFRETinterval(newInterval, obj.TrackArray);
            end
            
        end
        
        function obj = set.MitosisInterval(obj, newInterval)
            
            if isempty(obj.DataTable)
                obj.DataTable = obj.DataTable.setMitosisInterval(newInterval);
            else
                obj.DataTable = obj.DataTable.setMitosisInterval(newInterval, obj.TrackArray);
            end
            
        end
        
        function obj = set.MaxFRETinterval(obj, newInterval)
            
            if isempty(obj.DataTable)
                obj.DataTable = obj.DataTable.setMaxFRETinterval(newInterval);
            else
                obj.DataTable = obj.DataTable.setMaxFRETinterval(newInterval, obj.TrackArray);
            end
            
        end
        
        function obj = setCDK2PostMitosis(obj, newInterval, newCDK2Threshold)

            if isempty(obj.DataTable)
                obj.DataTable = obj.DataTable.setCDK2PostMitosis(newInterval, newCDK2Threshold);
            else
                obj.DataTable = obj.DataTable.setCDK2PostMitosis(newInterval, newCDK2Threshold, obj.TrackArray);
            end
            
        end
        
        function obj = setCDK2PreMitosis(obj, newInterval)
            
            if isempty(obj.DataTable)
                obj.DataTable = obj.DataTable.setCDK2PreMitosis(newInterval);
            else
                obj.DataTable = obj.DataTable.setCDK2PreMitosis(newInterval, obj.TrackArray);
            end
            
        end
        
        function obj = setCDK2Window(obj, newInterval)
            
            if isempty(obj.DataTable)
                obj.DataTable = obj.DataTable.setCDK2Window(newInterval);
            else
                obj.DataTable = obj.DataTable.setCDK2Window(newInterval, obj.TrackArray);
            end
            
        end
        
    end
    
    methods
       
        function plotRecords(obj, recIdx, propToPlot, varargin)
            %PLOTRECORDS  Plot selected mitosis/non dividing tracks
            %
            %  PLOTRECORDS(W, R, P)
            %

            trackList = {obj.DataTable.getRecord(recIdx).TrackIdx};
            
            plotTrack(obj, trackList, propToPlot, varargin{:});
        end
        
        function plotTrack(obj, trackList, propToPlot, varargin)
            %PLOTTRACK  Plot a specified track
            %
            %  W.PLOTTRACK(I, P) will plot the property P from track with
            %  index I.
            %
            %  Examples:
            %
            %  %Plot NuclInt from Track 1
            %  W.PLOTTRACK(1, 'NuclInt')
            %
            %  If multiple
            %  W.PLOTTRACK(1, 'NuclInt', 'alignmitosis')
            %  W.PLOTTRACK(1, 'NuclInt', 'showmitosis')
                                    
            %Default settings
            timeUnit = 'Frames';
            timeMultiplier = 1;
            markMitosis = false;
            alignMitosis = false;
            dotColor = [1, 0, 0];
            smoothWindow = 0;   %0 = no smoothing
            normalizePlot = false;
            classifyCDK2 = false;
            
            
            %Parse variable input
            iV = 1;
            while iV <= numel(varargin)
                
                switch lower(varargin{iV})
                    
                    case 'alignmitosis'
                        
                        %Only adjust the time vector if there are two
                        %tracks
                        alignMitosis = true;
                        
                    case 'markmitosis'
                        
                        markMitosis = true;
                        
                    case 'seconds'
                        
                        timeUnit = 'Seconds';
                        timeMultiplier = obj.TrackArray.MeanFrameDeltaT;
                        
                    case 'minutes'
                        
                        timeMultiplier = obj.TrackArray.MeanFrameDeltaT / 60;
                        timeUnit = 'Minutes';
                        
                    case 'hours'
                        
                        timeMultiplier = obj.TrackArray.MeanFrameDeltaT / 3600;
                        timeUnit = 'Hours';
                        
                    case 'classifycdk2'
                        classifyCDK2 = true;
                        
                    case 'linecolor'
                        
                        if ~isnumeric(varargin{iV + 1}) && numel(varargin{iV + 1}) ~= 3
                            error('Invalid color specification. Expected 1x3 vector.')
                        end
                        
                        lineColor = varargin{iV + 1};
                        varargin(iV + 1) = [];
                        
                    case 'dotcolor'
                        
                        if ~isnumeric(varargin{iV + 1}) && numel(varargin{iV + 1}) ~= 3
                            error('Invalid color specification. Expected 1x3 vector.')
                        end
                        
                        dotColor = varargin{iV + 1};
                        varargin(iV + 1) = [];
                        
                    case 'smooth'
                        
                        if ~isnumeric(varargin{iV + 1})
                            error('Smooth window must be a single integer')
                        end
                        
                        smoothWindow = varargin{iV + 1};
                        varargin(iV + 1) = [];
                        
                    case 'normalize'
                        normalizePlot = true;
                        
                    case 'markmaxfret'
                        markMaxFRET = true;
                end
                
                iV = iV + 1;
            end            
            
            %Convert trackIdx into a cell (single entry inputs)
            if ischar(trackList) && strcmpi(trackList, 'all')
                trackList = {1:obj.NumTracks};
            elseif ~iscell(trackList)
                tempList = trackList;
                trackList = cell(1,numel(tempList));
                for ii = 1:numel(tempList)
                    trackList{ii} = tempList(ii);
                end
                clear tempList
            end
            
            trackList = fliplr(trackList);
            
            for iTr = 1:numel(trackList)
 
                %Get track data. trackData and frameVec will be cells
                [trackData, frameVec] = obj.getTrackData(trackList{iTr}, propToPlot);
                
                %Convert the frame vector into a time vector
                timeVec = cell2mat(frameVec);
                
                %Join the track data into a vector
                joinedTrackData = cat(1,trackData{:});
                
                %Smooth the data if requested
                if smoothWindow > 0
                    joinedTrackData = smooth(joinedTrackData,smoothWindow);
                end
                
                %Adjust the time vector if aligned to mitosis
                if alignMitosis && numel(frameVec) == 2
                    timeVec = timeVec - frameVec{1}(end);
                end
                
                %Convert the time vector to specified units
                timeVec = timeVec .* timeMultiplier;
                
                %Normalize the plots if set
                if normalizePlot 
                    switch propToPlot
                        case 'FRETratioCorr'
                            restingFRET = obj.DataTable.getByTrackIdx(trackList{iTr},'RestingFRETCorr');
                            joinedTrackData = joinedTrackData ./ restingFRET;
                        case 'FRETratio'
                            restingFRET = obj.DataTable.getByTrackIdx(trackList{iTr},'RestingFRET');
                            joinedTrackData = joinedTrackData ./ restingFRET;
                            
                        case {'CDK2ratioCorr', 'CDK2ratio', 'NuclInt'}
                            
                            joinedTrackData = joinedTrackData - prctile(joinedTrackData, 5);
                            
                        
                    end
                end
                
                if exist('lineColor','var')
                    plot(timeVec,joinedTrackData,'Color',lineColor);
                elseif classifyCDK2 && strcmpi(propToPlot, 'CDK2ratioCorr')
                    
                    cdk2Class = obj.DataTable.getByTrackIdx(trackList{iTr},'CDK2Classification');
                    
                    if cdk2Class > 0
                        switch cdk2Class
                            case 1 %Active
                                cdk2lineColor = [0 0 1];
                                
                            case 2 %Emerging
                                cdk2lineColor = [0 1 0];
                                
                            case 3 %Quiescent
                                cdk2lineColor = [1 0 0];

                        end
                        plot(timeVec,joinedTrackData,'Color',cdk2lineColor);
                    end
                    
                else
                    plot(timeVec,joinedTrackData);
                end
                
                hold on
                if markMitosis && numel(frameVec) == 2
                    plot(timeVec(numel(trackData{1})),joinedTrackData(numel(trackData{1})),'o',...
                        'MarkerFaceColor', dotColor, ...
                        'MarkerEdgeColor', dotColor);
                end
                               
            end
            hold off
            
            ylabel(propToPlot)
            xlabel(timeUnit)
            
        end
        
        function plotCellCount(obj, varargin)
            %PLOTCELLCOUNT  Plot the cell counts
            %
            %  W.PLOTCELLCOUNTS will plot the cell counts over time (in
            %  frames). 
            %
            %  You can specify the time units as an optional argument. For
            %  example, to plot the time in hours:
            %
            %  W.PLOTCELLCOUNTS('hours')
            
            tt = obj.TrackArray.FirstFrame:obj.TrackArray.LastFrame;
            timeUnit = 'Frames';
            
            if ~isempty(varargin)
                
                switch lower(varargin{1})
                    
                    case 'seconds'
                        
                        timeUnit = 'Seconds';
                        tt = tt .* obj.TrackArray.MeanFrameDeltaT;
                        
                    case 'minutes'
                        
                        tt = tt .* obj.TrackArray.MeanFrameDeltaT / 60;
                        timeUnit = 'Minutes';
                        
                    case 'hours'
                        
                        tt = tt .* obj.TrackArray.MeanFrameDeltaT / 3600;
                        timeUnit = 'Hours';
                    
                end
                
            end
            
            cellCnt = obj.TrackArray.getCellCounts;
            
            plot(tt,cellCnt);
            xlabel(timeUnit)
            ylabel('Cell Count')
            
        end
        
        function [trackData, frameVec] = getTrackData(obj, trackIdx, prop)
            %GETTRACKDATA  Get specified track data
            %
            %  D = W.GETTRACKDATA(I, P) gets the specified track data
            %
            %  [D, T] = W.GETTRACKDATA(I, P) gets the specified track data
            %  and the time information
            %
            %  D = W.GETTRACKDATA(I, P) if I is multiple will combine the
            %  data
            
            tracks = obj.TrackArray.getTrack(trackIdx);
            
            trackData = cell(1, numel(tracks));
            frameVec = cell(1, numel(tracks));
            
            for iTrack = 1:numel(tracks)
                trackData{iTrack} = tracks(iTrack).getData(prop);
                frameVec{iTrack} = tracks(iTrack).FirstFrame:tracks(iTrack).LastFrame;
            end
            
        end
        
    end
    
end