classdef ZNSTrackArray
    %ZNSTrackArray  Data class for holding multiple tracks
    %
    %  ZNSTrackArray is a data class, designed to hold information for
    %  multiple tracks.
    
    properties (Access = private)
        Tracks
    end
    
    properties (Dependent)
        
        FirstFrame
        LastFrame
        NumFrames
        MeanFrameDeltaT
        
        NumDividingMotherTracks
        NumNonDividingTracks
                
        hasFRETdata
        hasCDK2data
        hasBGCorrectedData
        
    end
    
    properties (SetAccess = private, Hidden)  %Save file metadata
        SaveSerialDate = now;   %Serial date of when this file was created
        Timestamps
        TimestampUnit
    end
    
    methods
        
        function numTracks = size(obj)
            
            numTracks = [1 numel(obj)];
            
        end
        
        function numTracks = numel(obj)
            
            numTracks = numel(obj.Tracks);
            
        end
        
        function firstFrame = get.FirstFrame(obj)
            %GET.FIRSTFRAME  Gets the first frame in all tracks
            
            firstFrame = min([obj.Tracks.FirstFrame]);
            
        end
        
        function lastFrame = get.LastFrame(obj)
            %GET.LASTFRAME  Gets the last frame in all tracks
            
            lastFrame = max([obj.Tracks.LastFrame]);           
                        
        end
        
        function numFrames = get.NumFrames(obj)
            %GET.NUMFRAMES  Get number of frames in the track array
            
            numFrames = obj.LastFrame - obj.FirstFrame + 1;
            
        end
        
        function hasFRETdata = get.hasFRETdata(obj)
            %Check if tracks have FRET data
            
            hasFRETdata = false;
            if numel(obj) > 0
                if ismember('FRETInt',obj.getTrack(1).DataFieldnames)
                    hasFRETdata = true;
                end
            end
            
        end
        
        function hasCDK2data = get.hasCDK2data(obj)
            
            hasCDK2data = false;
            if numel(obj) > 0
                if ismember('CDK2nucl',obj.getTrack(1).DataFieldnames)
                    hasCDK2data = true;
                end
            end
            
        end
        
        function hasBGCorrectedData = get.hasBGCorrectedData(obj)
            
            hasBGCorrectedData = false;
            if numel(obj) > 0
                if ismember('FRETratioCorr',obj.getTrack(1).DataFieldnames) || ismember('CDK2ratioCorr',obj.getTrack(1).DataFieldnames)
                    hasBGCorrectedData = true;
                end
            end
            
        end
        
        function obj = setTimestamp(obj, timestampData, timestampUnit)
            %SETTIMESTAMP  Set the timestamp data
            
            obj.Timestamps = timestampData;
            obj.TimestampUnit = timestampUnit;
            
        end
        
        function [tsData, tsUnit] = getTimestamp(obj)
            %SETTIMESTAMP  Set the timestamp data
            
            tsData = obj.Timestamps;
            tsUnit = obj.TimestampUnit;
            
        end
        
        function deltaT = get.MeanFrameDeltaT(obj)
            %GET.FRAMEDELTAT  Get time between frames
            %
            %  A.FRAMEDELTAT will return the mean difference of the time
            %  between frames. This should provide at least an approximate
            %  value of the frame rate, although the mean could be
            %  incorrect if there is a stop between frames.
            
            deltaT = mean(diff(obj.getTimestamp));
            
        end
        
        %----- Track functions -----%
        
        function [obj, idxNewTrack] = addTrack(obj,tFrame,varargin)
            %ADDTRACK  Create a new track
            %
            %  A = A.ADDTRACK(tFrame, frameData) will create a new track
            %  which starts at the tFrame specified. The track will be
            %  populated with frameData.
            
            %Create the track
            if isempty(obj.Tracks)
                obj.Tracks = ZNSTrack;
                if ~isempty(varargin)
                    obj.Tracks = obj.Tracks.addFrame(tFrame, varargin{:});
                end
                idxNewTrack = 1;
            else
                idxNewTrack = numel(obj.Tracks) + 1;
                obj.Tracks(idxNewTrack) = ZNSTrack;
                if ~isempty(varargin)
                    obj.Tracks(idxNewTrack) = obj.Tracks(idxNewTrack).addFrame(tFrame, varargin{:});
                end
            end
            
        end
        
        function obj = deleteTrack(obj,trackIdx)
            %DELETETRACK  Delete a track
            %
            %  A = A.DELETETRCK(I) delete track I from the array.
            
            obj.Tracks(trackIdx) = [];
            
        end
        
        function trackOut = getTrack(obj,idxTrack)
            %GETTRACK  Get a track from the array
            %
            %  T = A.GETTRACK(trackIdx) returns the track T specified by
            %  the index specified.
            
            trackOut = obj.Tracks(idxTrack);            
            
        end
        
        function obj = updateTrack(obj,idxTrack,tFrame,varargin)
            %UPDATETRACK  Update a track by adding frame to the end
            %
            %  A = A.UPDATETRACK(trackIdx, tFrame, frameData) will add a
            %  frame to the end of the track specified. tFrame has to be
            %  larger than the last frame index of the track.
            
            if tFrame > obj.Tracks(idxTrack).LastFrame
                obj.Tracks(idxTrack) = obj.Tracks(idxTrack).addFrame(tFrame,varargin{:});                
            else
                error('ZNSTrackArray:updateTrack:InvalidFrame',...
                    'tFrame has to be larger than the last frame of track specified.');
            end
                        
        end
        
        function obj = updateTrackProperty(obj, trackIdx, propertyName, propertyValue)
                        
            obj.Tracks(trackIdx).(propertyName) = propertyValue;
            
        end
        
        function obj = deleteFrame(obj,idxTrack,tFrame)
            %DELETEFRAME  Delete a frame from a track
            %
            %  A = A.DELETEFRAME(trackIdx, tFrame) will delete frame tFrame
            %  from the track specified
            
            obj.Tracks(idxTrack) = obj.Tracks(idxTrack).deleteFrame(tFrame);
            
        end
        
        function [obj, idxAppended] = appendTracks(obj, varargin)
            %APPENDTRACKS  Add new tracks to the end of the track array
            %
            %  [A, I] = A.APPENDTRACKS(T1, T2, ..., TN) will add tracks T1,
            %  T2, ..., TN, to the end of the track array. The indices
            %  corresponding the to these tracks are returned in I.
            %
            %  [A, I] = A.APPENDTRACKS(A2) will add all tracks from another
            %  ZNSTrackArray object to the end of the track array. The
            %  indices corresponding the to these tracks are returned in I.
            
            idxAppended = [];
            for iArg = 1:numel(varargin)
                
                switch class(varargin{iArg})
                    
                    case 'ZNSTrackArray'
                        
                        idxNewTrack = zeros(1, numel(varargin{iArg}));
                        for iTrack = 1:numel(varargin{iArg})
                            
                            idxNewTrack(iTrack) = numel(obj.Tracks) + 1;
                            
                            obj.Tracks(idxNewTrack(iTrack)) = varargin{iArg}.getTrack(iTrack);
                            
                        end                        
                        
                    case 'ZNSTrack'
                        
                        idxNewTrack = numel(obj.Tracks) + 1;
                        obj.Tracks(idxNewTrack) = varargin{iArg};
                    
                end
                
                idxAppended = [idxAppended, idxNewTrack];
                
            end
            
            
        end
        
        %----- Analysis functions -----%
        
        function describe(obj, tFrame)
            %DESCRIBE  Prints useful information about the tracks
            %
            %  A.DESCRIBE will show the mean, min, max, and std of the
            %  measured track properties.
            
            if isempty(obj.Tracks)
                fprintf('No tracks available.');                
            end
            
            %Get list of track properties
            trackProps = obj.Tracks(1).getPropertyList;
            
            dataOut(numel(trackProps)) = struct('Property',[],...
                'Mean', [], 'StDev', [], 'Min', [], 'Max', []);
            
            %Collect data - should this be mean of means? min of mins? etc.
            %Or per frame specified?
            for iP = 1:numel(trackProps)
                dataOut(iP).Property = trackProps{iP};
                
                currPropValues = cat(1,obj.Tracks.Data(tFrame).(trackProps{iP}));
                
                dataOut(iP).Mean = mean(currPropValues,1,'omitnan');
                dataOut(iP).StDev = std(currPropValues,1,'omitnan');
                dataOut(iP).Min = min(currPropValues,[],1,'omitnan');
                dataOut(iP).Max = max(currPropValues,[],1,'omitnan');
                
            end
            
            tblOut = struct2table(dataOut);
            disp(tblOut);
            
        end
        
        function cellCnt = getCellCounts(obj)
            %GETCELLCOUNTS  Get cell counts over time
            %
            %  C = A.GETCELLCOUNTS will return the number of objects in
            %  each frame as a row vector C.
            
            cellCnt = zeros(1, obj.NumFrames);
            
            for ii = 1:numel(obj)
                
                %Add a cell count for each frame that an object exists
                %Note: Missing frames are still counted as an object
                currTrack = obj.getTrack(ii);
                cellCnt(currTrack.FirstFrame:currTrack.LastFrame) = ...
                    cellCnt(currTrack.FirstFrame:currTrack.LastFrame) + 1;
                
            end
            
        end
        
        function numDividingTracks = get.NumDividingMotherTracks(obj)
            %GET.NUMDIVISIONEVENTS  Get number of division events
            %
            %  The number of division events is equal to the number of
            %  tracks which have a populated 'DaugherTrackIdxs' property
            %  (i.e. despite having two daughter tracks, each mitosis event
            %  division counts only once).
            %
            %  Note: This is actually the number of mother tracks
            
            hasDaughter = false(1,numel(obj));
            for ii = 1:numel(obj)
                if ~isnan(obj.Tracks(ii).DaughterIdxs)
                    hasDaughter(ii) = true;
                end
            end
            
            numDividingTracks = nnz(hasDaughter);
        end
        
        function motherIdxList = getMotherTrackList(obj)
            %GETMOTHERTRACKLIST  Get a list of mother tracks
            %
            %  V = A.GETMOTHERTRACKLIST returns a vector containing a list
            %  of track indices of cells which divided. The indices are for
            %  the mother tracks (i.e. tracks which have a populated
            %  'DaughterTrackIdxs' field).
            
            %daughterIdxs = {obj.Tracks.DaughterIdxs};
            
            %motherIdxList = find(~isnan({obj.Tracks.DaughterIdxs}));
            
            hasDaughter = false(1,numel(obj));
            for ii = 1:numel(obj)
                if ~isnan(obj.Tracks(ii).DaughterIdxs)
                    hasDaughter(ii) = true;
                end
            end
            
            motherIdxList = find(hasDaughter);
        end
        
        function numNonDividingTracks = get.NumNonDividingTracks(obj)
            %GET.NUMNONDIVIDINGTRACKS  Get number of non-dividing tracks
            %
            %  A non-dividing track is a track that never divides in its
            %  lifetime.
            
            neverDivided = false(1,numel(obj));
            for ii = 1:numel(obj)
                if any(isnan(obj.Tracks(ii).DaughterIdxs)) && isnan(obj.Tracks(ii).MotherIdx)
                    neverDivided(ii) = true;
                end
            end
                       
            numNonDividingTracks = nnz(neverDivided);
            
        end
        
        function ndivCellList = getNonDividingCells(obj)
            %GETNONDIVIDINGCELLS  Get a list of non-dividing cells
            %
            %  V = A.GETNONDIVIDINGCELLS returns a vector containing a list
            %  of track indices of cells which never divided (i.e.
            %  motherIdx and daughterIdxs are nans).
            
            neverDivided = false(1,numel(obj));
            for ii = 1:numel(obj)
                if any(isnan(obj.Tracks(ii).DaughterIdxs)) && isnan(obj.Tracks(ii).MotherIdx)
                    neverDivided(ii) = true;
                end
            end
            
            ndivCellList = find(neverDivided);
        end
    end
    
end