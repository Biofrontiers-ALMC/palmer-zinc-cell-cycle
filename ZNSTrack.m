classdef ZNSTrack
    %ZNSTRACK  Track data for the Zinc Signalling project
    %
    %  T = ZNSTrack(tFrame, frameData) will create a ZNSTrack object T with
    %  the first frame containing the data provided. tFrame should be the
    %  number of the frame.
    
    properties (Hidden)
        Data
        FrameIndex
    end
    
    properties
        MotherIdx = NaN;
        DaughterIdxs = NaN;
    end
    
    properties (Dependent)
       
        hasFRETdata
        hasCDK2data
        
        FirstFrame
        LastFrame
        NumFrames
        
        DataFieldnames
        
    end    
    
    methods
        
        function obj = ZNSTrack(varargin)
            
            if nargin > 0
                
                ip = inputParser;
                ip.addRequired('frameIndex', @(x) isnumeric(x) && isscalar(x));
                ip.addRequired('trackData', @(x) isstruct(x));
                ip.parse(varargin{:});
                
                obj = obj.addFrame(ip.Results.frameIndex, ip.Results.trackData);
                
            end
            
        end
        
        function firstFrame = get.FirstFrame(obj)
            %GET.FIRSTFRAME  Get the first frame index of the track
            %
            %  I = T.FIRSTFRAME will return the index of the first frame of
            %  the track. If the track is empty, I will be Inf.
            
            if isempty(obj.Data)
                firstFrame = Inf;            
            else
                firstFrame = obj.FrameIndex(1);
            end
            
        end
        
        function lastFrame = get.LastFrame(obj)
            %GET.LASTFRAME  Get the last frame index of the track
            %
            %  I = T.LASTFRAME will return the index of the first frame of
            %  the track. If the track is empty, I will be -Inf.
            
            if isempty(obj.Data)
                lastFrame = -Inf;               
            else
                lastFrame = obj.FrameIndex(end);                
            end
            
        end
        
        function numFrames = get.NumFrames(obj)
            %GET.NUMFRAMES  Returns the length of the track in frames
            %
            %  L = T.NUMFRAMES will return the length of the track in
            %  frames L.
            
            if isempty(obj.Data)
                numFrames = 0;                
            else            
                numFrames = obj.LastFrame - obj.FirstFrame + 1;
            end
            
        end
        
        function hasFRET = get.hasFRETdata(obj)
            %GET.HASFRETDATA  Returns true is the object has FRET data
            %
            %  B = T.HASFRETDATA returns a value of B = true if FRET data
            %  is present. Otherwise, B = false;
            
            if obj.NumFrames > 0
                
                %Look for FRET data fields in the Data property
                if ismember('FRETInt',obj.getPropertyList) ||...
                        ismember('FRETIntCorr',obj.getPropertyList)
                    hasFRET = true;
                    return
                end
                
            end
                
            hasFRET = false;                
            
        end
        
        function hasCDK2 = get.hasCDK2data(obj)
            %GET.HASCDK2DATA  Returns true is the object has CDK2 data
            %
            %  B = T.HASFRETDATA returns a value of B = true if FRET data
            %  is present. Otherwise, B = false;
            
            if obj.NumFrames > 0
                
                %Look for CDK2 data fields in the Data property
                if ismember('CDK2nucl',obj.getPropertyList) ||...
                    ismember('CDK2nuclCorr',obj.getPropertyList)
                    hasCDK2 = true;
                    return
                end
                
            end
            
            hasCDK2 = false;
            
        end
        
        function fnOut = get.DataFieldnames(obj)
            %GET.DATAFIELDNAMES  Get names of data fields
            
            if obj.NumFrames == 0
                fnOut = '';
            else
                fnOut = fieldnames(obj.Data(1));            
            end
                        
        end
        
        
        function obj = addFrame(obj, tFrame, frameData)
            %ADDFRAME  Add data for a frame
            %
            %  T = T.ADDFRAME(f, dataStruct) adds a new frame at index f to
            %  the start or the end of the track. The frame data should be
            %  in a structure, with the fieldnames of the structure
            %  corresponding to the measured data property name.
            %
            %  If the new frame data has a new property that was not
            %  present in the previous frames, the value for the missing
            %  data will be empty ([]).
            %
            %  Example:
            %
            %    T = ZNSTrack(1, struct('Area', 5));
            %
            %    %In frame 2, 'Area' is no longer measured, but 'Centroid'
            %    %is
            %    T = T.ADDFRAME(2, struct('Centroid', [10 20]));
            %
            %    %These are the expected outputs:
            %    T.Data(1).Area = 2
            %    T.Data(1).Centroid = []
            %
            %    T.Data(2).Area = []
            %    T.Data(2).Centroid = [10 20]
            %
            %  See also: ZNSTrack.updateTrack
            
            %Validate the frame number
            if ~isnumeric(tFrame)
                error('ZNSTrack:addFrame:frameIndexNotNumeric',...
                    'Expected the frame index to be a number.');
                
            elseif ~isscalar(tFrame)
                error('ZNSTrack:addFrame:frameIndexNotScalar',...
                    'Expected the frame index to be a scalar number.');
                
            else
                if ~(tFrame < obj.FirstFrame || tFrame > obj.LastFrame)
                    
                    error('ZNSTrack:addFrame:frameIndexInvalid',...
                        'The frame index should be < %d or > %d.',...
                        obj.FirstFrame, obj.LastFrame);
                end
            end
            
            %Valide the input data
            if ~isstruct(frameData)
                error('ZNSTrack:addFrame:dataNotStruct',...
                    'Expected data to be a struct.');
            end
            
            %Add the frame to the track
            if tFrame > obj.LastFrame
                
                if isinf(obj.FirstFrame) && isinf(obj.LastFrame)
                    %If both start and end frames are infinite, then this
                    %is the first frame to be added
                    obj.Data = frameData;
                    obj.FrameIndex = tFrame;
                    
                else
                    %Calculate the number of frames to add
                    numFramesToAdd = tFrame - obj.LastFrame;
                    
                    %Add the frame to the end of the array
                    obj.Data(end + numFramesToAdd) = frameData;
                    
                    %Update the frame indices
                    obj.FrameIndex = obj.FirstFrame:tFrame;
                    
                end
                
            elseif tFrame < obj.FirstFrame
                %Overwrite the Data property with new frame data, then move
                %the old data to the end of the structure.
                oldData = obj.Data;         %Save a copy of the old data
                obj.Data = frameData;       %Overwrite the Data property
                
                %Move the old data to the end of the structure
                dataInd = obj.FirstFrame - tFrame + 1;
                obj.Data(dataInd:dataInd + numel(oldData) - 1) = oldData;
                
                %Update the frame indices
                obj.FrameIndex = tFrame:obj.LastFrame;
                
            end
           
        end
        
        function obj = deleteFrame(obj, frameIndex)
            %DELETEFRAME  Deletes the specified frame
            %
            %  T = T.DELETEFRAME(f, frameIndex) deletes the specified
            %  frame(s) from the track.
            %
            %  Examples:
            %
            %    %Create a track with four frames
            %    trackObj = ZNSTrack;
            %    trackObj = trackObj.addFrame(1, struct('Area',5));
            %    trackObj = trackObj.addFrame(2, struct('Area',10));
            %    trackObj = trackObj.addFrame(3, struct('Area',20));
            %    trackObj = trackObj.addFrame(4, struct('Area',40));
            %
            %    %Delete frame 2
            %    trackObj = trackObj.deleteFrame(2);
            %
            %    %Delete frames 1 and 4
            %    trackOb = trackObj.deleteFrame([1, 4]);
            %
            %  See also: ZNSTrack.updateTrack
            
            %Validate the frame index input
            if isnumeric(frameIndex)
                if ~(all(frameIndex >= obj.FirstFrame & frameIndex <= obj.LastFrame))
                    error('ZNSTrack:deleteFrame:frameIndexInvalid',...
                        'The frame index should be between %d and %d.',...
                        obj.FirstFrame, obj.LastFrame);
                end
                
                %Convert the frame index into the index for the data array
                dataInd = frameIndex - obj.FirstFrame + 1;
                
            elseif islogical(frameIndex)
                if (numel(frameIndex) ~= obj.NumFrames) || (~isvector(frameIndex))
                    error('ZNSTrack:deleteFrame:frameIndexInvalidSize',...
                        'If the frame index is a logical array, it must be a vector with the same number of elements as the number of frames.');
                end
                
                %If it is a logical array, the usual deletion syntax should
                %work
                dataInd = frameIndex;
                
                %Calculate the frame indices to delete
                frameIndex = obj.FirstFrame + find(dataInd) - 1;
                
            elseif ischar(frameIndex)
                
                if any(strcmpi(frameIndex,{'last','end'}))
                    dataInd = numel(obj.Data);
                    frameIndex = obj.LastFrame;
                    
                else
                    error('ZNSTrack:deleteFrame:frameIndexCharInvalid',...
                        'Expected the frame index to be a number, a logical array, or ''last''.');
                end
                
            else
                error('ZNSTrack:deleteFrame:frameIndexNotNumericOrLogical',...
                    'Expected the frame index to be a number or a logical array.');
            end
            
            %Remove the frame(s)
            obj.Data(dataInd) = [];
            
            for iF = 1:numel(frameIndex)
                if frameIndex(iF) == obj.FirstFrame
                    obj.FrameIndex(1) = [];
                    
                elseif frameIndex(iF) == obj.LastFrame
                    obj.FrameIndex(end) = [];
                    
                else
                    obj.FrameIndex = obj.FirstFrame:obj.FirstFrame + numel(obj.Data) - 1;
                end
            end
            
        end
        
        function valueOut = getLatestData(obj, varargin)
            %GETLATESTDATA  Get data from last frame of the specified property
            %
            %  V = T.GETLATESTDATA(dataPropertyName) returns the last frame
            %  of the specified data property.           
            
            if nargin > 1            
                if ~ismember(varargin{1},fieldnames(obj.Data))
                    error('ZNSTrack:getLatestData:InvalidDataPropertyName', ...
                        '%s is not a tracked property name.', varargin{1});
                end
                
                valueOut = obj.Data(end).(varargin{1});
            else
                valueOut = obj.Data(end);
            end
            
        end

        function frameData = getFrame(obj, tFrame)
            %GETFRAME  Get measured data from a frame
            %
            %  S = T.GETFRAME(tFrame) will get the data from the frame
            %  specified, returning it as a structure S.
            
            %Check that the frame exists
            if tFrame < obj.FirstFrame || tFrame > obj.LastFrame
                error('ZNSTrack:getFrameData:InvalidFrame',...
                    'Frame %d does not exist.',tFrame);
            end
            
            frameData = obj.Data(tFrame - obj.FirstFrame + 1);
            
        end        
        
        function propList = getPropertyList(obj)
            %GETPROPERTYLIST  Get list of measured property names
            %
            %  C = T.GETPROPERTYLIST returns a list of fieldnames of the
            %  Data property as a cell string array A.
            
            propList = fieldnames(obj.Data);
            
        end
        
        function describe(obj)
            %DESCRIBE  Prints useful information about the track
            %
            %  A.DESCRIBE will show the mean, min, max, and std of the
            %  measured track properties.
            
            %Get list of track properties
            trackProps = obj.getPropertyList;
            
            dataOut(numel(trackProps)) = struct('Property',[],...
                'Mean', [], 'StDev', [], 'Min', [], 'Max', []);
            
            %Collect data
            for iP = 1:numel(trackProps)
                dataOut(iP).Property = trackProps{iP};
                
                currPropValues = cat(1,obj.Data.(trackProps{iP}));
                
                dataOut(iP).Mean = mean(currPropValues,1,'omitnan');
                dataOut(iP).StDev = std(currPropValues,1,'omitnan');
                dataOut(iP).Min = min(currPropValues,[],1,'omitnan');
                dataOut(iP).Max = max(currPropValues,[],1,'omitnan');
                
            end
            
            tblOut = struct2table(dataOut);
            disp(tblOut);
            
        end
        
        function data = getData(obj, dataPropertyName,varargin)
            %GETDATA  Get specified data 
            %
            %  D = T.GETDATA(P) gets data property P from the track,
            %  returning it as a vector in D.
            %
            %  D = T.GETDATA(P,'last') gets data property P in the last
            %  frame of the track, returning it as a vector in D.
            %
            %  D will be the same length as T.NumFrames, and will contain
            %  NaNs.
            
            if isempty(varargin)
                data = nan(obj.NumFrames, size(obj.Data(1).(dataPropertyName),2));
                for ii = 1:obj.NumFrames
                    if ~isnan(obj.Data(ii).(dataPropertyName))
                        data(ii,:) = obj.Data(ii).(dataPropertyName);
                    end
                end
            else
                
                if strcmpi(varargin{1},'last')
                    data = obj.Data(end).(dataPropertyName);
                end
                  
            end
        end
        
    end
    
end