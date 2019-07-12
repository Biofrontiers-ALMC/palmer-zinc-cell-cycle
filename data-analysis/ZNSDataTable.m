classdef ZNSDataTable
    %ZNSDATATABLE  Calculates data from track array
    %
    %  The data in this table will be used for filtering operations
    
    properties (SetAccess = private) %Changing these properties should trigger a recalculation on the affected columns
        
        RestingFRETinterval = [-8, -3];   %In frames
        MaxFRETinterval = [1, 5];
        MitosisInterval = 5; %Number of frames pre and post mitosis when calculating the diffs
        
        %Pre-mitosis thresholds
        CDK2PreMitosisWindow = [-5, -1];
        
        %Post-mitosis thresholds
        CDK2ratioThreshold = 0.5;
        CDK2PostMitosisWindow = [3, 8];    %In frames
        
        CDK2Window = [20 30];  %Window for arbitrary CDK2 levels
        
        CDK2DiffWindow = [-10 10]; %Window for CDK2 drop
        
    end
    
    properties (SetAccess = private) %Data columns
        
        Data = struct('TrackIdx', {},...
            'MotherIdx',{}, ...
            'DaughterIdx',{}, ...
            'DidDivide', {}, ...
            'NuclDiff', {}, ...
            'NuclMaxRatio', {}, ...
            'RestingFRET', {}, ...
            'RestingFRETCorr', {}, ...
            'FRETDiff', {}, ...
            'FRETDiffCorr', {}, ...
            'FRETmaxOffset', {}, ...
            'CFPInitial', {}, ...
            'FirstFrame', {}, ...
            'LastFrame', {}, ...
            'MitosisFrame', {},...
            'Generation', {}, ...
            'IntermitoticTime', {},...
            'MaxFRETratio', {},...
            'MinFRETratio', {}, ...
            'CDK2Classification', {});
        
    end
    
    methods
        
        function obj = ZNSDataTable(varargin)
            %Constructor
            
            if ~isempty(varargin)
                
                if isa(varargin{1},'ZNSTrackArray')
                    
                    trackData = varargin{1};
                    
                elseif exist(varargin{1},'file')
                    
                    load(varargin{1},'trackData');      
                    
                end
                
                obj = obj.analyzeTrackArray(trackData);
%             else
%                 %Open a GUI to select file if TrackArrayObj is not present
%                 [filename, pathname] = uigetfile({'*.mat','MAT file (*.mat)'},...
%                     'Select file tracks to load');
%                 
%                 if filename == 0
%                     return;                  
%                 end
%                 
%                 data = load(fullfile(pathname,filename),'trackData');
%                 obj = obj.analyzeTrackArray(data.trackData);
            end
        end
        
        function recOut = getRecord(obj, recIdx)
            %Get data record
            
            recOut = obj.Data(recIdx);
            
        end
        
        function numRec = numel(obj)
            
            numRec = numel(obj.Data);
            
        end
        
        function [dataOut, isValid] = get(obj, reqData, varargin)
            %GET  Returns requested data from the data column
            %
            %  D = DT.GET('NuclDiff') returns all the nuclear difference
            %  D = DT.GET('TrackIdx', 'NuclDiff > 300') returns all the
            %  track indices for nuclear differences > 300
            
            %Make sure that the requested data is a field in the Data
            %property
            if ~isfield(obj.Data,reqData)
                error('ZNSDataTable:get:InvalidRequest',...
                    '%s is not a data column.', reqData);
            end
            
            if strcmp(reqData,'TrackIdx') ||strcmp(reqData,'CDK2Classification')
                %TrackIdx is the only column with non-uniform number of
                %cols
                dataOut = {obj.Data.(reqData)};
            else
                dataOut = cat(1,obj.Data.(reqData));
            end
            
            %Parse varargin if it is populated
            if ~isempty(varargin)
                
                %Split string into operations at AND and OR operator
                [expList, opList] = strsplit(varargin{1}, {'&','|','AND','OR'});
                
                %Parse each expression
                for iExp = 1:numel(expList)
                    
                    %Tokens will be {data column} {Operator} {value}
                    strToParse = regexprep(expList{iExp},'\s','');
                    
                    [propName, matches] = strsplit(strToParse,{'>','<','=','~'});
                    
                      tokens = regexp(strToParse,'(\w*)([><=~]{1,2})([0-9]*)','tokens');
                    
                    if isempty(propName)
                        error('Invalid string %s',strToParse);
                    end
                    
                    comparisonStr = [matches{1},propName{2}];                    
                    
                    %Concatenate the data
                    temp = nan(1, numel(obj.Data));
                    for ii = 1:numel(obj.Data)
%                         if numel(obj.Data(ii).(propName{1})) ~= 1
%                             disp(numel(obj.Data(ii).(propName{1})))
%                             keyboard
%                             
%                         end
                        
                        if ~isempty(obj.Data(ii).(propName{1}))
                            temp(ii) = obj.Data(ii).(propName{1});
                        end
                        
                    end
                    
                    filtCol = cat(1, obj.Data.(propName{1}));      
                    
                    if numel(filtCol) ~= numel(obj.Data)
                        keyboard
                    end
                    
                    %filtCol = temp;
                    
                    %Check the relation
                    if ~exist('isValid','var')
                        isValid = eval(sprintf('filtCol %s',comparisonStr));
                    else
                        %Run the comparison depending on the operator
                        switch opList{iExp - 1}
                            
                            case {'&', 'AND', '&&'}
                                try
                                isValid = isValid & ...
                                    eval(sprintf('filtCol %s',comparisonStr));
                                catch
                                    keyboard
                                end
                            case {'|', 'OR', '||'}
                                isValid = isValid | ...
                                    eval(sprintf('filtCol %s',comparisonStr));
                            
                        end
                    end
                end
                                                    
                %Remove rows which are not valid
                dataOut(~isValid) = [];
            else
                isValid = true(numel(dataOut), 1);
            end
            
        end
        
        function isObjEmpty = isempty(obj)
            
            isObjEmpty = isempty(obj.Data);
            
        end
        
        function obj = setRestingFRETinterval(obj, restingInterval, varargin)
            %SETRESTINGFRETINTERVAL  Change the resting FRET interval
            %
            %  DT = DT.SETRESTINGFRETINTERVAL([M, N]) will set the resting
            %  FRET interval property of the object, as long as the object
            %  is empty. M and N must be negative numbers, with M < N,
            %  specifying the window in number of frames before mitosis.
            %
            %  DT = DT.SETRESTINGFRETINTERVAL([M, N], A) If the data table
            %  is not empty, the track array is required to recalculate the
            %  affected data (resting FRET and FRET diff).
            
            %Validate the inputs
            if numel(restingInterval) ~= 2
                error('ZNSDataTable:setRestingFRETinterval:InsufficientIntervalInputs',...
                    'The interval must be specified as a 1x2 vector.')
            else
                if any(restingInterval > 0) || restingInterval(1) > restingInterval(2)
                    error('ZNSDataTable:setRestingFRETinterval:InvalidInterval',...
                        'The interval values must be < 0 and the first value must be less than the second.')
                end
            end
           
            if ~isempty(obj) 
                if isempty(varargin)
                    error('ZNSDataTable:setRestingFRETinterval:TrackArrayNeeded',...
                        'Track Array needed since data table is not empty.')
                elseif ~isa(varargin{1},'ZNSTrackArray')
                    error('ZNSDataTable:setRestingFRETinterval:InputNotTrackArray',...
                        'Expected the input to be a ZNSTrackArray object.')
                end
            end
            
            obj.RestingFRETinterval = restingInterval;
            if ~isempty(obj)
                
                for ii = 1:numel(obj.Data)
                    
                    if obj.Data(ii).DidDivide  %Update the resting FRET and FRET diff values
                        
                        %Get the mother track FRET ratio corrected
                        motherFRETint = varargin{1}.getTrack(obj.Data(ii).MotherIdx).getData('FRETratio');
                        daughterFRETint = varargin{1}.getTrack(obj.Data(ii).DaughterIdx).getData('FRETratio');

                        [obj.Data(ii).RestingFRET, obj.Data(ii).FRETDiff] = obj.measureResting(motherFRETint, daughterFRETint);
                        
                        motherFRETint = varargin{1}.getTrack(obj.Data(ii).MotherIdx).getData('FRETratioCorr');
                        daughterFRETint = varargin{1}.getTrack(obj.Data(ii).DaughterIdx).getData('FRETratioCorr');
                        
                        [obj.Data(ii).RestingFRETCorr, obj.Data(ii).FRETDiffCorr, obj.Data(ii).FRETmaxOffset] = obj.measureResting(motherFRETint, daughterFRETint);
                        
                    end
                    
                end
                
            end
            
        end
        
        function obj = setMaxFRETinterval(obj, maxFRETinterval, varargin)
            %SETRESTINGFRETINTERVAL  Change the resting FRET interval
            %
            %  DT = DT.SETRESTINGFRETINTERVAL([M, N]) will set the resting
            %  FRET interval property of the object, as long as the object
            %  is empty. M and N must be negative numbers, with M < N,
            %  specifying the window in number of frames before mitosis.
            %
            %  DT = DT.SETRESTINGFRETINTERVAL([M, N], A) If the data table
            %  is not empty, the track array is required to recalculate the
            %  affected data (resting FRET and FRET diff).
            
            %Validate the inputs
            if numel(maxFRETinterval) ~= 2
                error('ZNSDataTable:setRestingFRETinterval:InsufficientIntervalInputs',...
                    'The interval must be specified as a 1x2 vector.')
            else
                if maxFRETinterval(2) < maxFRETinterval(1)
                    error('ZNSDataTable:setRestingFRETinterval:InvalidInterval',...
                        'The interval values must be < 0 and the first value must be less than the second.')
                end
            end
            
            if ~isempty(obj)
                if isempty(varargin)
                    error('ZNSDataTable:setRestingFRETinterval:TrackArrayNeeded',...
                        'Track Array needed since data table is not empty.')
                elseif ~isa(varargin{1},'ZNSTrackArray')
                    error('ZNSDataTable:setRestingFRETinterval:InputNotTrackArray',...
                        'Expected the input to be a ZNSTrackArray object.')
                end
            end
            
            obj.MaxFRETinterval = maxFRETinterval;
            if ~isempty(obj)
                
                for ii = 1:numel(obj.Data)
                    
                    if obj.Data(ii).DidDivide  %Update the resting FRET and FRET diff values
                        
                        %Get the mother track FRET ratio corrected
                        motherFRETint = varargin{1}.getTrack(obj.Data(ii).MotherIdx).getData('FRETratio');
                        daughterFRETint = varargin{1}.getTrack(obj.Data(ii).DaughterIdx).getData('FRETratio');
                        
                        [obj.Data(ii).RestingFRET, obj.Data(ii).FRETDiff] = obj.measureResting(motherFRETint, daughterFRETint);
                        
                        %Get the mother track FRET ratio corrected
                        motherFRETint = varargin{1}.getTrack(obj.Data(ii).MotherIdx).getData('FRETratioCorr');
                        daughterFRETint = varargin{1}.getTrack(obj.Data(ii).DaughterIdx).getData('FRETratioCorr');
                        
                        [obj.Data(ii).RestingFRETCorr, obj.Data(ii).FRETDiffCorr, obj.Data(ii).FRETmaxOffset] = obj.measureResting(motherFRETint, daughterFRETint);
                        
                        
                    end
                    
                end
                
            end
            
        end
        
        function obj = setMitosisInterval(obj, mitosisInterval, varargin)
            %SETMITOSISINTERVAL  Change the mitosis interval
            %
            %  The NuclDiff property is computed by taking the difference
            %  between the nuclear intensity at the end of the mother track
            %  with the minimum value N data points on either side of 
            %  mitosis.
            %
            %  If the object is empty, DT = DT.SETMITOSISINTERVAL(N) will
            %  set the number of frames N before and after mitosis used to
            %  calculate the NuclDiff property.
            %
            %  If the data table is not empty, the track array A is
            %  required to recalculate the affected data (NuclDiff). In
            %  this case, use DT = DT.SETRESTINGFRETINTERVAL([M, N], A).
            
            %Validate the inputs
            if numel(mitosisInterval) ~= 1
                error('ZNSDataTable:setMitosisInterval:InvalidNumberOfInputs',...
                    'Interval value should be a single number.')
            else
                if mitosisInterval <= 0
                    error('ZNSDataTable:setMitosisInterval:InvalidInterval',...
                        'Expected interval to be > 0.')
                end
            end
           
            if ~isempty(obj) 
                if isempty(varargin)
                    error('ZNSDataTable:setMitosisInterval:TrackArrayNeeded',...
                        'Track Array needed since data table is not empty.')
                elseif ~isa(varargin{1},'ZNSTrackArray')
                    error('ZNSDataTable:setMitosisInterval:InputNotTrackArray',...
                        'Expected the input to be a ZNSTrackArray object.')
                end
            end
            
            %Update the object            
            obj.MitosisInterval = mitosisInterval;
            
            if ~isempty(obj)
                for ii = 1:numel(obj.Data)
                    
                    if obj.Data(ii).DidDivide  %Update the resting FRET and FRET diff values
                        
                        %Get the mother track FRET ratio corrected
                        motherNuclInt = varargin{1}.getTrack(obj.Data(ii).MotherIdx).getData('NuclInt');
                        daughterNuclInt = varargin{1}.getTrack(obj.Data(ii).DaughterIdx).getData('NuclInt');
                        
                        obj.Data(ii).NuclDiff = obj.measureDiff(motherNuclInt,daughterNuclInt);
                        
                    end
                    
                end
            end
        end
        
        
        function obj = setCDK2PostMitosis(obj, cdk2Interval, cdk2threshold, varargin)
            %SETRESTINGFRETINTERVAL  Change the resting FRET interval
            %
            %  DT = DT.SETRESTINGFRETINTERVAL([M, N]) will set the resting
            %  FRET interval property of the object, as long as the object
            %  is empty. M and N must be negative numbers, with M < N,
            %  specifying the window in number of frames before mitosis.
            %
            %  DT = DT.SETRESTINGFRETINTERVAL([M, N], A) If the data table
            %  is not empty, the track array is required to recalculate the
            %  affected data (resting FRET and FRET diff).
            
            %Validate the inputs
            if numel(cdk2Interval) ~= 2
                error('ZNSDataTable:setCDK2Interval:InsufficientIntervalInputs',...
                    'The interval must be specified as a 1x2 vector.')
            else
                if any(cdk2Interval <= 0) || cdk2Interval(2) < cdk2Interval(1)
                    error('ZNSDataTable:setCDK2Interval:InvalidInterval',...
                        'The interval values must be > 0 and the first value must be less than the second.')
                end
            end
            
            if ~isempty(obj)
                if isempty(varargin)
                    error('ZNSDataTable:setCDK2Interval:TrackArrayNeeded',...
                        'Track Array needed since data table is not empty.')
                elseif ~isa(varargin{1},'ZNSTrackArray')
                    error('ZNSDataTable:setCDK2Interval:InputNotTrackArray',...
                        'Expected the input to be a ZNSTrackArray object.')
                end
            end
            
            obj.CDK2PostMitosisWindow = cdk2Interval;
            obj.CDK2ratioThreshold = cdk2threshold;
            
            if ~isempty(obj)
                for ii = 1:numel(obj.Data)
                    
                    if obj.Data(ii).DidDivide  %Update the resting FRET and FRET diff values
                        
                        cdk2Daughter = varargin{1}.getTrack(obj.Data(ii).DaughterIdx).getData('CDK2ratioCorr');
                        
                        if obj.CDK2PostMitosisWindow(2) > numel(cdk2Daughter)
                            obj.Data(ii).MaxCDK2AfterMitosis = NaN;
                        else
                            obj.Data(ii).MaxCDK2AfterMitosis = max(cdk2Daughter(obj.CDK2PostMitosisWindow));
                        end
                        
                        obj.Data(ii).CDK2Classification = ...
                            obj.classifyCDK2(varargin{1}.getTrack(obj.Data(ii).DaughterIdx).getData('CDK2ratioCorr'));
                                                
                    end
                end
                
            end
            
        end
        
        function obj = setCDK2PreMitosis(obj, cdk2Interval, varargin)
            %SETRESTINGFRETINTERVAL  Change the resting FRET interval
            %
            %  DT = DT.SETRESTINGFRETINTERVAL([M, N]) will set the resting
            %  FRET interval property of the object, as long as the object
            %  is empty. M and N must be negative numbers, with M < N,
            %  specifying the window in number of frames before mitosis.
            %
            %  DT = DT.SETRESTINGFRETINTERVAL([M, N], A) If the data table
            %  is not empty, the track array is required to recalculate the
            %  affected data (resting FRET and FRET diff).
            
            %Validate the inputs
            if numel(cdk2Interval) ~= 2
                error('ZNSDataTable:setCDK2Interval:InsufficientIntervalInputs',...
                    'The interval must be specified as a 1x2 vector.')
            else
                if any(cdk2Interval > 0) || cdk2Interval(2) < cdk2Interval(1)
                    error('ZNSDataTable:setCDK2Interval:InvalidInterval',...
                        'The interval values must be > 0 and the first value must be less than the second.')
                end
            end
            
            if ~isempty(obj)
                if isempty(varargin)
                    error('ZNSDataTable:setCDK2Interval:TrackArrayNeeded',...
                        'Track Array needed since data table is not empty.')
                elseif ~isa(varargin{1},'ZNSTrackArray')
                    error('ZNSDataTable:setCDK2Interval:InputNotTrackArray',...
                        'Expected the input to be a ZNSTrackArray object.')
                end
            end
            
            obj.CDK2PreMitosisWindow = cdk2Interval;
            
            if ~isempty(obj)
                for ii = 1:numel(obj.Data)
                    
                    if obj.Data(ii).DidDivide  %Update the Max CDK2 before mitosis
                        
                        cdk2Mother = varargin{1}.getTrack(obj.Data(ii).MotherIdx).getData('CDK2ratioCorr');
                        
                        if (numel(cdk2Mother) + obj.CDK2PreMitosisWindow(1)) <= 0
                            obj.Data(ii).MaxCDK2BeforeMitosis = NaN;
                            obj.Data(ii).MeanCDK2BeforeMitosis = NaN;
                        else
                            obj.Data(ii).MaxCDK2BeforeMitosis = max(cdk2Mother( (numel(cdk2Mother) + obj.CDK2PreMitosisWindow(1)):(numel(cdk2Mother) +  obj.CDK2PreMitosisWindow(2))));
                            obj.Data(ii).MeanCDK2BeforeMitosis = mean(cdk2Mother( (numel(cdk2Mother) + obj.CDK2PreMitosisWindow(1)):(numel(cdk2Mother) +  obj.CDK2PreMitosisWindow(2))));;
                        end
                        
                    end
                end
                
            end
            
        end
        
        
        function obj = setCDK2Window(obj, cdk2Interval, varargin)
            %SETCDK2WINDOW  Change window to calculate mean CDK2 value
            %
            %  DT = DT.SETRESTINGFRETINTERVAL([M, N]) will set the resting
            %  FRET interval property of the object, as long as the object
            %  is empty. M and N must be negative numbers, with M < N,
            %  specifying the window in number of frames before mitosis.
            %
            %  DT = DT.SETRESTINGFRETINTERVAL([M, N], A) If the data table
            %  is not empty, the track array is required to recalculate the
            %  affected data (resting FRET and FRET diff).
            
            %Validate the inputs
            if numel(cdk2Interval) ~= 2
                error('ZNSDataTable:setCDK2Interval:InsufficientIntervalInputs',...
                    'The interval must be specified as a 1x2 vector.')
            else
                if any(cdk2Interval < 0) || cdk2Interval(2) < cdk2Interval(1)
                    error('ZNSDataTable:setCDK2Interval:InvalidInterval',...
                        'The interval values must be > 0 and the first value must be less than the second.')
                end
            end
            
            if ~isempty(obj)
                if isempty(varargin)
                    error('ZNSDataTable:setCDK2Interval:TrackArrayNeeded',...
                        'Track Array needed since data table is not empty.')
                elseif ~isa(varargin{1},'ZNSTrackArray')
                    error('ZNSDataTable:setCDK2Interval:InputNotTrackArray',...
                        'Expected the input to be a ZNSTrackArray object.')
                end
            end
            
            obj.CDK2Window = cdk2Interval;
            
            if ~isempty(obj)
                for recIdx = 1:numel(obj.Data)
                    
                    if obj.Data(recIdx).DidDivide 
                        
                        cdk2Mother = varargin{1}.getTrack(obj.Data(recIdx).MotherIdx).getData('CDK2ratioCorr');
                        cdk2Daughter = varargin{1}.getTrack(obj.Data(recIdx).DaughterIdx).getData('CDK2ratioCorr');
                        
                        cdk2 = [cdk2Mother; cdk2Daughter];
                        
                        %Compute CDK2 ratio within a certain time window
                        %Calculate the offset from frames to indicies
                        startI = obj.CDK2Window(1) - varargin{1}.getTrack(obj.Data(recIdx).MotherIdx).FirstFrame + 1;
                        endI = startI + (obj.CDK2Window(2) - obj.CDK2Window(1));
                        
                        if startI <= 0 || endI > numel(cdk2)
                            obj.Data(recIdx).MeanCDK2inWindow = NaN;
                            obj.Data(recIdx).MaxCDK2inWindow = NaN;
                            obj.Data(recIdx).MinCDK2inWindow = NaN;
                        else
                            obj.Data(recIdx).MeanCDK2inWindow = mean(cdk2(startI:endI));
                            obj.Data(recIdx).MaxCDK2inWindow = max(cdk2(startI:endI));
                            obj.Data(recIdx).MinCDK2inWindow = min(cdk2(startI:endI));
                        end
                    else
                        cdk2 = varargin{1}.getTrack(obj.Data(recIdx).TrackIdx).getData('CDK2ratioCorr');               
                   
                        %Compute CDK2 ratio within a certain time window
                        %Calculate the offset from frames to indicies
                        startI = obj.CDK2Window(1) - varargin{1}.getTrack(obj.Data(recIdx).TrackIdx).FirstFrame + 1;
                        endI = startI + (obj.CDK2Window(2) - obj.CDK2Window(1));
                        
                        if startI <= 0 || endI > numel(cdk2)
                            obj.Data(recIdx).MeanCDK2inWindow = NaN;
                            obj.Data(recIdx).MaxCDK2inWindow = NaN;
                            obj.Data(recIdx).MinCDK2inWindow = NaN;
                        else
                            obj.Data(recIdx).MeanCDK2inWindow = mean(cdk2(startI:endI));
                            obj.Data(recIdx).MaxCDK2inWindow = max(cdk2(startI:endI));
                            obj.Data(recIdx).MinCDK2inWindow = min(cdk2(startI:endI));
                        end

                        
                    end
                    
                end
                
            end
            
        end
                
        
        function obj = analyzeTrackArray(obj, trackArrayObj)
            %ANALYZETRACKARRAY  Analyze a track array object
            %
            %  DT = DT.ANALYZETRACKARRAY(A) will analyze the ZNSTrackArray
            %  object A to build up the data tables required.
            
            %Validate the input
            if ~isa(trackArrayObj, 'ZNSTrackArray')
                error('ZNSDataTable:analyzeTrackArray:InputNotTrackArray',...
                    'Expected the input to be a ZNSTrackArray object.');
            end
            
            %Initialize the size of the Data property. The MAXIMUM size
            %this structure can be is twice the number of dividing tracks +
            %the number of non-dividing tracks.
            numRecords = 2 * trackArrayObj.NumDividingMotherTracks + ... 
                trackArrayObj.NumNonDividingTracks;                      
            
            obj.Data(numRecords).DidDivide = false;
            
            %----- Calculate the data values -----
            ctrData = 0;
            for iTr = 1:numel(trackArrayObj)
                
                currTrack = trackArrayObj.getTrack(iTr);

                if ~all(isnan(currTrack.DaughterIdxs))
                    %Dividing track: Create records for the mother-daughter
                    %pairs
                    ctrData = ctrData + 1;
                    obj = obj.calcTrackInfo(ctrData, trackArrayObj, ...
                        [iTr, currTrack.DaughterIdxs(1)]);
                    
                    if numel(currTrack.DaughterIdxs) == 2
                        %Second entry only created if there is a second
                        %daughter track (if mitosis fixed, there might not
                        %be)
                        ctrData = ctrData + 1;
                        obj = obj.calcTrackInfo(ctrData, trackArrayObj, ...
                            [iTr, currTrack.DaughterIdxs(2)]);
                    end
                    
                elseif isnan(currTrack.MotherIdx)
                    
                    %Non-dividing track (must have both mother & daughter
                    %track = nans)
                    ctrData = ctrData + 1;
                    obj = obj.calcTrackInfo(ctrData, trackArrayObj, ...
                        iTr);
                end
                
            end
            
            %Delete unused data records
            if ctrData < numRecords
                obj.Data((ctrData+1):end) = [];
            end
            
        end
        
        function dataOut = getByTrackIdx(obj, trackIdx, prop)
            %GETBYTRACKIDX  Get data by track index
            %
            %  D = DT.GETBYTRACKIDX(I, Property)
            allTrackIdxs = {obj.Data.TrackIdx};

            A = cellfun(@(x) all(x == trackIdx), allTrackIdxs);
            dataOut = obj.Data(A).(prop);        
        end
        
        
    end
    
    methods (Access = private)
        
        function obj = calcTrackInfo(obj, recIdx, trackArrayObj, trackIdx)
            %CALCTRACKINFO  Populate the data table by calculating track
            %info
            %
            % Syntax: DT.CALCTRACKINFO(newRecIdx, trackArrayObj, I)
            %
            %  I is either a single track Idx for a non-dividing track or
            %  1x2 vector for dividing track
            
            if numel(trackIdx) > 2
                error('Track index cannot be more than 2')                
            end
            
            %Get track data
            track = cell(1,numel(trackIdx));
            for iT = 1:numel(trackIdx)
                track{iT} = trackArrayObj.getTrack(trackIdx(iT));
            end
            
            %Populate common data columns
            obj.Data(recIdx).TrackIdx = trackIdx;
            obj.Data(recIdx).FirstFrame = track{1}.FirstFrame;
            obj.Data(recIdx).LastFrame = track{end}.LastFrame;
            obj.Data(recIdx).TrackLen = track{end}.NumFrames;
            
            obj.Data(recIdx).Rejected = false;  %Rejected track
            
            if trackArrayObj.hasFRETdata
                
                tempCFP = track{1}.getData('CFPIntCorr');
                obj.Data(recIdx).CFPInitial = tempCFP(1);
                
                %Combine the FRETratioCorr
                tempFRETratioCorr = [track{1}.getData('FRETratioCorr');track{end}.getData('FRETratioCorr')];
                
                obj.Data(recIdx).MaxFRETratio = max(tempFRETratioCorr);
                obj.Data(recIdx).MinFRETratio = min(tempFRETratioCorr);
               
            end
            
            obj.Data(recIdx).NuclMaxRatio= max(track{1}.getData('NuclInt'),[],'omitnan') ./ mean(track{1}.getData('NuclInt'),'omitnan');
            
            %Calculate the data values which change
            switch numel(trackIdx)
                
                case 1 %Non-dividing cell
                    obj.Data(recIdx).MotherIdx = nan;
                    obj.Data(recIdx).DaughterIdx = nan;
                    obj.Data(recIdx).DidDivide = false;
                    obj.Data(recIdx).MitosisFrame = nan;
                    
                    obj.Data(recIdx).NuclDiff = nan;
                    
                    obj.Data(recIdx).Generation = nan;
                    obj.Data(recIdx).IntermitoticTime = NaN;
                    
                    if trackArrayObj.hasFRETdata
                        obj.Data(recIdx).RestingFRET = mean(track{1}.getData('FRETratio'));
                        obj.Data(recIdx).FRETDiff = nan;
                        
                        obj.Data(recIdx).RestingFRETCorr = mean(track{1}.getData('FRETratioCorr'));
                        obj.Data(recIdx).FRETDiffCorr = nan;
                        
                         obj.Data(recIdx).FRETmaxOffset = 0;
                    end
                    
                    if trackArrayObj.hasCDK2data
                        obj.Data(recIdx).CDK2Classification = 0;
                        
                        cdk2 = getData(track{1}, 'CDK2ratioCorr');
                        
                        obj.Data(recIdx).MaxCDK2 = max(cdk2);
                        obj.Data(recIdx).MinCDK2 = min(cdk2);
                        
                        obj.Data(recIdx).MeanCDK2BeforeMitosis = NaN;
                        obj.Data(recIdx).MaxCDK2BeforeMitosis = NaN;
                        obj.Data(recIdx).MaxCDK2AfterMitosis = NaN;
                        
                        obj.Data(recIdx).CDK2Diff = NaN;
                        
                        %Compute CDK2 ratio within a certain time window
                        %Calculate the offset from frames to indicies
                        startI = obj.CDK2Window(1) - track{1}.FirstFrame + 1;
                        endI = startI + (obj.CDK2Window(2) - obj.CDK2Window(1));
                        
                        if startI <= 0 || endI > numel(cdk2)
                            obj.Data(recIdx).MeanCDK2inWindow = NaN;
                            obj.Data(recIdx).MaxCDK2inWindow = NaN;
                            obj.Data(recIdx).MinCDK2inWindow = NaN;
                        else
                            obj.Data(recIdx).MeanCDK2inWindow = mean(cdk2(startI:endI));
                            obj.Data(recIdx).MaxCDK2inWindow = max(cdk2(startI:endI));
                            obj.Data(recIdx).MinCDK2inWindow = min(cdk2(startI:endI));
                        end

                        
                    end
                    
                case 2 %Mother-daughter pair
                    
                    obj.Data(recIdx).MotherIdx = trackIdx(1);
                    obj.Data(recIdx).DaughterIdx = trackIdx(2);
                    obj.Data(recIdx).DidDivide = true;
                    obj.Data(recIdx).MitosisFrame = track{1}.LastFrame;
                    obj.Data(recIdx).Generation = obj.getGeneration(trackIdx(1), recIdx);

                    try
                    obj.Data(recIdx).NuclDiff = obj.measureDiff(track{1}.getData('NuclInt'), track{2}.getData('NuclInt'));
                    catch
                        keyboard
                    end
                    
                    
                    if trackArrayObj.hasFRETdata
                        [obj.Data(recIdx).RestingFRET, obj.Data(recIdx).FRETDiff]...
                            = obj.measureResting(track{1}.getData('FRETratio'), track{2}.getData('FRETratio'));
                        
                        [obj.Data(recIdx).RestingFRETCorr, obj.Data(recIdx).FRETDiffCorr, obj.Data(recIdx).FRETmaxOffset]...
                            = obj.measureResting(track{1}.getData('FRETratioCorr'), track{2}.getData('FRETratioCorr'));
                    end
                                        
                    if obj.Data(recIdx).Generation > 1
                        obj.Data(recIdx).IntermitoticTime = track{1}.NumFrames;
                    else
                        obj.Data(recIdx).IntermitoticTime = NaN;
                    end
                    
                    if trackArrayObj.hasCDK2data
                        
                        cdk2Mother = track{1}.getData('CDK2ratioCorr');
                        cdk2Daughter = track{2}.getData('CDK2ratioCorr');
                        
                        if isempty(cdk2Mother) || isempty(cdk2Daughter)
                            keyboard
                            
                        end
                        
                        
                        cdk2 = [cdk2Mother; cdk2Daughter];
                        
                        obj.Data(recIdx).MaxCDK2 = max(cdk2);
                        obj.Data(recIdx).MinCDK2 = min(cdk2);
                        
                        if (numel(cdk2Mother) + obj.CDK2PreMitosisWindow(1)) <= 0
                            obj.Data(recIdx).MaxCDK2BeforeMitosis = NaN;
                            obj.Data(recIdx).MeanCDK2BeforeMitosis = NaN;
                        else
                            obj.Data(recIdx).MaxCDK2BeforeMitosis = max(cdk2Mother( (numel(cdk2Mother) + obj.CDK2PreMitosisWindow(1)):(numel(cdk2Mother) +  obj.CDK2PreMitosisWindow(2))));
                            obj.Data(recIdx).MeanCDK2BeforeMitosis = mean(cdk2Mother( (numel(cdk2Mother) + obj.CDK2PreMitosisWindow(1)):(numel(cdk2Mother) +  obj.CDK2PreMitosisWindow(2))));
                        end
                        
                        if obj.CDK2PostMitosisWindow(2) > numel(cdk2Daughter)
                            obj.Data(recIdx).MaxCDK2AfterMitosis = NaN;
                        else
                            obj.Data(recIdx).MaxCDK2AfterMitosis = max(cdk2Daughter(obj.CDK2PostMitosisWindow));
                        end
                        
                        obj.Data(recIdx).CDK2Classification = classifyCDK2(obj, cdk2Daughter);
                        
                        %Compute CDK2 ratio within a certain time window
                        %Calculate the offset from frames to indicies
                        startI = obj.CDK2Window(1) - track{1}.FirstFrame + 1;
                        endI = startI + (obj.CDK2Window(2) - obj.CDK2Window(1));
                        
                        %Compute the drop in CDK2 ratio within the
                        %specified window
                        if numel(cdk2Mother) > abs(obj.CDK2DiffWindow(1))
                            
                            tempCDK2Mother = cdk2Mother((end + obj.CDK2DiffWindow(1)):end);
      
                        else
                            tempCDK2Mother = cdk2Mother;                            
                        end
                        
                        if numel(cdk2Daughter) > obj.CDK2DiffWindow(2)
                            tempCDK2Daughter = cdk2Daughter(1:obj.CDK2DiffWindow(2));
                        else
                            tempCDK2Daughter = cdk2Daughter;
                        end
                        
                        obj.Data(recIdx).CDK2Diff = max(tempCDK2Mother) - min(tempCDK2Daughter);
                        
                        if startI <= 0 || endI > numel(cdk2)
                            obj.Data(recIdx).MeanCDK2inWindow = NaN;
                            obj.Data(recIdx).MaxCDK2inWindow = NaN;
                            obj.Data(recIdx).MinCDK2inWindow = NaN;
                        else
                            obj.Data(recIdx).MeanCDK2inWindow = mean(cdk2(startI:endI));
                            obj.Data(recIdx).MaxCDK2inWindow = max(cdk2(startI:endI));
                            obj.Data(recIdx).MinCDK2inWindow = min(cdk2(startI:endI));
                        end
                        
                    end
                    
            end
            
        end
        
        function valOut = measureDiff(obj, motherInt, daughterInt)
            %MEASUREDIFF  Measures the difference in the signal at mitosis
            %
            %  S = D.MEASUREDIFF(motherInt, daughterInt)
            
            %If the tracks are too short, then set the value to NaN
            if (length(motherInt) - 1) < obj.MitosisInterval || length(daughterInt) < obj.MitosisInterval
                valOut = NaN;
                return;
            end
            
            %Measure the difference in intensity both pre and post mitosis
            preMitosisDiff = motherInt(end) - motherInt(end - obj.MitosisInterval:end - 1);
            postMitosisDiff = motherInt(end) - daughterInt(1:obj.MitosisInterval);

            %Calculate the maximum difference for each of the two sides
            diffs = [max(preMitosisDiff); max(postMitosisDiff)];
            
            %Throw away values less than zero
            diffs(diffs <= 0) = [];
            
            %The final difference value is the minimum value between
            %the two sides. If the diffs were all negative, then return
            %NaN.
            if ~isempty(diffs)
                valOut = min(diffs);
            else
                valOut = NaN;
            end
            
        end
        
        function [restingVal, FRETdiff, FRETMaxOffset] = measureResting(obj, motherInt, daughterInt)
            %MEASURERESTING  Measure resting FRET value
            
            %If track is too short, return NaNs
            if length(motherInt) - abs(obj.RestingFRETinterval(1)) <= 0 || length(daughterInt) < obj.MaxFRETinterval(2)
                restingVal = NaN;
                FRETdiff = NaN;
                FRETMaxOffset = NaN;
                return;
            end
            try

            %Compute the resting FRET            
            restingVal = mean(motherInt((end + obj.RestingFRETinterval(1)):(end + obj.RestingFRETinterval(2))));
                        
            %Find the max FRET ratio
            [maxFRET, FRETMaxOffset] = max(daughterInt(obj.MaxFRETinterval(1):obj.MaxFRETinterval(2)));            
            
            %Measure the difference in intensity both pre and post mitosis
            FRETdiff = maxFRET - restingVal;            
                        catch
                keyboard
            end
        end
               
        function genNum = getGeneration(obj, motherTrackIdx, currRec)
            %GETGENERATION  Identify the generation of the track
            %
            %  G = M.GETGENERATION(I) will identify the generation of track
            %  I by searching through the list of daughter track indices
            %  for a matching track index. It is expected that MitosisEvent
            %  array will be built in increasing order of track IDs.
            
            %Check if the track exists in the daughter indices list
            %To speed things up, we only check preceeding records.
            
            for iRec = 1:(currRec - 1)
                
                if ismember(motherTrackIdx,obj.getRecord(iRec).DaughterIdx)
                    genNum = obj.Data(iRec).Generation + 1;
                    return;
                end
                
            end
            
            genNum = 1;
        end
       
        function cdk2Class = classifyCDK2(obj, cdk2Ratio)
            
            %cdkRatio = smooth(cdk2Ratio, 5);
            
            if numel(cdk2Ratio) < obj.CDK2PostMitosisWindow(2)
                cdk2Class = 0;
            elseif any(cdk2Ratio(obj.CDK2PostMitosisWindow(1):obj.CDK2PostMitosisWindow(2)) > obj.CDK2ratioThreshold)
                cdk2Class = 1;                
            elseif any(cdk2Ratio(obj.CDK2PostMitosisWindow(1):end) > obj.CDK2ratioThreshold)
                cdk2Class = 2;
            else
                cdk2Class = 3;
            end
            
        end
        
    end
        
end