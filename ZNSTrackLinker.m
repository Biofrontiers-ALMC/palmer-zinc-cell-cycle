classdef ZNSTrackLinker
    %ZNSTRACKLINKER  Links objects by position to form tracks
        
    properties
        
        options = struct('MaxLinkDistance',100,...
            'MaxTrackAge',2,...
            'TrackMitosis',true,...
            'MinAgeSinceMitosis',2,...
            'MaxMitosisDistance',30,...
            'MaxMitosisAreaChange',0.3,...
            'MaxMitosisIntensityDiff', 0.1, ...
            'LAPSolver','lapjv');
        
    end
    
    properties (Hidden) %Active tracks
        
        TrackArray = ZNSTrackArray; %Use this to hold the celltrack object
        
        trackingInfo = struct('Index',{},...
            'Position', {},...
            'Area',{},...
            'Age',{},...
            'AgeSinceDivision',{});

    end
    
    methods
        
        function obj = ZNSTrackLinker(varargin)
            %TRACKLINKER  Constuctor function for the TrackLinker class
            %
            %  O = TRACKLINKER will create a TrackLinker object O.
            
            if nargin > 0
                obj = obj.initializeTracks(varargin{:});
            end
            
        end
        
        function numTracks = numel(obj)
            
            numTracks = numel(obj.trackingInfo);
            
        end
        
        function disp(obj)
            
            fprintf('Cell tracker object\n')
            fprintf('Number of tracked objects: %d\n',numel(obj))
            fprintf('Total tracks: %d\n',numel(obj.TrackArray))
            
        end
        
        function obj = setOptions(obj, varargin)
            %SETOPTIONS  Set options for the linker
            %
            %  linkerObj = linkerObj.SETOPTIONS(parameter, value) will set
            %  the parameter to value.
            %
            %  linkerObj = linkerObj.SETOPTIONS(O) where O is a data object
            %  with the same property names as the options will work.
            %
            %  linkerObj = linkerObj.SETOPTIONS(S) where S is a struct
            %  with the same fieldnames as the options will also work.
            %
            %  Non-matching parameter names will be ignored.
                                   
            if numel(varargin) == 1 && isstruct(varargin{1})
                %Parse a struct as input
                
                inputParameters = fieldnames(varargin{1});
                
                for iParam = 1:numel(inputParameters)
                    if ismember(inputParameters{iParam},fieldnames(obj.options))
                        obj.options.(inputParameters{iParam}) = ...
                            varargin{1}.(inputParameters{iParam});
                    else
                        %Just skip unmatched options
                    end
                    
                end
                
            elseif numel(varargin) == 1 && isobject(varargin{1})
                %Parse an object as input
                
                inputParameters = properties(varargin{1});
                
                for iParam = 1:numel(inputParameters)
                    if ismember(inputParameters{iParam},fieldnames(obj.options))
                        obj.options.(inputParameters{iParam}) = ...
                            varargin{1}.(inputParameters{iParam});
                    else
                        %Just skip unmatched options
                    end
                    
                end
                
            else
                if rem(numel(varargin),2) ~= 0
                    error('Input must be Property/Value pairs.');
                end
                inputArgs = reshape(varargin,2,[]);
                for iArg = 1:size(inputArgs,2)
                    if ismember(inputArgs{1,iArg},fieldnames(obj.options))
                        
                        obj.options.(inputArgs{1,iArg}) = inputArgs{2,iArg};
                    else
                        %Just skip unmatched options
                    end
                end
            end
        end
        
    end
    
    methods
        
        function [obj, newTrackIdx] = addTrackedObject(obj, tFrame, trackParam, didDivide)
            %ADDTRACKEDOBJECT  Add a new track and start tracking
            %
            %  [T, I] = T.ADDTRACKEDOBJECT(tFrame, trackData, divFlag) adds
            %  a new track to the TrackArray object. The index of the new
            %  track is returned in I.
            %
            %  divFlag indicates whether the cell divided recently. If so,
            %  the track will not be considered for a second mitosis until
            %  a certain amount of time has passed. divFlag should be true
            %  if the new track is created due to a mitosis event, and
            %  false otherwise.
            
            %Create the new track
            [obj.TrackArray, newTrackIdx] = obj.TrackArray.addTrack(tFrame,trackParam);
            
            %Add it to the list of active tracks
            idxNewObject = numel(obj.trackingInfo) + 1;
            
            obj.trackingInfo(idxNewObject).Index = newTrackIdx;

            %Get location of Position and Area, then use these as the
            %values
            obj.trackingInfo(idxNewObject).Position = trackParam.Centroid;
            obj.trackingInfo(idxNewObject).Area = trackParam.Area;
            obj.trackingInfo(idxNewObject).Age = 0;
            
            %Check if cell division flag was enabled
            if ~exist('didDivide','var')
                didDivide = false;
            end
            
            if didDivide
                obj.trackingInfo(idxNewObject).AgeSinceDivision = 0;
            else
                obj.trackingInfo(idxNewObject).AgeSinceDivision = Inf;
            end
            
        end
        
        function obj = ageTrackedObject(obj,idx)
            
            obj.trackingInfo(idx).age = obj.trackingInfo(idx).age + 1;
            
        end
        
        function obj = removeOldTracks(obj)
            %Removes tracks which have not been updated in awhile. The
            %number of frames allowed is set in the options 'MaxTrackAge'.
            trackAges = [obj.trackingInfo.Age];
            
            obj.trackingInfo(trackAges > obj.options.MaxTrackAge) = [];
            
        end
        
        function costMatrix = makeLinkMatrix(obj, newPositions)
            %Calculates the cost matrix for the linear assignment operation
            
            currPositions = cat(1,obj.trackingInfo.Position);
            
            %Calculate the cost to link based on object distances
            costToLink = obj.calcDistance(currPositions,newPositions);
            costToLink(costToLink > obj.options.MaxLinkDistance) = Inf;
            
            maxCostToLink = max(costToLink(costToLink < Inf));
            
            %Costs for stopping
            stopCost = diag(1.05 * maxCostToLink * ones(1,numel(obj)));
            stopCost(stopCost == 0) = Inf;
            
            %Cost to start a new segment
            segStartCost = diag(1.05 * maxCostToLink * ones(1,size(newPositions,1)));
            segStartCost(segStartCost == 0) = Inf;
            
            %Auxiliary matrix
            auxMatrix = costToLink';
            auxMatrix(auxMatrix < Inf) = min(costToLink(costToLink < Inf));
            
            costMatrix = [costToLink, stopCost; segStartCost, auxMatrix];
            
        end
                
        function obj = assignToTrack(obj,tFrame,varargin)
            %ASSIGNTOTRACK  Assign data to tracks
            %
            %  L = L.ASSIGNTOTRACK(frameNumber, dataIn) will assign the
            %  detected objects in the frame specified either to new tracks
            %  or to existing tracks.
            %
            %  dataIn should be parameter/value lists, with each new
            %  detection in a new row.
           
            if isstruct(varargin{1})
                frameData = varargin{1};
            else
                frameData = ZNSTrackLinker.paramlistToStruct(varargin{:});
            end            
            
            %If the object is empty, initialize new tracks
            if numel(obj) == 0
                obj = obj.initializeTracks(tFrame, frameData);       
                return;
            end
                      
            %Calculate the costMatrix
            costMat = obj.makeLinkMatrix(cat(1,frameData.Centroid));
            
            %Solve the assignment problem
            switch obj.options.LAPSolver
                case 'munkres'
                    assignments = obj.munkres(costMat);
                    
                case 'lapjv'
                    assignments = obj.lapjv(costMat);
            end
            
            %-----Process the assignments-----%
            
            nExistingTracks = numel(obj.trackingInfo);
            nNewDetections = numel(frameData);
            
            %First set of numbers are associations with the current tracks
            for iM = 1:nExistingTracks
                
                assignedValue = assignments(iM);
                
                if assignedValue > 0 && assignedValue <= nNewDetections
                    
                    %Update the track
                    obj = obj.updateTrackedObject(iM,tFrame,frameData(assignedValue));
                    
                else
                    
                    %Age the tracks which were not updated
                    obj.trackingInfo(iM).Age = obj.trackingInfo(iM).Age + 1;
                    
                end
            end
            
            %Stop track segments that have not been updated in awhile
            obj = obj.removeOldTracks;
            
            %Second set of assignments are 'start segments'
            for iN = 1:nNewDetections
                
                assignedValue = assignments(nExistingTracks + iN);
                
                if assignedValue > 0 && assignedValue <= nNewDetections
                    
                    %Test for cell division
                    if obj.options.TrackMitosis
                        [isMitosis, motherIdx] = obj.testForMitosis(frameData(assignedValue));
                        
                        if isMitosis
                            
                            %Get the mother track index
                            motherTrackIdx = obj.trackingInfo(motherIdx).Index;
                            
                            %If the "mother cell" was created at the same frame, then it
                            %is not a mitosis event
                            if obj.TrackArray.getTrack(motherTrackIdx).FirstFrame == tFrame
                                isMitosis = false;
                            end
                        end
                        
                        if isMitosis

                            %Make two new tracks, with assigned mother cells
                            [obj, daughterIdx1] = obj.addTrackedObject(tFrame,frameData(assignedValue),1);
                            
                            currMotherTrack = obj.TrackArray.getTrack(motherTrackIdx);
                            [obj, daughterIdx2] = obj.addTrackedObject(tFrame,currMotherTrack.getLatestData,1);
                            
                            %Update mother track indices
                            obj.TrackArray = obj.TrackArray.updateTrackProperty(daughterIdx1,'MotherIdx',motherTrackIdx);
                            obj.TrackArray = obj.TrackArray.updateTrackProperty(daughterIdx2,'MotherIdx',motherTrackIdx);
                            
                            %Remove the last frame from the mother track
                            obj.TrackArray = obj.TrackArray.deleteFrame(motherTrackIdx, tFrame);
                            obj.TrackArray = obj.TrackArray.updateTrackProperty(motherTrackIdx,'DaughterIdxs',[daughterIdx1,daughterIdx2]);
                            
                            %Stop tracking mother cell
                            obj.trackingInfo(motherIdx) = [];
                            
                        else
                            %If not mitosis event, create a new track
                            obj = obj.addTrackedObject(tFrame,frameData(assignedValue),0);
                        end
                        
                    else
                        %Create a new track
                        obj = obj.addTrackedObject(tFrame,frameData(assignedValue),0);
                    end
                    
                else
                    
                    %Do nothing, it's assigned to the auxiliary matrix
                    
                end
                
                
            end
            
        end
        
        function obj = updateTrackedObject(obj,idx,tFrame,frameDataIn)
            
            %Update the celltrack object
            trackToUpdate = obj.trackingInfo(idx).Index;
            obj.TrackArray = obj.TrackArray.updateTrack(trackToUpdate,tFrame,frameDataIn);
            
            %Update tracking info
            obj.trackingInfo(idx).Position = frameDataIn.Centroid;
            obj.trackingInfo(idx).Area = frameDataIn.Area;
            obj.trackingInfo(idx).Age = 0;
            obj.trackingInfo(idx).AgeSinceDivision = obj.trackingInfo(idx).AgeSinceDivision + 1;
            
        end
                
        function trackOut = getTrack(obj,trackIdx)
            %GETTRACK  Get a track from the track array
            
            trackOut = obj.TrackArray.getTrack(trackIdx);
            
        end
        
        function [isMitosis, motherIdx] = testForMitosis(obj,newData)
            %TESTFORMITOSIS  Test for mitosis events            
            
            isMitosis = false;  %Default value
            motherIdx = NaN;
            
            currPos = newData.Centroid;
            currArea = newData.Area;
            
            validTrackedPos = cat(1,obj.trackingInfo.Position);
            
            %Exclude cells which were not updated in the current
            excludedCells = 1:numel(obj);
            excludedCells = excludedCells([obj.trackingInfo.Age] > 0);
            
            for iEx = excludedCells
                validTrackedPos(iEx, :) = [Inf Inf];
            end
            
            %Test to see if there is a particle nearby
            if numel(validTrackedPos) == 0
                return;
            end
            
            %Calculate the distance of the objects, blocking invalid
            %entries
            distances = obj.calcDistance(validTrackedPos,currPos);
            distances(distances == 0) = Inf;
            distances(distances > obj.options.MaxMitosisDistance) = Inf;
            
            if all(isinf(distances))
                return;
            end
            
            %Given the group of nearby neighbors, look for one that fits
            %the size criteria
            if ~all(isnan(distances))
                for nn = (find(~isnan(distances) & ~isinf(distances)))'
                    
                    %Is the nearby neighbour particle of similar size?
                    areaNN = obj.trackingInfo(nn).Area;
                    
                    if abs(areaNN - currArea)/currArea < obj.options.MaxMitosisAreaChange
                        
                        %If nuclear intensity is present, check the
                        %intensity difference
                        if isfield(newData,'NuclInt')
                            
                            intensityNN = obj.getTrack(obj.trackingInfo(nn).Index).Data(end).NuclInt;
                            
                            currIntensity = newData.NuclInt;
                            
                            if (abs(intensityNN - currIntensity)/ currIntensity) > obj.options.MaxMitosisIntensityDiff
                                continue;                                
                            end
                            
                        end
                        
                        %If nearest neighbour divided recently, then don't allow it to
                        %be classified as a division
                        if obj.trackingInfo(nn).AgeSinceDivision < obj.options.MinAgeSinceMitosis
                            %Is there another particle that fits nearby?
                            continue;
                        end
                        
                        isMitosis = true;
                        motherIdx = nn;
                        
                        return
                    end
                end
            end
    
            
        end
                
        function trackArrayOut = getTrackArray(obj)
            %GETTRACKARRAY  Gets the track array
            %
            %  A = L.GETTRACKARRAY will return the TrackArray object A
            %  (i.e. all track data).
            
            trackArrayOut = obj.TrackArray;
            
        end
        
    end
    
    methods (Access = private)
        
        function obj = initializeTracks(obj,tFrame,varargin)
            %INITIALIZETRACKS  Initialize tracks when the linker is empty
            %
            %  L = L.INITIALIZETRACKS(tFrame, inputData) will initialize
            %  the TrackLinker object L with tracks from the inputData.
            %
            %  inputData can be a structure or Parameter/Value pairs. If
            %  Parameter/Value pairs are used, the data for each new object
            %  should be in a new row.
            
            if nargin < 3
                error('ZNSTrackLinker:initializeTracks:InsufficientInputs',...
                    'Input arguments must include the frame time and the track data.');
            end
            
            if ~isstruct(varargin{1})
                
                %TODO: Convert Parameter/Value pairs into a struct
                error('Invalid input. Must be struct.');
                
            else
                trackData = varargin{1};
            end
           
            for iTrack = 1:numel(trackData)
                obj = obj.addTrackedObject(tFrame,trackData(iTrack),0);
            end
            
        end
        
    end
        
    methods (Access = private, Hidden = true, Static)
        
        function distancesOut = calcDistance(originalPos, newPos)
            %CALCDISTANCE  Calculates the Euclidean distance between two sets of points
            %
            %Returns N x M matrix where N is the number of original
            %positions and M is the number of new positions
            
            distancesOut = zeros(size(originalPos,1),size(newPos,1));
            
            for iR = 1:size(originalPos,1)
                
                currPos = originalPos(iR,:);
                
                orX = currPos(1);
                orY = currPos(2);
                
                newPosX = newPos(:,1);
                newPosY = newPos(:,2);
                
                distancesOut(iR,:) = sqrt((orX - newPosX).^2 + (orY - newPosY).^2);
                
            end
            
        end
        
        function [optimAssign, optimCost, unassigned_cols] = munkres(costMatrix)
            %MUNKRES  Munkres (Hungarian) linear assignment
            %
            %  [I, C] = MUNKRES(M) returns the column indices I assigned to each row,
            %  and the minimum cost C based on the assignment. The cost of the
            %  assignments are given in matrix M, with workers along the rows and tasks
            %  along the columns. The matrix optimizes the assignment by minimizing the
            %  total cost.
            %
            %  The code can deal with partial assignments, i.e. where M is not a square
            %  matrix. Unassigned rows (workers) will be given a value of 0 in the
            %  output I. [I, C, U] = MUNKRES(M) will give the index of unassigned
            %  columns (tasks) in vector U.
            %
            %  The algorithm attempts to speed up the process in the case where values
            %  of a row or column are all Inf (i.e. impossible link). In that case, the
            %  row or column is excluded from the assignment process; these will be
            %  automatically unassigned in the result.
            
            %This code is based on the algorithm described at:
            %http://csclab.murraystate.edu/bob.pilgrim/445/munkres.html
            
            %Get the size of the matrix
            [nORows, nOCols] = size(costMatrix);
            
            %Check for rows and cols which are all infinity, then remove them
            validRows = ~all(costMatrix == Inf,2);
            validCols = ~all(costMatrix == Inf,1);
            
            nRows = sum(validRows);
            nCols = sum(validCols);
            
            nn = max(nRows,nCols);
            
            if nn == 0
                error('Invalid cost matrix: Cannot be all Inf.')
            elseif any(isnan(costMatrix(:))) || any(costMatrix(:) < 0)
                error('Invalid cost matrix: Expected costs to be all positive numbers.')
            end
            
            %Make a new matrix
            tempCostMatrix = ones(nn) .* (10 * max(max(costMatrix(costMatrix ~= Inf))));
            tempCostMatrix(1:nRows,1:nCols) = costMatrix(validRows,validCols);
            
            tempCostMatrix(tempCostMatrix == Inf) = realmax;
            
            %Get the minimum values of each row
            rowMin = min(tempCostMatrix,[],2);
            
            %Subtract the elements in each row with the corresponding minima
            redMat = bsxfun(@minus,tempCostMatrix,rowMin);
            
            %Mask matrix (0 = not a zero, 1 = starred, 2 = primed)
            mask = zeros(nn);
            
            %Vectors of column and row numbers
            rowNum = 1:nn;
            colNum = rowNum;
            
            %Row and column covers (1 = covered, 0 = uncovered)
            rowCover = zeros(1,nn);
            colCover = rowCover;
            
            %Search for unique zeros (i.e. only one starred zero should exist in each
            %row and column
            for iRow = rowNum(any(redMat,2) == 0)
                for iCol = colNum(any(redMat(iRow,:) == 0))
                    if (redMat(iRow,iCol) == 0 && rowCover(iRow) == 0 && colCover(iCol) == 0)
                        mask(iRow,iCol) = 1;
                        rowCover(iRow) = 1;
                        colCover(iCol) = 1;
                    end
                end
            end
            
            %Clear the row cover
            rowCover(:) = 0;
            
            %The termination condition is when each column has a single starred zero
            while ~all(colCover)
                
                %---Step 4: Prime an uncovered zero---%
                %Find a non-covered zero and prime it.
                %If there is no starred zero in the row containing this primed zero,
                %proceed to step 5.
                %Otherwise, cover this row and uncover the column contianing the
                %starred zero.
                %Continue until there are no uncovered zeros left. Then get the minimum
                %value and proceed to step 6.
                
                stop = false;
                
                %Find an uncovered zero
                for iRow = rowNum( (any(redMat == 0,2))' & (rowCover == 0) )
                    for iCol = colNum(redMat(iRow,:) == 0)
                        
                        if (redMat(iRow,iCol) == 0) && (rowCover(iRow) == 0) && (colCover(iCol) == 0)
                            mask(iRow,iCol) = 2;    %Prime the zero
                            
                            if any(mask(iRow,:) == 1)
                                rowCover(iRow) = 1;
                                colCover(mask(iRow,:) == 1) = 0;
                            else
                                
                                %Step 5: Augment path algorithm
                                currCol = iCol; %Initial search column
                                storePath = [iRow, iCol];
                                
                                %Test if there is a starred zero in the current column
                                while any(mask(:,currCol) == 1)
                                    %Get the (row) index of the starred zero
                                    currRow = find(mask(:,currCol) == 1);
                                    
                                    storePath = [storePath; currRow, currCol];
                                    
                                    %Find the primed zero in this row (there will
                                    %always be one)
                                    currCol = find(mask(currRow,:) == 2);
                                    
                                    storePath = [storePath; currRow, currCol];
                                end
                                
                                %Unstar each starred zero, star each primed zero in the
                                %searched path
                                indMask = sub2ind([nn,nn],storePath(:,1),storePath(:,2));
                                mask(indMask) = mask(indMask) - 1;
                                
                                %Erase all primes
                                mask(mask == 2) = 0;
                                
                                %Uncover all rows
                                rowCover(:) = 0;
                                
                                %Step 3: Cover the columns with stars
                                colCover(:) = any((mask == 1),1);
                                
                                stop = true;
                                break;
                            end
                        end
                        
                        %---Step 6---
                        
                        %Find the minimum uncovered value
                        minUncVal = min(min(redMat(rowCover == 0,colCover== 0)));
                        
                        %Add the value to every element of each covered row
                        redMat(rowCover == 1,:) = redMat(rowCover == 1,:) + minUncVal;
                        
                        %Subtract it from every element of each uncovered column
                        redMat(:,colCover == 0) = redMat(:,colCover == 0) - minUncVal;
                    end
                    
                    if (stop)
                        break;
                    end
                end
                
            end
            
            %Assign the outputs
            optimAssign = zeros(nORows,1);
            optimCost = 0;
            
            unassigned_cols = 1:nCols;
            
            validRowNum = 1:nORows;
            validRowNum(~validRows) = [];
            
            validColNum = 1:nOCols;
            validColNum(~validCols) = [];
            
            %Only assign valid workers
            for iRow = 1:numel(validRowNum)
                
                assigned_col = colNum(mask(iRow,:) == 1);
                
                %Only assign valid tasks
                if assigned_col > numel(validColNum)
                    %Assign the output
                    optimAssign(validRowNum(iRow)) = 0;
                else
                    optimAssign(validRowNum(iRow)) = validColNum(assigned_col);
                    
                    %         %Calculate the optimized (minimized) cost
                    optimCost = optimCost + costMatrix(validRowNum(iRow),validColNum(assigned_col));
                    
                    unassigned_cols(unassigned_cols == assigned_col) = [];
                end
            end
        end
        
        function [rowsol,cost,v,u,costMat] = lapjv(costMat,resolution)
            % LAPJV  Jonker-Volgenant Algorithm for Linear Assignment Problem.
            %
            % [ROWSOL,COST,v,u,rMat] = LAPJV(COSTMAT, resolution) returns the optimal column indices,
            % ROWSOL, assigned to row in solution, and the minimum COST based on the
            % assignment problem represented by the COSTMAT, where the (i,j)th element
            % represents the cost to assign the jth job to the ith worker.
            % The second optional input can be used to define data resolution to
            % accelerate speed.
            % Other output arguments are:
            % v: dual variables, column reduction numbers.
            % u: dual variables, row reduction numbers.
            % rMat: the reduced cost matrix.
            %
            % For a rectangular (nonsquare) costMat, rowsol is the index vector of the
            % larger dimension assigned to the smaller dimension.
            %
            % [ROWSOL,COST,v,u,rMat] = LAPJV(COSTMAT,resolution) accepts the second
            % input argument as the minimum resolution to differentiate costs between
            % assignments. The default is eps.
            %
            % Known problems: The original algorithm was developed for integer costs.
            % When it is used for real (floating point) costs, sometime the algorithm
            % will take an extreamly long time. In this case, using a reasonable large
            % resolution as the second arguments can significantly increase the
            % solution speed.
            %
            % See also munkres, Hungarian
            
            % version 1.0 by Yi Cao at Cranfield University on 3rd March 2010
            % version 1.1 by Yi Cao at Cranfield University on 19th July 2010
            % version 1.2 by Yi Cao at Cranfield University on 22nd July 2010
            % version 2.0 by Yi Cao at Cranfield University on 28th July 2010
            % version 2.1 by Yi Cao at Cranfield University on 13th August 2010
            % version 2.2 by Yi Cao at Cranfield University on 17th August 2010
            % version 3.0 by Yi Cao at Cranfield University on 10th April 2013
            
            % This Matlab version is developed based on the orginal C++ version coded
            % by Roy Jonker @ MagicLogic Optimization Inc on 4 September 1996.
            % Reference:
            % R. Jonker and A. Volgenant, "A shortest augmenting path algorithm for
            % dense and spare linear assignment problems", Computing, Vol. 38, pp.
            % 325-340, 1987.
            
            %
            % Examples
            % Example 1: a 5 x 5 example
            %{
                    [rowsol,cost] = lapjv(magic(5));
                    disp(rowsol); % 3 2 1 5 4
                    disp(cost);   %15
            %}
            % Example 2: 1000 x 1000 random data
            %{
                    n=1000;
                    A=randn(n)./rand(n);
                    tic
                    [a,b]=lapjv(A);
                    toc                 % about 0.5 seconds
            %}
            % Example 3: nonsquare test
            %{
                    n=100;
                    A=1./randn(n);
                    tic
                    [a,b]=lapjv(A);
                    toc % about 0.2 sec
                    A1=[A zeros(n,1)+max(max(A))];
                    tic
                    [a1,b1]=lapjv(A1);
                    toc % about 0.01 sec. The nonsquare one can be done faster!
                    %check results
                    disp(norm(a-a1))
                    disp(b-b)
            %}
            
            if nargin<2
                maxcost=min(1e16,max(max(costMat)));
                resolution=eps(maxcost);
            end
            % Prepare working data
            [rdim,cdim] = size(costMat);
            M=min(min(costMat));
            if rdim>cdim
                costMat = costMat';
                [rdim,cdim] = size(costMat);
                swapf=true;
            else
                swapf=false;
            end
            dim=cdim;
            costMat = [costMat;2*M+zeros(cdim-rdim,cdim)];
            costMat(costMat~=costMat)=Inf;
            maxcost=max(costMat(costMat<Inf))*dim+1;
            if isempty(maxcost)
                maxcost = Inf;
            end
            costMat(costMat==Inf)=maxcost;
            % free = zeros(dim,1);      % list of unssigned rows
            % colist = 1:dim;         % list of columns to be scaed in various ways
            % d = zeros(1,dim);       % 'cost-distance' in augmenting path calculation.
            % pred = zeros(dim,1);    % row-predecessor of column in augumenting/alternating path.
            v = zeros(1,dim);         % dual variables, column reduction numbers.
            rowsol = zeros(1,dim)-1;  % column assigned to row in solution
            colsol = zeros(dim,1)-1;  % row assigned to column in solution
            
            numfree=0;
            free = zeros(dim,1);      % list of unssigned rows
            matches = zeros(dim,1);   % counts how many times a row could be assigned.
            % The Initilization Phase
            % column reduction
            for j=dim:-1:1 % reverse order gives better results
                % find minimum cost over rows
                [v(j), imin] = min(costMat(:,j));
                if ~matches(imin)
                    % init assignement if minimum row assigned for first time
                    rowsol(imin)=j;
                    colsol(j)=imin;
                elseif v(j)<v(rowsol(imin))
                    j1=rowsol(imin);
                    rowsol(imin)=j;
                    colsol(j)=imin;
                    colsol(j1)=-1;
                else
                    colsol(j)=-1; % row already assigned, column not assigned.
                end
                matches(imin)=matches(imin)+1;
            end
            
            % Reduction transfer from unassigned to assigned rows
            for i=1:dim
                if ~matches(i)      % fill list of unaasigned 'free' rows.
                    numfree=numfree+1;
                    free(numfree)=i;
                else
                    if matches(i) == 1 % transfer reduction from rows that are assigned once.
                        j1 = rowsol(i);
                        x = costMat(i,:)-v;
                        x(j1) = maxcost;
                        v(j1) = v(j1) - min(x);
                    end
                end
            end
            
            % Augmenting reduction of unassigned rows
            loopcnt = 0;
            while loopcnt < 2
                loopcnt = loopcnt + 1;
                % scan all free rows
                % in some cases, a free row may be replaced with another one to be scaed next
                k = 0;
                prvnumfree = numfree;
                numfree = 0;    % start list of rows still free after augmenting row reduction.
                while k < prvnumfree
                    k = k+1;
                    i = free(k);
                    % find minimum and second minimum reduced cost over columns
                    x = costMat(i,:) - v;
                    [umin, j1] = min(x);
                    x(j1) = maxcost;
                    [usubmin, j2] = min(x);
                    i0 = colsol(j1);
                    if usubmin - umin > resolution
                        % change the reduction of the minmum column to increase the
                        % minimum reduced cost in the row to the subminimum.
                        v(j1) = v(j1) - (usubmin - umin);
                    else % minimum and subminimum equal.
                        if i0 > 0 % minimum column j1 is assigned.
                            % swap columns j1 and j2, as j2 may be unassigned.
                            j1 = j2;
                            i0 = colsol(j2);
                        end
                    end
                    % reassign i to j1, possibly de-assigning an i0.
                    rowsol(i) = j1;
                    colsol(j1) = i;
                    if i0 > 0 % ,inimum column j1 assigned easier
                        if usubmin - umin > resolution
                            % put in current k, and go back to that k.
                            % continue augmenting path i - j1 with i0.
                            free(k)=i0;
                            k=k-1;
                        else
                            % no further augmenting reduction possible
                            % store i0 in list of free rows for next phase.
                            numfree = numfree + 1;
                            free(numfree) = i0;
                        end
                    end
                end
            end
            
            % Augmentation Phase
            % augment solution for each free rows
            for f=1:numfree
                freerow = free(f); % start row of augmenting path
                % Dijkstra shortest path algorithm.
                % runs until unassigned column added to shortest path tree.
                d = costMat(freerow,:) - v;
                pred = freerow(1,ones(1,dim));
                collist = 1:dim;
                low = 1; % columns in 1...low-1 are ready, now none.
                up = 1; % columns in low...up-1 are to be scaed for current minimum, now none.
                % columns in up+1...dim are to be considered later to find new minimum,
                % at this stage the list simply contains all columns.
                unassignedfound = false;
                while ~unassignedfound
                    if up == low    % no more columns to be scaned for current minimum.
                        last = low-1;
                        % scan columns for up...dim to find all indices for which new minimum occurs.
                        % store these indices between low+1...up (increasing up).
                        minh = d(collist(up));
                        up = up + 1;
                        for k=up:dim
                            j = collist(k);
                            h = d(j);
                            if h<=minh
                                if h<minh
                                    up = low;
                                    minh = h;
                                end
                                % new index with same minimum, put on index up, and extend list.
                                collist(k) = collist(up);
                                collist(up) = j;
                                up = up +1;
                            end
                        end
                        % check if any of the minimum columns happens to be unassigned.
                        % if so, we have an augmenting path right away.
                        for k=low:up-1
                            if colsol(collist(k)) < 0
                                endofpath = collist(k);
                                unassignedfound = true;
                                break
                            end
                        end
                    end
                    if ~unassignedfound
                        % update 'distances' between freerow and all unscanned columns,
                        % via next scanned column.
                        j1 = collist(low);
                        low=low+1;
                        i = colsol(j1); %line 215
                        x = costMat(i,:)-v;
                        h = x(j1) - minh;
                        xh = x-h;
                        k=up:dim;
                        j=collist(k);
                        vf0 = xh<d;
                        vf = vf0(j);
                        vj = j(vf);
                        vk = k(vf);
                        pred(vj)=i;
                        v2 = xh(vj);
                        d(vj)=v2;
                        vf = v2 == minh; % new column found at same minimum value
                        j2 = vj(vf);
                        k2 = vk(vf);
                        cf = colsol(j2)<0;
                        if any(cf) % unassigned, shortest augmenting path is complete.
                            i2 = find(cf,1);
                            endofpath = j2(i2);
                            unassignedfound = true;
                        else
                            i2 = numel(cf)+1;
                        end
                        % add to list to be scaned right away
                        for k=1:i2-1
                            collist(k2(k)) = collist(up);
                            collist(up) = j2(k);
                            up = up + 1;
                        end
                    end
                end
                % update column prices
                j1=collist(1:last+1);
                v(j1) = v(j1) + d(j1) - minh;
                % reset row and column assignments along the alternating path
                while 1
                    i=pred(endofpath);
                    colsol(endofpath)=i;
                    j1=endofpath;
                    endofpath=rowsol(i);
                    rowsol(i)=j1;
                    if (i==freerow)
                        break
                    end
                end
            end
            rowsol = rowsol(1:rdim);
            u=diag(costMat(:,rowsol))-v(rowsol)';
            u=u(1:rdim);
            v=v(1:cdim);
            cost = sum(u)+sum(v(rowsol));
            costMat=costMat(1:rdim,1:cdim);
            costMat = costMat - u(:,ones(1,cdim)) - v(ones(rdim,1),:);
            if swapf
                costMat = costMat';
                t=u';
                u=v';
                v=t;
            end
            if cost>maxcost
                cost=Inf;
            end
            
        end
        
        function structOut = paramlistToStruct(varargin)
            %PARAMLISTTOSTRUCT  Convert a Parameter/Value list to struct
            
            if rem(numel(varargin),2) ~= 0
                error('ZNSTrackLinker:paramlistToStruct:InvalidNumberOfInputs',...
                    'Inputs must be matched Parameter/Value pairs.');                
            end
            
            %Reshape the input
            params = reshape(varargin, 2, []);
            
            %Get number of entries (= number of rows in value)
            numParams = size(params{2,1},2);
            
            for iParam = 1:size(params,1)
                structOut(numParams).(params{iParam,1}) = []; %#ok<AGROW>
                
                for iEntries = 1:numParams
                    structOut(iEntries).(params{iParam,1}) = params{iParam,2};
                end
            end
            
        end
        
    end
    
end