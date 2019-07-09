classdef FRETtracker < handle
    %FRETtracker  Tracks cells
    %
    %  This is the main class for the cell tracking module. Use this class
    %  to process a video to segment and track cells.
    %
    %  FRETtracker Properties:
    %  
    %  FRETtracker Methods:
    %      processVideo - Segment and track cells from a dataset
    %
    %  Example:
    %
    %  %Intialize a new FRETtracker object
    %  T = FRETtracker;
    %
    %  %Process a specified file
    %  T.processVideo(filename, outputDir, (optional) Parameter/Value)
    %
    %    
    %  See also: FRETtracker.processVideo
    %
    %  Copyright 2017 University of Colorado Boulder

    properties
        
        %Image properties
        nuclearChan = 'Cy5';
        
        measureFRET = true;
        FRETlocation = 'cyto'; %Or nucl
        FRETChan = 'CFP-YFP Emi_s';
        CFPChan = 'CFP_s';
        
        measureCDK2 = true;
        CDK2Chan = 'mCherry';
        
        CellType = 'DHBmCherry';
        
        %Output options
        OutputMovie = true;
        ExportMasks = false;
        
        %Processing options
        FrameRange = Inf;
        ROI = 'full';
        ParallelProcess = false;
        
        %Maximum and minimum nuclear area
        MinNucleusArea = 15;
        MaxNucleusArea = 5000;
        
        %Segmentation options
        AutoAdjustNuclImg = false;
        
        ThresholdLvl = NaN;
        CytoRingSize = 2;
        
        CorrectLocalBackground = true;
        NumBackgroundBlocks = [11 11];
        BgPrctile = 5;
        
        %Track linker options
        MaxLinkDistance = 200;
        MaxTrackAge = 2;
        TrackMitosis = true;
        MinAgeSinceMitosis = 2;
        MaxMitosisDistance = 30;
        MaxMitosisAreaChange = 0.3;
        MaxMitosisIntensityDiff = 0.1;
        LAPSolver = 'lapjv';
        
    end
    
    properties (Hidden)
        
        EnableDebugMode = false; %Set to true to enable debug mode
        
    end
    
    methods %Public methods
        
        function obj = FRETtracker
            %FRETtracker  Constructor function
            %
            %  T = FRETTRACKER will create a FRETtracker object T. This
            %  function will check that the BioformatsImage toolbox is
            %  installed, prompting the user to download it otherwise.
            
            if ~exist('BioformatsImage','file')
                error('FRETtracker:BioformatsImageNotInstalled',...
                    'A required toolbox ''BioformatsImage'' is not installed. Visit https://biof-git.colorado.edu/core-code/bioformats-image-toolbox/wikis/home to get the latest version.');
            end
            
        end
           
        function processFiles(obj, varargin)
            %PROCESSFILES  Process the selected file(s)
            %
            %  T.PROCESSFILES(filelist) will process all the files
            %  specified.
            %
            %  T.PROCESSFILES will provide a pop-up box which will
            %  allow you to select multiple files.
            
            ip = inputParser;
            ip.addOptional('Filelist', '',@(x) ischar(x) || iscellstr(x));
            ip.addOptional('OutputDir', '',@(x) ischar(x));
            ip.KeepUnmatched = true;
            ip.parse(varargin{:});
            
           
            %----- Get input files ----%
            %Get a directory if none was supplied
            if isempty(ip.Results.Filelist)

                [filelist,dataDir] = uigetfile(...
                    {'*.nd2; *.tif', 'Image files (*.nd2, *.tif)'},...
                    'Choose files to process',...
                    'MultiSelect', 'on');
                
                if ~iscellstr(filelist)
                    if filelist == 0
                        return;
                    end
                end
            else
                if ischar(ip.Results.Filelist)
                    filelist = {ip.Results.Filelist};
                else
                    filelist = ip.Results.Filelist;
                end
                
                %Get the data directory
                dataDir = fileparts(filelist{1});
            end
            
            obj.checkForSettingsFile(dataDir);
            
            %----- Check for output directory -----%
            if isempty(ip.Results.OutputDir)
                outputDir = uigetdir(dataDir,'Select output directory');
                
                if outputDir == 0
                    %Cancelled
                   return; 
                end
                
            else
                outputDir = ip.Results.OutputDir;
            end
            
            %Make the output directory if it doesn't exist
            if ~exist(outputDir,'dir')
                mkdir(outputDir)
            end
            
            %If masks are exported, create a sub-directory to store them
            if obj.ExportMasks
                if ~exist(fullfile(outputDir,'Masks'),'dir')
                    mkdir(fullfile(outputDir,'Masks'));
                end
            end
            
            %----- Current settings -----%
            if ~isempty(ip.Unmatched)
                obj.setOptions(ip.Unmatched);
            end
            
            obj.checkSettings;
            
            %----- Process files -----%
            
            %Enable parallel processing (if selected and more than 1 file is running)
            if obj.ParallelProcess && numel(filelist) > 1
                nParWorkers = Inf;
                
                %Starts a parallel pool object one is not already running.
                %Otherwise it returns the current parallel pool.
                currParpool = gcp;
                
                %Attach the BioformatsImage object to the parallel pool
                %(otherwise it complains about not being able to locate the
                %JAR files).
                currParpool.addAttachedFiles({'BioformatsImage','TIFFreader','ZNSTrackLinker','ZNSTrackArray', 'ZNSTrack'});
                
            else
                nParWorkers = 0;
            end
                        
            optionsStruct = obj.exportOptionsToStruct;
            
            if ~iscell(filelist) && ischar(filelist)
               filelist =  {filelist};
            end
            
            diary(fullfile(outputDir,'log.txt'));      %Start a log file
            
            parfor (iF = 1:numel(filelist), nParWorkers)
                try
                    fprintf('%s %s: Starting processing.\n', datestr(now), filelist{iF});
                    %Call the processing function
                    FRETtracker.analyzeFile(fullfile(dataDir,filelist{iF}),outputDir, optionsStruct);
                    fprintf('%s %s: Completed.\n', datestr(now), filelist{iF});
                catch ME
                    fprintf('%s %s: An error occured:\n', datestr(now), filelist{iF});
                    fprintf('%s \n',getReport(ME,'extended','hyperlinks','off'));
                    
                    
                    
                end
            end
            
            diary off
            
            %Save the processing settings to the save directory
            obj.exportSettings(fullfile(outputDir,'trackerSettings.txt'));
            
            %Ask user whether to open save directory
            openDir = input('Processing completed. Change current path to save directory (Y = Yes)? ','s');
            if strcmpi(openDir,'y')
                cd(outputDir);
            end
            
        end
        
        function setOptions(obj, varargin)
            %SETOPTIONS  Set options file
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
                    if ismember(inputParameters{iParam},properties(obj))
                        obj.(inputParameters{iParam}) = ...
                            varargin{1}.(inputParameters{iParam});
                    else
                        %Just skip unmatched options
                    end
                    
                end
                
            elseif numel(varargin) == 1 && isobject(varargin{1})
                %Parse an object as input
                
                inputParameters = properties(varargin{1});
                
                for iParam = 1:numel(inputParameters)
                    if ismember(inputParameters{iParam},properties(obj))
                        obj.(inputParameters{iParam}) = ...
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
                    if ismember(inputArgs{1,iArg},properties(obj))
                        
                        obj.(inputArgs{1,iArg}) = inputArgs{2,iArg};
                    else
                        %Just skip unmatched options
                    end
                end
            end
        end
                
        function importSettings(obj, filename)
            %IMPORTSETTINGS  Import settings from file
            %
            %  S = FRETtracker.IMPORTSETTINGS(filename) will import
            %  settings from the file specified. The file should be a txt
            %  file.
            
            fid = fopen(filename,'r');
            
            if fid == -1
                error('FRETtracker:importSettings:ErrorReadingFile',...
                    'Could not open file %s for reading.',filename);
            end
            
            while ~feof(fid)
                currLine = strtrim(fgetl(fid));
                
                if isempty(currLine)
                    %Empty lines should be skipped
                    
                elseif strcmpi(currLine(1),'%') || strcmpi(currLine(1),'#')
                    %Lines starting with '%' or '#' are comments, so ignore
                    %those
                    
                else
                    %Expect the input to be PARAM_NAME = VALUE
                    parsedLine = strsplit(currLine,'=');
                    
                    %Get parameter name (removing spaces)
                    parameterName = strtrim(parsedLine{1});
                    
                    %Get value name (removing spaces)
                    value = strtrim(parsedLine{2});
                    
                    if isempty(value)
                        %If value is empty, just use the default
                    else
                        obj.setOptions(parameterName,eval(value));
                    end
                    
                    %Handle the hidden parameter EnableDebugMode
                    if strcmpi(parameterName,'EnableDebugMode')
                        obj.EnableDebugMode = eval(value);
                    end
                    
                end
                
            end
            
            fclose(fid);
            
        end
        
        function exportSettings(obj, exportFilename)
            %EXPORTSETTINGS  Export settings to a txt file
            %
            %  FRETtrackerOptions.EXPORTSETTINGS(filename) will export the
            %  settings to a txt file. If the filename is not provided, a
            %  dialog box will pop-up asking the user to select a location
            %  to save the file.
            
            if ~exist('exportFilename','var')
                
                [filename, pathname] = uiputfile({'*.txt','Text file (*.txt)'},...
                    'Select output file location');
                
                if filename == 0
                    %If cancel button was pressed, then do nothing
                    return;                    
                end
                
                exportFilename = fullfile(pathname,filename);
                
            end
            
            fid = fopen(exportFilename,'w');
            
            if fid == -1
                error('FRETtrackerOptions:exportSettings:CouldNotOpenFile',...
                    'Could not open file to write')
            end
            
            %Write a header
            fprintf(fid,'%%FRETtracker Settings\r\n');
            
            propertyList = properties(obj);
                        
            %Write output data depending on the datatype of the value
            for ii = 1:numel(propertyList)
                
                if ischar(obj.(propertyList{ii}))
                    fprintf(fid,'%s = ''%s'' \r\n',propertyList{ii}, ...
                        obj.(propertyList{ii}));
                    
                elseif isnumeric(obj.(propertyList{ii}))
                    fprintf(fid,'%s = %s \r\n',propertyList{ii}, ...
                        mat2str(obj.(propertyList{ii})));
                    
                elseif islogical(obj.(propertyList{ii}))
                    
                    if obj.(propertyList{ii})
                        fprintf(fid,'%s = true \r\n',propertyList{ii});
                    else
                        fprintf(fid,'%s = false \r\n',propertyList{ii});
                    end
                    
                end
                
            end
            
            fclose(fid);
            
        end
        
    end
    
    methods (Access = private)
        
        function structOut = exportOptionsToStruct(obj)
            %EXPORTOPTIONSTOSTRUCT  Export options to struct
            %
            %  Exports the object properties as a struct
            
            propList = properties(obj);
            for iP = 1:numel(propList)
                structOut.(propList{iP}) = obj.(propList{iP});
            end
            
            %Add the hidden property EnableDebugMode
            structOut.EnableDebugMode = obj.EnableDebugMode;
        end
        
        function checkSettings(obj)
            %CHECKSETTINGS  Confirm current settins with user
            %
            %  O.CHECKSETTINGS will display the current object properties
            %  and ask the user to confirm they are correct. If the user
            %  indicates they are not, the program will prompt for a
            %  settings file.
            
            %Display current settings and ask user to confirm whether to
            %continue
            promptForSettings = true;
            while promptForSettings
                
                disp(obj);
                
                startTrack = input('Verify that the settings are correct. Type ''(S)tart'' to start processing, ''(L)oad'' to load a settings file, or anything else to cancel: ','s');
                
                if (strcmpi(startTrack,'l') || strcmpi(startTrack,'load'))
                    
                        [filename, pathname] = uigetfile({'*.*','All files'},...
                            'Select the settings file');
                        
                        obj.importSettings(fullfile(pathname,filename));
                   
                elseif (strcmpi(startTrack,'y') || strcmpi(startTrack,'yes') ||...
                        strcmpi(startTrack,'s') || strcmpi(startTrack,'start'))
                    %Proceed with the rest of the code
                    promptForSettings = false;
                    
                else
                    error('FRETtracker:processDir:UserCancelled',...
                        'User cancelled operation.');
                end
            end
            
        end
        
        function checkForSettingsFile(obj, dataDir)
            %CHECKFORSETTINGSFILE  Check for settings file
            %
            %  CHECKFORSETTINGSFILE(dataDir) will search for a settings
            %  file in the data directory specified.
            
            filelist = dir(fullfile(dataDir,'*.txt'));
            
            if ~isempty(filelist)
                
                for iFile = 1:numel(filelist)
                    fid = fopen(fullfile(dataDir,filelist(iFile).name),'r');
                    
                    if fid == -1
                        %On read error, just continue
                        continue;
                    end
                    
                    %Check if the file is a settings file by looking for a
                    %header
                    firstLine = fgetl(fid);
                    
                    if ~isempty(regexp(firstLine,'(FRETtracker Settings)*','once'))
                        isSetting = true;
                        settingFile = filelist(iFile).name;
                        break;
                    end
                
                    fclose(fid);
                end
                
                if isSetting
                    %Prompt user whether to load settings
                    loadSetting = input(sprintf('Found settings file ''%s''. Do you want to load it (Y = Yes)? ',settingFile),'s');
                    
                    if strcmpi('y',loadSetting)
                        obj.importSettings(fullfile(dataDir,settingFile));                        
                    end
                    
                end
                
            end
            
        end
        
    end
        
    methods (Static, Access = private)  %Segmentation and tracking functions
        
        function analyzeFile(filename, outputDir, optionsStruct)
            %ANALYZEFILE  Segment and track cells in the specified file
            %
            %  FRETtracker.ANALYZEFILE(filename, outputDir, optionsStruct)
            %  will run the segmentation and tracking algorithm on the
            %  specified file. Output files will be saved to the outputDir
            %  specified. 
            %
            %  Segmentation and tracking options have to be supplied using
            %  the optionsStruct structure.
            %
            %  Note: This method has been made Static to reduce memory
            %  usage. See
            %  https://stackoverflow.com/questions/45041783/referring-to-class-method-in-parfor-loop-significant-memory-usage
            %  Test case 2.
            
            %Get an appropriate reader object for the images
            imgReader = FRETtracker.getImageReader(filename);
            
            %Convert the filename to a well format
            [~,saveFName] = fileparts(imgReader.filename);
            
            %Initialize the VideoWriter object if the option is set to save
            %a movie
            if optionsStruct.OutputMovie
                
%                 %Video Writer for the segmentation movie
%                 vidObjMask = VideoWriter(fullfile(outputDir,[fname,'_masks.avi']));
%                 vidObjMask.Quality = 100;
%                 vidObjMask.FrameRate = 10;
%                 open(vidObjMask);
                
                %Video Writer for the track movie
                vidObjTracks = VideoWriter(fullfile(outputDir,[saveFName,'_tracks.avi']));
                vidObjTracks.Quality = 100;
                vidObjTracks.FrameRate = 10;
                open(vidObjTracks);
                
            end
            
            %--- Begin processing code ---%
            %Determine the range of frames to process
            if isinf(optionsStruct.FrameRange)
                optionsStruct.FrameRange = 1:imgReader.sizeT;
            end
            
            numSkippedFrames = 0;   %Number of frames where no cells were found
            for iT = optionsStruct.FrameRange
                
                %Load the nuclear image
                nuclImg = imgReader.getPlane(1, optionsStruct.nuclearChan, iT);
                
                %Segment and label the cell nuclei
                nuclLabels = FRETtracker.getNuclLabels(nuclImg, optionsStruct);
                 
                %Get frame data
                frameData = FRETtracker.getFrameData(iT, imgReader, nuclLabels, nuclImg, optionsStruct);
                                
                %Link the data to tracks                
                if iT == optionsStruct.FrameRange(1)
                    %Initialize the tracker
                    trackerObj = ZNSTrackLinker;
                    trackerObj = trackerObj.setOptions(optionsStruct);
                    trackerObj = trackerObj.assignToTrack(iT,frameData);
                else
                    try
                        %Assign data to tracks
                        trackerObj = trackerObj.assignToTrack(iT,frameData);
                        
                        %Reset number of skipped frames
                        numSkippedFrames = 0;
                    catch ME
                        %If an error occurs during track assignment, try to
                        %continue if the issue is that no cell was found.
                        %If the error occurs three times in a row, then
                        %stop processing.
                        if all(nuclLabels(:) == 0)
                            numSkippedFrames = numSkippedFrames + 1;
                            
                            if numSkippedFrames > 3
                                error('FRETtracker:analyzeFile:NoCellsFound',...
                                    '%s(%d): No cells were found for three consecutive frames.',...
                                    saveFName,iT)                               
                            end
                        
                            warning('FRETtracker:analyzeFile:NoCellsFound',...
                                '%s(%d): No cells were found.',saveFName,iT)
                        else
                            if optionsStruct.EnableDebugMode && ~optionsStruct.ParallelProcess
                                %Debug mode (only when non-parallel otherwise
                                %execution will freeze)
                                keyboard;
                            end
                            
                            error('FRETtracker:analyzeFile:ErrorDuringAssignment',...
                                '%s(%d): Track assignment error. Original error message:\n%s\n',...
                                saveFName, iT, ME.message);
                        end
                    end
                end
                
                %Draw the movie frame (if set)
                if optionsStruct.OutputMovie
                    %Draw the current frame if saving a movie
                    fh = figure;
                    set(fh,'Visible','off')
                    FRETtracker.plotTracks(nuclImg,trackerObj,nuclLabels,iT);
                    set(fh,'units','normalized','outerposition',[0 0 1 1],'innerposition',[0 0 1 1],'Visible','off')
                    vidObjTracks.writeVideo(getframe(fh));
                    close(fh);
                   
%                     %Generate the cell mask movie. This movie has the
%                     %highest resolution but does not contain tracks or cell
%                     %numbers. Dividing cells are colored red in the movie
%                     %for five frames after division (will be division + 1).
%                     imgOut = FRETtracker.makeMovieFrame(nuclImg,trackerObj,nuclLabels,iT);
%                     vidObjMask.writeVideo(imgOut);
                end
                
                %If option is set, export the cell nuclei labels as a TIF
                %file in the sub-directory 'Masks'.
                if optionsStruct.ExportMasks
                    imwrite(nuclLabels,...
                        fullfile(outputDir,'Masks',sprintf('%s_%d.tif',saveFName,iT)));
                end
                
            end

            if optionsStruct.OutputMovie
                %Close the VideoWriter object
%                 close(vidObjMask)
                close(vidObjTracks)
            end
                        
            %--- Save data ----%
            trackData = trackerObj.getTrackArray;
            
            %Append the timestamp information
            [tsData, tsUnit] = imgReader.getTimestamps(1,optionsStruct.nuclearChan,'TimeRange',optionsStruct.FrameRange);
            trackData = trackData.setTimestamp(tsData, tsUnit); %#ok<NASGU>
            
            metadata = optionsStruct; %#ok<NASGU>
            save(fullfile(outputDir,[saveFName,'.mat']),'trackData','metadata');
            
        end
               
        function plotTracks(nuclImg,trackerObj,nuclLabels,iT)
            %PLOTTRACKS  Shows an image with cell tracks overlaid
            %
            %  FRETtracker.PLOTTRACKS(nuclImg, trackerObj, nuclLabels, iT)

            warning off
            FRETtracker.showoverlay(FRETtracker.normalizeimg(imadjust(nuclImg)),bwperim(nuclLabels > 0),[0 1 0]);
            warning on
            
            hold on
            for iTrack = 1:numel(trackerObj)
                
                trackIdx = trackerObj.trackingInfo(iTrack).Index;
                
                currTrack = trackerObj.getTrack(trackIdx);
                
                if currTrack.FirstFrame <= iT && currTrack.LastFrame >= iT
                    %If the cell exists in the current frame, label it by
                    %drawing the outline and inserting the text
                    
                    labelled = false;
                    if ~isnan(currTrack.MotherIdx)
                        currMother = trackerObj.getTrack(currTrack.MotherIdx);
                        
                        if (iT - currMother.LastFrame) <= 5
                            text(currTrack.Data(end).Centroid(1),currTrack.Data(end).Centroid(2),int2str(trackIdx),'Color',[240,113,113]./255);
                            labelled = true;
                            
                        end
                    end
                    
                    if ~labelled
                        text(currTrack.Data(end).Centroid(1),currTrack.Data(end).Centroid(2),int2str(trackIdx),'Color','w');
                    end
                    
                end
                
            end
            hold off
            
        end
        
        function imgOut = makeMovieFrame(nuclImg, trackerObj, nuclLabels, iT)
            %MAKEMOVIEFRAME  Make an annotated movie frame
            %
            %  F = FRETtracker.MAKEMOVIEFRAME(I, T, L, T) draws a frame of
            %  the movie from the nuclear image I, the TrackLink object T,
            %  the nuclear label L, and the timepoint T. The function
            %  returns the labelled image F as a matrix.
            %
            %  The outlines will be blue if the nucleus divided recently
            %  (within the last 5 frames), or green otherwise.
           
            %Overlay the nuclei outlines on the image as green lines
            nuclImg = double(nuclImg);
            imgOut = FRETtracker.showoverlay(FRETtracker.normalizeimg(imadjust(nuclImg)),bwperim(nuclLabels > 0),[0 1 0]);
            
            %If the nucleus divided recently (within the last 5 frames)
            %label the outline of the image red.
            for iTrack = 1:numel(trackerObj)
                trackIdx = trackerObj.trackingInfo(iTrack).Index;
                currTrack = trackerObj.getTrack(trackIdx);
                
                if currTrack.FirstFrame <= iT && currTrack.LastFrame >= iT
                    if ~isnan(currTrack.MotherIdx)
                        currMother = trackerObj.getTrack(currTrack.MotherIdx);
                        if (iT - currMother.LastFrame) <= 5
                            imgOut = FRETtracker.showoverlay(imgOut,bwperim(nuclLabels == currTrack.MaskID(end)),[1 0 0]);
                        end
                    end
                end
            end
            
        end
        
        function imgReader = getImageReader(filename)
            %GETIMAGEREADER  Return a reader for the image file
            %
            %  R = FRETtracker.GETIMAGEREADER returns a reader object R for
            %  image filename specified. R will be of class BioformatsImage
            %  if the file extension is .ND2, and a FRETtiffs object is
            %  the extension is .TIFF or .TIF.
            %
            %  See also: BioformatsImage, FRETtriffs
            
            [~,fname,ext] = fileparts(filename);
            switch lower(ext)
                
                case '.nd2'
                    imgReader = BioformatsImage(filename);
                    
                case {'.tiff','.tif'}
                    imgReader = TIFFreader(filename);
                    
                otherwise
                    error('FRETtracker:analyzeFile:UnsupportedFileExtension',...
                        '%s: Unsupported file extension %s.',...
                        fname, ext)
                    
            end
            
        end
        
        function frameData = getFrameData(iT, imageReader, nuclLabels, nuclImg, optionsStruct)
            %GETFRAMEDATA  Get frame data from the image reader
            %
            %  S = FRETtracker.GETFRAMEDATA(iT, imageReader, nuclLabels,
            %  optionsStruct) returns data from the current frame iT to the
            %  structure S.
            
            nuclData = regionprops(nuclLabels,...
                {'Area','Centroid','PixelIdxList'});
            
            %Make the cytoplasm ring if necessary
            if optionsStruct.measureCDK2 || ...
                    (optionsStruct.measureFRET && strcmpi(optionsStruct.FRETlocation, 'cyto'))
                
                %Make the cyto ring mask
                cytoRingMask = nuclLabels;
                cytoRingMask = imdilate(cytoRingMask, strel('disk',optionsStruct.CytoRingSize));
                cytoRingMask(cat(1,nuclData.PixelIdxList)) = 0;
                
            end
            
            
            if numel(nuclData) <= 0
                frameData = [];
                return;                
            end
            
               
            %% Initialize the data structure to the maximum number is can be
            %(i.e. the size of the regionprops output)
            frameData(numel(nuclData)) = struct('Centroid', [],...
                'Area', [],...
                'MaskID', [], ...
                'NuclInt', []);
            
            if optionsStruct.measureFRET
                frameData(numel(nuclData)).FRETInt = [];
                frameData(numel(nuclData)).CFPInt = [];
                frameData(numel(nuclData)).FRETratio = [];
                
                if optionsStruct.CorrectLocalBackground
                    frameData(numel(nuclData)).FRETIntCorr = [];
                    frameData(numel(nuclData)).CFPIntCorr = [];
                    frameData(numel(nuclData)).FRETratioCorr = [];
                end
            end
            
            if optionsStruct.measureCDK2
                frameData(numel(nuclData)).CDK2nucl = [];
                frameData(numel(nuclData)).CDK2cyto = [];
                frameData(numel(nuclData)).CDK2ratio = [];
                  
                if optionsStruct.CorrectLocalBackground
                    frameData(numel(nuclData)).CDK2nuclCorr = [];
                    frameData(numel(nuclData)).CDK2cytoCorr = [];
                    frameData(numel(nuclData)).CDK2ratioCorr = [];
                end
            end                                   
            
            %% Load the images we need
            if optionsStruct.measureFRET
                
                %Load the FRET and CFP image
                YFPImg = imageReader.getPlane(1, optionsStruct.FRETChan, iT);
                CFPImg = imageReader.getPlane(1, optionsStruct.CFPChan, iT);
                
                %Correct the local background (if options is set)
                if optionsStruct.CorrectLocalBackground
                    
                    YFPImgCorrected = FRETtracker.correctLocalBackground(YFPImg,...
                        optionsStruct.NumBackgroundBlocks,optionsStruct.BgPrctile);
                    CFPImgCorrected = FRETtracker.correctLocalBackground(CFPImg,...
                        optionsStruct.NumBackgroundBlocks,optionsStruct.BgPrctile);
                    
                end
            end
            
            if optionsStruct.measureCDK2
                
                %Load the CDK2 image
                CDK2Img = imageReader.getPlane(1, optionsStruct.CDK2Chan, iT);
                
                %Correct the local background (if options is set)
                if optionsStruct.CorrectLocalBackground
                    CDK2ImgCorrected = FRETtracker.correctLocalBackground(CDK2Img,...
                        optionsStruct.NumBackgroundBlocks,optionsStruct.BgPrctile);
                end
                
            end
                
            %% Populate the frame data structure
            idxToDelete = false(1,numel(nuclData));
            
            for iNucl = 1:numel(nuclData)
                
                %If the cell area is zero, skip it
                if nuclData(iNucl).Area == 0
                    idxToDelete(iNucl) = true;
                    continue;
                end
                
                %% Add the common data
                frameData(iNucl).Centroid = nuclData(iNucl).Centroid;
                frameData(iNucl).Area= nuclData(iNucl).Area;
                frameData(iNucl).NuclInt = mean(nuclImg(nuclData(iNucl).PixelIdxList));
                frameData(iNucl).MaskID = iNucl;
                
                %% Measure FRET data
                if optionsStruct.measureFRET
                    
                    switch optionsStruct.FRETlocation
                        
                        case 'nucl'
                            %Measure NLS FRET intensity
                            frameData(iNucl).FRETInt = mean(YFPImg(nuclData(iNucl).PixelIdxList));
                            frameData(iNucl).CFPInt = mean(CFPImg(nuclData(iNucl).PixelIdxList));
                            frameData(iNucl).FRETratio = frameData(iNucl).FRETInt / frameData(iNucl).CFPInt;
                            
                            %Background
                            if optionsStruct.CorrectLocalBackground
                                frameData(iNucl).FRETIntCorr = mean(YFPImgCorrected(nuclData(iNucl).PixelIdxList));
                                frameData(iNucl).CFPIntCorr = mean(CFPImgCorrected(nuclData(iNucl).PixelIdxList));
                                frameData(iNucl).FRETratioCorr = frameData(iNucl).FRETIntCorr / frameData(iNucl).CFPIntCorr;
                            end
                            
                        case 'cyto'
                            
                            %Measure NES FRET intensity
                            currMask = (cytoRingMask == iNucl);
                            
                            frameData(iNucl).FRETInt = mean(mean(YFPImg(currMask)));
                            frameData(iNucl).CFPInt = mean(mean(CFPImg(currMask)));
                            frameData(iNucl).FRETratio = frameData(iNucl).FRETInt / frameData(iNucl).CFPInt;
                            
                            %Measure the background corrected data
                            if optionsStruct.CorrectLocalBackground
                                frameData(iNucl).FRETIntCorr = mean(mean(YFPImgCorrected(currMask)));
                                frameData(iNucl).CFPIntCorr = mean(mean(CFPImgCorrected(currMask)));
                                frameData(iNucl).FRETratioCorr = frameData(iNucl).FRETIntCorr / frameData(iNucl).CFPIntCorr;
                            end

                    end
                    
                end
                %-- End measureFRET --%
                
                %% Measure CDK2 data if set
                
                if optionsStruct.measureCDK2
                    
                    %Measure data for the CDK2 cells
                    currMask = (cytoRingMask == iNucl);
                    frameData(iNucl).CDK2nucl = mean(CDK2Img(nuclData(iNucl).PixelIdxList));
                    frameData(iNucl).CDK2cyto = mean(mean(CDK2Img(currMask)));
                    frameData(iNucl).CDK2ratio = frameData(iNucl).CDK2cyto / frameData(iNucl).CDK2nucl;
                    
                    %Measure the background corrected data
                    if optionsStruct.CorrectLocalBackground
                        frameData(iNucl).CDK2nuclCorr = mean(CDK2ImgCorrected(nuclData(iNucl).PixelIdxList));
                        frameData(iNucl).CDK2cytoCorr = mean(mean(CDK2ImgCorrected(currMask)));
                        frameData(iNucl).CDK2ratioCorr = frameData(iNucl).CDK2cytoCorr / frameData(iNucl).CDK2nuclCorr;
                    end
                end
                %-- End measure CDK2 --%
                
            end
            
            %% Prune the data
            % There are empty datasets where the nuclear labels were
            % deleted during segmentation. Delete these.
            frameData(idxToDelete) = [];
                
        end
        
    end
    
    methods (Static)
        
        function nuclLabels = getNuclLabels(nuclImg, varargin)
            %GETNUCLLABELS  Labels the nuclear image
            %
            %  L = GETNUCLLABELS(I) segments and labels the foreground
            %  objects (i.e. cell nuclei) in image I.
            %
            %  L = GETNUCLLABELS(I, true) will set auto-contrast adjustment
            %  of the image before segmentation and labelling. The
            %  adjustment is carried out using 'imadjust', with no
            %  additional parameters.
            
            ip = inputParser;
            ip.addParameter('AutoAdjustNuclImg', false, @(x) islogical(x) || (isscalar(x) && isnumeric(x)) );
            ip.addParameter('NuclAreaRange', '', @(x) numel(x) == 2 && x(2) > x(1));
            ip.addParameter('MinNucleusArea', 0);
            ip.addParameter('MaxNucleusArea', Inf);
            ip.addParameter('ThresholdLvl', nan);
            ip.KeepUnmatched = true;
            ip.parse(varargin{:});
                        
            if ip.Results.AutoAdjustNuclImg
                %Automatically contast-adjust the image if set
                nuclImg = imadjust(nuclImg);
            end
                 
            if isempty(ip.Results.NuclAreaRange)
                nuclAreaRange = [ip.Results.MinNucleusArea, ip.Results.MaxNucleusArea];
            else
                nuclAreaRange = ip.Results.NuclAreaRange;
            end
            
            %Get threshold
            if isnan(ip.Results.ThresholdLvl)
                thLvl = FRETtracker.getThreshold(nuclImg);
            else
                thLvl = ip.Results.ThresholdLevel;
            end
            
            %Make the foreground mask
            binMask = nuclImg > thLvl;

            %binMask = activecontour(imadjust(nuclImg),binMask);
            binMask = imopen(binMask,strel('disk',2));
            binMask = imfill(binMask,'holes');
            binMask = bwareaopen(binMask,50);

            %Generate the distance transform
            dd = -bwdist(~binMask);
            dd(~binMask) = -Inf;
            
            %Surpress minima that are above the threshold level
            dd = imhmin(dd,1.2);
            
            %Run the watershed algorithm
            nuclLabels = watershed(dd);
            
            %Remove labels which are intersecting the image border, or are
            %too small/large
            nuclLabels = imclearborder(nuclLabels);
            
            maskValidSize = bwareaopen(nuclLabels, nuclAreaRange(1));
            mask_tooLarge = bwareaopen(nuclLabels, nuclAreaRange(2));
            maskValidSize(mask_tooLarge) = 0;
                        
            nuclLabels(~maskValidSize) = 0;
            
%             FRETtracker.showoverlay(FRETtracker.normalizeimg(imadjust(nuclImg)),bwperim(nuclLabels),[0 1 0]);
%             keyboard
            
        end
                
        function thLvl = getThreshold(imageIn)
            %GETTHRESHOLD  Get a threshold for the image
            %
            %  T = FRETtracker.GETTHRESHOLD(I) gets a greyscale threshold
            %  level T for the image I.
            %
            %  Threshold is determined by looking at image histogram, then
            %  looking for the greyscale value where the maximum count
            %  drops to at least 20%.
            
            %Get the image intensity histogram
            binEdges = linspace(0,double(max(imageIn(:))),200);
            [nCnts, binEdges] = histcounts(imageIn(:),binEdges);
            binCenters = diff(binEdges) + binEdges(1:end-1);
            
            nCnts = smooth(nCnts,5);
            
            %Find the background peak count
            [bgCnt,bgLoc] = findpeaks(nCnts,'Npeaks',1,'SortStr','descend');
            
            %Find where the histogram counts drops to at least 20% of this value
            thLoc = find(nCnts(bgLoc:end) <= bgCnt * 0.01,1,'first');
            
            if isempty(thLoc)
                error('FRETtracker:getThreshold:CouldNotGetThreshold',...
                    'Auto thresholding failed to find a suitable threshold level. Try specifying one manually.');                
            end
            
            thLvl = binCenters(thLoc + bgLoc);
            
%             plot(binCenters,nCnts,binCenters(bgLoc),bgCnt,'x',[thLvl, thLvl],ylim,'r--');
%             keyboard
%             
%             
            
        end
        
        function [correctedImg, maskedImg, bgImage] = correctLocalBackground(imageIn,numBlocks,bgPrctile)
            %CORRECTLOCALBACKGROUND  Apply local background correction
            %
            %  C = FRETtracker.CORRECTLOCALBACKGROUND(I, N, P) runs the
            %  local background correction algorithm on the image I,
            %  returning the background-corrected image C. N is the number
            %  of blocks (N = [Nrows, Ncols]) to divide the image into. P
            %  is the percentile used to estimate the background.
            %
            %  During image acquisition, uneven illumination due to
            %  microscope optics makes the objects at the center of the
            %  image appear brighter than objects at the edge. 
            %
            %  To correct for this, we use a method similar to the Spencer
            %  lab: The image is divided into a number of blocks. The local
            %  background is estimated by calculating the lowest percentile
            %  of intensities in the block.
                        
            %Take the median filter of the image to remove hotspots
            imageInTemp = double(medfilt2(imageIn,[3 3]));
            
            %Calculate the indices for each block
            blockHeight = floor(size(imageInTemp,1)/numBlocks(1));
            blockWidth = floor(size(imageInTemp,2)/numBlocks(2));
            
            blockRowIdxs = 1:blockHeight:size(imageInTemp,1);
            blockRowIdxs(end) = size(imageInTemp,1);
            
            blockColIdxs = 1:blockWidth:size(imageInTemp,2);
            blockColIdxs(end) = size(imageInTemp,2);
            
            maskedImg = zeros(size(imageInTemp));
            
            %Make the background image
            temp = zeros(size(imageIn));
            
            bgImage = zeros(size(imageInTemp));
            for iRow = 1:numBlocks(1)
                for iCol = 1:numBlocks(2)
                    
                    croppedImg = imageInTemp(blockRowIdxs(iRow):blockRowIdxs(iRow+1),...
                        blockColIdxs(iCol):blockColIdxs(iCol+1));
                    
                    temp(blockRowIdxs(iRow):blockRowIdxs(iRow+1),...
                        blockColIdxs(iCol):blockColIdxs(iCol+1)) = croppedImg;
                    
                    bgValue = prctile(croppedImg(:),bgPrctile);
                    
                    bgImage(blockRowIdxs(iRow):blockRowIdxs(iRow+1),...
                        blockColIdxs(iCol):blockColIdxs(iCol+1)) = ...
                        ones(size(croppedImg)) .* bgValue;
                    
                    maskedImg(blockRowIdxs(iRow):blockRowIdxs(iRow+1),...
                        blockColIdxs(iCol):blockColIdxs(iCol+1)) = ...
                        croppedImg <= bgValue;
                end
            end
            
            correctedImg = double(imageIn) - bgImage;
        end
        
        function varargout = showoverlay(baseimage, mask, color, varargin)
            %SHOWOVERLAY    Plot an overlay mask on an image
            %
            %  FRETtracker.SHOWOVERLAY(IMAGE,MASK,COLOR) will plot an
            %  overlay specified by a binary MASK on the IMAGE. The color
            %  of the overlay is specified using a three element vector
            %  COLOR.
            %
            %  O = FRETtracker.SHOWOVERLAY(IMAGE,MASK,COLOR) will return
            %  the overlaid image to the variable O.
            
            if ~exist('color','var')
                color = [1 1 1]; %Default color of the overlay
            end
            
            if size(baseimage,3) == 3
                red = baseimage(:,:,1);
                green = baseimage(:,:,2);
                blue = baseimage(:,:,3);
                
            elseif size(baseimage,3) == 1
                red = baseimage;
                green = baseimage;
                blue = baseimage;
                
            else
                error('Image should be either NxNx1 (greyscale) or NxNx3 (rgb)')
            end
            
            %Make sure the mask is binary (anything non-zero becomes true)
            mask = (mask ~= 0);
            
            if isinteger(baseimage)
                maxInt = intmax(class(baseimage));
            else
                maxInt = 1;
            end
            
            red(mask) = color(1) .* maxInt;
            green(mask) = color(2) .* maxInt;
            blue(mask) = color(3) .* maxInt;
            
            %Concatenate the output
            outputImg = cat(3,red,green,blue);
            
            if nargout == 0
                imshow(outputImg,[])
            else
                varargout{1} = outputImg;
            end
        
        end
    
        function imageOut = normalizeimg(imageIn,varargin)
            %NORMALIZEIMG   Linear dynamic range expansion for contrast enhancement
            %   N = NORMALIZEIMG(I) expands the dynamic range (or contrast) of image I
            %   linearly to maximize the range of values within the image.
            %
            %   This operation is useful when enhancing the contrast of an image. For
            %   example, if I is an image with uint8 format, with values ranging from
            %   30 to 100. Normalizing the image will expand the values so that they
            %   fill the full dynamic range of the format, i.e. from 0 to 255.
            %
            %   The format of the output image N depends on the format of the input
            %   image I. If I is a matrix with an integer classs (i.e. uint8, int16), N
            %   will returned in the same format. If I is a double, N will be
            %   normalized to the range [0 1] by default.
            %
            %   N = NORMALIZEIMG(I,[min max]) can also be used to specify a desired
            %   output range. For example, N = normalizeimg(I,[10,20]) will normalize
            %   image I to have values between 10 and 20. In this case, N will be
            %   returned in double format regardless of the format of I.
            %
            %   In situations where most of the interesting image features are
            %   contained within a narrower band of values, it could be useful to
            %   normalize the image to the 5 and 95 percentile values.
            %
            %   Example:
            %       I = imread('cameraman.tif');
            %
            %       %Calculate the values corresponding to the 5 and 95 percentile of
            %       %values within the image
            %       PRC5 = prctile(I(:),5);
            %       PRC95 = prctile(I(:),95);
            %
            %       %Threshold the image values to the 5 and 95 percentiles
            %       I(I<PRC5) = PRC5;
            %       I(I>PRC95) = PRC95;
            %
            %       %Normalize the image
            %       N = normalizeimg(I);%
            %
            %       %Display the normalized image
            %       imshow(N)
            
            %Define default output value range
            outputMin = 0;
            outputMax = 1;
            
            %Check if the desired output range is set. If it is, make sure it contains
            %the right number of values and format, then update the output minimum and
            %maximum values accordingly.
            if nargin >= 2
                if numel(varargin{1}) ~= 2
                    error('The input parameter should be [min max]')
                end
                
                outputMin = varargin{1}(1);
                outputMax = varargin{1}(2);
            else
                %If the desired output range is not set, then check if the image is an
                %integer class. If it is, then set the minimum and maximum values
                %to match the range of the class type.
                if isinteger(imageIn)
                    inputClass = class(imageIn);
                    
                    outputMin = 0;
                    outputMax = double(intmax(inputClass)); %Get the maximum value of the class
                    
                end
            end
            
            %Convert the image to double for the following operations
            imageIn = double(imageIn);
            
            %Calculate the output range
            outputRange = outputMax - outputMin;
            
            %Get the maximum and minimum input values from the image
            inputMin = min(imageIn(:));
            inputMax = max(imageIn(:));
            inputRange = inputMax - inputMin;
            
            %Normalize the image values to fit within the desired output range
            imageOut = (imageIn - inputMin) .* (outputRange/inputRange) + outputMin;
            
            %If the input was an integer before, make the output image the same class
            %type
            if exist('inputClass','var')
                eval(['imageOut = ',inputClass,'(imageOut);']);
            end
            
        end
                
    end
    
end