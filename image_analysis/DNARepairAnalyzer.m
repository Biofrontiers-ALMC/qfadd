classdef DNARepairAnalyzer < handle
    %DNAREPAIRANALYZER  Class for analyzing DNA repair movies
    %
    %  The DNAREPAIRANALYZER class contains the code for processing and
    %  analyzing DNA repair movies.
    %
    %  The following is an example to run this code:
    %
    %     %Initialize a DNAREPAIRANALYZER object
    %     analyzerObj = DNAREPAIRANALYZER;
    %
    %     %You can run the processing on a single file or a whole directory
    %     data = analyzerObj.processFiles('filename.nd2');

    properties
        
        NuclearChannel = 'TRITC';
        NumFramesIrradiation = 6;
        ROIgrowByPixels = [-10,0];
        
        CorrectPhotobleaching = true;
        ThresholdLevel = 'auto';
        BackgroundCorrection = true;
        
        AdditionalROIs = 2;
        
        AutosaveData = true;
        SaveVideo = true;
        SaveROIdata = true;
        
        CellsCanTouchSides = false;
        
    end
    
    methods
        
        function processFiles(obj, files, outputDir)
            %PROCESSFILES  Collect intensity data from file(s)
            %
            %  PROCESSFILES(OBJ) will open a dialog box for the user to
            %  select the file(s) to process. Multiple files can be
            %  selected.
            %
            %  PROCESSFILES(OBJ, FILES, OUTPUTDIR) is an alternative
            %  way to call this method. FILES can be a cell of strings
            %  containing full pathnames to process.
            
            %---Parse inputs---
            %Get files to process
            if ~exist('files', 'var')
                [files, fpath] = uigetfile({'*.nd2; *.tif; *.tiff', 'Supported images (*.nd2; *.tif; *.tiff)'},...
                    'MultiSelect', 'on');
                
                if fpath == 0
                    return;                    
                end
                
                if ~iscell(files)
                    files = {files};
                end
                
                for iF = 1:numel(files)
                    files{iF} = fullfile(fpath, files{iF});
                end
                
            else
                if ~iscell(files)
                    files = {files};
                end
            end
            
            if ~exist('outputDir','var')
                
                %Prompt the user to select an input directory
                outputDir = uigetdir(fileparts(files{1}),'Choose the output data directory');
                
                if outputDir == 0
                    return;
                end
            else
                if ~exist(outputDir,'dir')
                    mkdir(outputDir);
                end
                
            end
            
            %---Start processing---
            %Save the settings file
            exportSettings(obj, fullfile(outputDir, 'settings.txt'));
            
            opts = props2struct(obj);
            
            for iFile = 1:numel(files)
                
                fprintf('%s: Processing %s\n', datestr(now), files{iFile});
                
                try
                    DNARepairAnalyzer.processVideo(files{iFile}, outputDir, opts);
                    fprintf('%s: Completed %s\n', datestr(now), files{iFile});
                catch ME
                    fprintf('%s: Error processing %s\n', datestr(now), files{iFile})
                    fprintf('%s\n', getReport(ME));
                    continue;
                end
 
            end
            
        end
        
        function importSettings(obj, varargin)
            %IMPORTSETTINGS  Import settings from file
            %
            %  IMPORTSETTINGS(OBJ, FILENAME) will load the options from the
            %  file specified.
            %
            %  IMPORTSETTINGS(OBJ) will open a dialog box for the user to
            %  select the option file.
            
            %Get the options file
            if isempty(varargin)                
                [fname, fpath] = uigetfile({'*.txt','Text file (*.txt)';...
                    '*.*','All files (*.*)'},...
                    'Select settings file');
                
                optionsFile = fullfile(fpath,fname);
                
            elseif numel(varargin) == 1
                optionsFile = varargin{1};
                
            else
                error('DNARepairAnalyzer:TooManyInputs', 'Too many input arguments.');
                
            end
            
            fid = fopen(optionsFile,'r');
            
            if isequal(fid,-1)
                error('DNARepairAnalyzer:ErrorReadingFile',...
                    'Could not open file %s for reading.',fname);
            end
            
            ctrLine = 0;
            while ~feof(fid)
                currLine = strtrim(fgetl(fid));
                ctrLine = ctrLine + 1;
                
                if isempty(currLine) || strcmpi(currLine(1),'%') || strcmpi(currLine(1),'#')
                    %Empty lines should be skipped
                    %Lines starting with '%' or '#' are comments, so ignore
                    %those
                    
                else
                    
                    parsedLine = strsplit(currLine,'=');
                    
                    %Check for errors in the options file
                    if numel(parsedLine) < 2 
                        error('DNARepairAnalyzer:ErrorReadingOption',...
                            'Error reading <a href="matlab:opentoline(''%s'', %d)">file %s (line %d)</a>',...
                            optionsFile, ctrLine, optionsFile, ctrLine);
                    elseif isempty(parsedLine{2})
                        error('DNARepairAnalyzer:ErrorReadingOption',...
                            'Missing value in <a href="matlab:opentoline(''%s'', %d)">file %s (line %d)</a>',...
                            optionsFile, ctrLine, optionsFile, ctrLine);
                    end
                    
                    %Get parameter name (removing spaces)
                    parameterName = strtrim(parsedLine{1});
                    
                    %Get value name (removing spaces)
                    value = strtrim(parsedLine{2});
                    
                    if isempty(value)
                        %If value is empty, just use the default
                    else
                        obj.(parameterName) = eval(value);
                    end
                    
                end
                
            end
            
            fclose(fid);
        end  
        
        function exportSettings(obj, varargin)
            %EXPORTSETTINGS  Export object properties into a text file
            %
            %  EXPORTSETTINGS(OBJ, FILENAME) will export the settings into
            %  a text file.
            
            if isempty(varargin)
                [fname, fpath] = uiputfile({'.txt','Text file (.txt)'; '*.*', 'All files (*.*)'});
                
                if fname == 0
                    return;
                end
                
                fn = fullfile(fpath, fname);
            else
                fn = varargin{1};
            end
          
            props = properties(obj);
            
            fid = fopen(fn,'w');
            
            fprintf(fid, '%%DNARepairAnalyzer Settings\r\n%%%s\r\n\r\n',datestr(now));
            
            for iP = 1:numel(props)
                
                if ischar(obj.(props{iP}))
                    fprintf(fid,'%s = ''%s'' \r\n',props{iP}, ...
                        obj.(props{iP}));
                    
                elseif isnumeric(obj.(props{iP}))
                    fprintf(fid,'%s = %s \r\n',props{iP}, ...
                        mat2str(obj.(props{iP})));
                    
                elseif islogical(obj.(props{iP}))
                    
                    if obj.(props{iP})
                        fprintf(fid,'%s = true \r\n',props{iP});
                    else
                        fprintf(fid,'%s = false \r\n',props{iP});
                    end
                    
                end
                
            end
            
            fclose(fid);
        end
        
    end

    methods (Static)
        
        function structOut = processVideo(filename, outputDir, opts)
            %PROCESSVIDEO  Process the specified file
            %
            %  DNARepairAnalyzer.PROCESSVIDEO(FILENAME, OUTPUTDIR, OPTS)
            %  will process the file specified, saving the data in
            %  OUTPUTDIR. OPTS should be a struct containing the settings.
            %  OPTS can be generated using the hidden function
            %  props2struct(OBJ).
            %
            %  This command supports .ND2 and .TIF.
            
            %Load images into memory. The output is a cell containing
            %matrices with image data. The number of cell elements should
            %equal the number of channels, while the size of the 3rd dim of
            %each matrix equals the number of timepoints.
            [~,~,fExt] = fileparts(filename);
            switch lower(fExt)
                
                case {'.tif', '.tiff'}
                    %TIF files are expected to have channels interleaved
                    
                    tiffInfo = imfinfo(filename);
                    
                    imgHeight = tiffInfo.Height;
                    imgWidth = tiffInfo.Width;
                    
                    sizeT = numel(tiffInfo) / 2;
                    
                    sizeC = 2; %Expect only two channels
                    
                    imgData = cell(1, sizeC);
                    
                    channelNames = cell(1, sizeC);
                    for iC = 1:sizeC
                        
                        imgData{iC} = zeros(imgHeight, imgWidth, sizeT);
                        
                        for iT = 1:sizeT
                            %Calculate the imgInd
                            imgInd = (iT - 1) * 2 + iC;
                            imgData{iC}(:,:,iT) = imread(filename, imgInd);
                        end
                        
                        channelNames{iC} = sprintf('Channel%d',iC);
                    end
                                    
                    dataOut.timestamps = ((1:sizeT) - 1) * 0.140; %Obtained from LIF file
                    dataOut.timeUnit = 's';
                    
                case '.nd2'  %Nikon files
                    
                    %Load all the images into memory
                    bfr = BioformatsImage(filename);
                    
                    imgData = cell(1, bfr.sizeC);
                    
                    for iC = 1:bfr.sizeC
                        
                        imgData{iC} = zeros(bfr.height, bfr.width, bfr.sizeT, 'uint16');
                        
                        for iT = 1:bfr.sizeT
                            imgData{iC}(:,:,iT) = bfr.getPlane(1, iC, iT);
                        end
                    end
                    
                    channelNames = regexprep(bfr.channelNames,'\W', '');
                    
                    %If the channel name in opts was a string, convert it
                    %into an index
                    if ischar(opts.NuclearChannel)
                        chanInd = find(strcmpi(opts.NuclearChannel, channelNames));               
                        
                        if isempty(chanInd)
                            error('Could not find channel %s.', opts.NuclearChannel);
                        else
                            opts.NuclearChannel = chanInd;
                        end
                    end
                    
                    %Get timestamps
                    [dataOut.timestamps,dataOut.timeUnit] = bfr.getTimestamps(1,1);
                    clear bfr 
                    
                otherwise
                    error('%s is not a supported file format.', fExt);
            end
                        
            %--- Collect data ---%
            
            %Initialize the data output structure
            dataOut.filename = filename;
            
            for iC = 1:numel(channelNames)
                dataOut.(channelNames{iC}).RawNucleusIntensity = zeros(1, size(imgData{1},3));
            end
            
            %Load the ROI
            roiImg = imread([filename(1:end-4),'_ROI.tif']);
            
            %Grow the ROI by the size specified
            roiImg = DNARepairAnalyzer.growROI(roiImg,opts.ROIgrowByPixels);
            
            %Generate a common save filename and start the video (if set)
            if opts.AutosaveData
                [~,saveFN] = fileparts(filename);
                
                if opts.SaveVideo
                    vidWriter = VideoWriter(fullfile(outputDir,[saveFN,'.avi']));
                    vidWriter.FrameRate = 20;
                    open(vidWriter);
                end
            end
            
            %---Begin processing---%
            
            %Redefine image properties
            sizeC = numel(imgData);
            sizeT = size(imgData{1}, 3);
            
            %Get the nuclear mask from the first frame
            nuclImg = imgData{opts.NuclearChannel}(:, :, 1);
            nuclMask = DNARepairAnalyzer.getNuclMask(nuclImg,'ThresholdLevel',opts.ThresholdLevel, 'CellsCanTouchSides', opts.CellsCanTouchSides);
            
            %----- Determine normalization values -----%
            
            %Calculate the background value
            if opts.BackgroundCorrection
                
                %Initializet output data
                dataOut.backgroundLvl = zeros(1, sizeC);
                
                %Calculate the background level for each channel
                for iC = 1:sizeC
                    bgValues = zeros(1,opts.NumFramesIrradiation);
                    
                    for iFrame = 1:opts.NumFramesIrradiation
                        
                        currImg = imgData{iC}(:,:,1);
                        
                        %Segment the background using the first frame
                        if iFrame == 1
                            [nCnts, xEdges] = histcounts(currImg(:),'BinLimits',[0, max(currImg(:))]);
                            nCnts = smooth(nCnts,3);
                            xCenters = diff(xEdges) + xEdges(1:end-1);
                            
                            %Find the first peak
                            [pkCnts,bgPkLoc] = findpeaks(nCnts, 'NPeaks',1);
                            
                            %Find where the value drops to about 67% (1/e)
                            %of the peak
                            tempCnts = nCnts(bgPkLoc:end);
                            [~,bgThresholdLoc] = min(abs(tempCnts - pkCnts/exp(1)));
                            bgThresholdLoc = bgThresholdLoc + bgPkLoc;
                            
                            %Segment the background from the first frame
                            bgMask = currImg < xCenters(bgThresholdLoc);
                            
                        end
                        
                        %Calculate the mean intensity in the background ROI
                        bgValues(iFrame) =  mean(currImg(bgMask));
                    end
                    
                    %Store the background value data in the output
                    %data structure
                    dataOut.backgroundLvl(iC) = mean(bgValues);
                end
            end
            
            if opts.CorrectPhotobleaching
                %To correct for photobleaching, we will take the
                %average of the total nuclear intensity for the
                %first seven frames, then use this value to
                %normalize the total nuclear intensity of all
                %subsequent frames.
                
                nuclInt = zeros(sizeC, opts.NumFramesIrradiation);
                
                for iFrame = 1:opts.NumFramesIrradiation
                    for iC = 1:sizeC
                        currImg = imgData{iC}(:,:,iFrame);
                        nuclInt(iC,iFrame) =  sum(sum(currImg(nuclMask)));
                    end
                end
                
                dataOut.phbCorrection = mean(nuclInt,2);
                
            else
                %If no corrections, set the correction factor to 1
                dataOut.phbCorrection = ones(bfr.sizeC,1);
                
            end
            
            %------ Measure the data  -------
            for iT = 1:sizeT                
                for iC = 1:sizeC
                    
                    %Get the image plane
                    currImgOriginal = imgData{iC}(:,:,iT);
                    
                    %Subtract the background (if set)
                    if opts.BackgroundCorrection
                        currImgOriginal = currImgOriginal - dataOut.backgroundLvl(iC);
                    end
                    
                    %Photobleaching correction
                    if opts.CorrectPhotobleaching
                        %Store the sum raw intensity in the nucleus
                        dataOut.(channelNames{iC}).RawNucleusIntensity(iT) = sum(sum(currImgOriginal(nuclMask)));
                        
                        normFactor = dataOut.phbCorrection(iC) ./ dataOut.(channelNames{iC}).RawNucleusIntensity(iT);
                        currImg = currImgOriginal .* normFactor ;
                        
                        %Store the normalized intensity (this should hopefully
                        %be 1)
                        dataOut.(channelNames{iC}).NormNucleusIntensity(iT) = sum(sum(currImg(nuclMask)));
                    else
                        currImg = currImgOriginal;
                    end
                    
                    %Store the sum raw intensity in the ROI
                    dataOut.(channelNames{iC}).RawIntensity(iT) = sum(sum(currImg(roiImg)));
                    
                    if opts.AutosaveData && opts.SaveVideo
                        imgOut = DNARepairAnalyzer.showoverlay(DNARepairAnalyzer.normalizeimg(currImg),bwperim(nuclMask),[0 1 0]);
                        imgOut = DNARepairAnalyzer.annotateImage(imgOut, roiImg, 'ROI 0');
                    end
                    
                    if opts.AdditionalROIs > 0
                        %For each additional ROI:
                        % 1. Move the main ROI up and down
                        % 2. Measure the intensity and append it to the
                        % approriate field in the output struct
                        %TODO: Check that the number of ROIs do not exceed
                        %the image size
                        
                        roiVert = any(roiImg,2);
                        roiHeight = sum(roiVert);
                        
                        for iAR = 1:opts.AdditionalROIs
                            %Each ROI should be the same size as the main ROI.
                            %So the offset should be the same as the number of
                            %ROI rows + 1.
                            %
                            %In the future, if we need automatic determination
                            %of displacement direction, then it could be the
                            %shorter axis (or user specified).
                            roiOffset = iAR * roiHeight;
                            
                            roiShiftUp = circshift(roiImg,[-roiOffset,0]);
                            dataOut.(channelNames{iC}).ROIAbove(iAR).RawIntensity(iT) = sum(sum(currImg(roiShiftUp)));
                            
                            roiShiftDown = circshift(roiImg,[roiOffset,0]);
                            dataOut.(channelNames{iC}).ROIBelow(iAR).RawIntensity(iT) = sum(sum(currImg(roiShiftDown)));
                            
                            if opts.AutosaveData && opts.SaveVideo
                                imgOut = DNARepairAnalyzer.annotateImage(imgOut, roiShiftUp, sprintf('ROI %d',iAR));
                                imgOut = DNARepairAnalyzer.annotateImage(imgOut, roiShiftDown, sprintf('ROI %d',-iAR));
                            end
                            
                        end
                    end
                end
                
                if opts.AutosaveData && opts.SaveVideo
                    
                    imgOut = double(imgOut);
                    imgOut = imgOut ./ max(imgOut(:));
                    
                    vidWriter.writeVideo(imgOut);
                end
            end
                     
            %---Normalize the intensity to the base fluorescence level (before irradiation)---
            for iC = 1:sizeC
                
                baseValue = mean(dataOut.(channelNames{iC}).RawIntensity(1:opts.NumFramesIrradiation));
                
                dataOut.(channelNames{iC}).NormIntensity = ...
                    dataOut.(channelNames{iC}).RawIntensity./baseValue;
                
                if opts.AdditionalROIs > 0
                    for iAR = 1:opts.AdditionalROIs
                        dataOut.(channelNames{iC}).ROIAbove(iAR).NormIntensity = ...
                            dataOut.(channelNames{iC}).ROIAbove(iAR).RawIntensity./...
                            mean(dataOut.(channelNames{iC}).ROIAbove(iAR).RawIntensity(1:opts.NumFramesIrradiation));
                        
                        dataOut.(channelNames{iC}).ROIBelow(iAR).NormIntensity = ...
                            dataOut.(channelNames{iC}).ROIBelow(iAR).RawIntensity./...
                            mean(dataOut.(channelNames{iC}).ROIBelow(iAR).RawIntensity(1:opts.NumFramesIrradiation));
                        
                    end
                end
            end
            
            %--- Save data (if set) ---%
            if opts.AutosaveData
                
                fh = figure('Visible','off');
                plot(dataOut.timestamps,dataOut.(channelNames{end}).NormIntensity)
                
                if opts.AdditionalROIs > 0
                    for iAR = 1:opts.AdditionalROIs
                        hold on
                        plot(dataOut.timestamps,dataOut.(channelNames{end}).ROIAbove(iAR).NormIntensity)
                        plot(dataOut.timestamps,dataOut.(channelNames{end}).ROIBelow(iAR).NormIntensity)
                        hold off
                    end
                end
                
                title(channelNames{end});
                roiRange = [-opts.AdditionalROIs:-1, 1:opts.AdditionalROIs];
                
                %Show the legend for each of the ROIs
                legendStr = cell(numel(roiRange)+1,1);
                legendStr{1} = 'ROI 0';
                for iROI = 1:numel(roiRange)
                    legendStr{iROI + 1} = sprintf('ROI %d',roiRange(iROI));
                end
                legend(legendStr)
                
                xlabel('Time')
                ylabel('Normalized intensity')
                
                print(fh,fullfile(outputDir,[saveFN,'.png']),'-dpng');
                close(fh)
                
                DNARepairAnalyzer.exportToCSV(fullfile(outputDir,[saveFN,'.csv']),dataOut);
                
                DNARepairAnalyzer.exportInfo(outputDir, saveFN, roiImg, nuclMask)
                
                save(fullfile(outputDir,[saveFN,'.mat']), 'dataOut');
                
                if opts.SaveVideo
                    close(vidWriter);
                end
                
                fprintf('Processed data saved successfully to ''%s''.\n',outputDir);
                
            end
            
            dataOut = [];
            
            if nargout > 0
                structOut = dataOut;
            end
            
        end
        
        function filename = getROIfile(initPath)
            %GETROIFILE  Prompts the user to select the ROI file
            
            if ~exist('initPath','var')
                initPath = '';
            end
            
            %Prompt the user so they know that the ROI file was not found
            btn = questdlg('ROI file was not found. Do you want to search for the file yourself?',...
                'ROI file not found','Yes','No','Yes');
            
            if strcmpi(btn,'yes')
                [filename, pathname] = uigetfile( ...
                    {'*.TIF, *.TIFF, *.GIF, *.JPG, *.JPEG', ...
                    'Image file (*.tiff,*.gif,*.jpeg)'}, ...
                    'Select the ROI file',initPath);
                
                filename = fullfile(pathname,filename);
            else
                filename = 0;
            end
        end
        
        function irradFrame = getIrradiationFrame(bfReaderIn)
            %GETIRRADIATIONFRAMELENGTH  Try to determine the irradiation
            %frame
            %
            %There is a slight pause in the movie while the nucleus is
            %irradiated (damaged). This function looks for the longest time
            %delta between frames and chooses that as the irradiation
            %frame.
            
            if ~isa(bfReaderIn,'BioformatsImage')
                error('Expected input to be a BioformatsImage.')
            end
            
            timeVec = bfReaderIn.getTimestamps(1,1);
            
            dTime = [0, diff(timeVec)];
            
            [~,irradFrame] = max(dTime);
            
        end
        
        function roiOut = growROI(roiIn,growSpec)
            %GROWROI  Grow the initial ROI
            %
            %  M = DNARepairAnalyzer.GROWROI(R, [width, height]) changes
            %  the original mask R by the width and height specified,
            %  returned in the logical array M.
            %
            %  The width and height should be specified in pixels. Negative
            %  values will shrink the mask in the specified dimension.
            %
            %  This function assumes that input mask R is rectangular. The
            %  code gets the coordinates of the rectangle, then makes a new
            %  one with the additional width and height specified.
            %
            %  Example:
            %
            %    %Shrink the mask by 10 pixels in height, and grow it by 5
            %    %pixels in width.
            %    roiOut = DNARepairAnalyzer(roiIn, [5, -10]);
            
            rowStart = find(any(roiIn,2),1,'first') - growSpec(2);
            rowEnd = find(any(roiIn,2),1,'last') + growSpec(2);
            
            colStart = find(any(roiIn,1),1,'first') - growSpec(1);
            colEnd = find(any(roiIn,1),1,'last') + growSpec(1);
            
            roiOut = false(size(roiIn));
            roiOut(rowStart:rowEnd,colStart:colEnd) = true;
            
        end
        
        function [nuclMask, thLvl] = getNuclMask(nuclImg,varargin)
            %GETNUCLMASK  Segment the nucleus
            %
            %  M = DNARepairAnalyzer.GETNUCLMASK(I) will segment the image
            %  I, returning a logical array M.
            %
            %  The threshold level is calculated using the intensity at the
            %  center of the image.
            %
            %  Optional arguments:
            %     'stdLevel' = Number of standard deviations for threshold
            %     (Default = 3)
            %
            %     'circRadius' = Radius near center to calculate values
            %     (Default = 5 px)
            
            ip = inputParser;
            ip.addParameter('ThresholdLevel','auto',@(x) ischar(x) || isnumeric(x));
            ip.addParameter('CircRadius',20,@(x) isscalar(x));
            ip.addParameter('CellsCanTouchSides', false, @(x) islogical(x));
            ip.parse(varargin{:});
            
            %Convert the image to double
            nuclImg = double(nuclImg);

            if strcmpi(ip.Results.ThresholdLevel,'auto')
                tempNuclImg = medfilt2(nuclImg, [30 30]);
                
                %Determine the threshold level as mean - 0.5 * mean of the
                %intensity in a 5 pixel circle around the center of the image
                centerMask = false(size(tempNuclImg));
                xx = 1:size(tempNuclImg,2);
                xx = xx - median(xx);
                yy = 1:size(tempNuclImg,1);
                yy = yy - median(yy);
                [xx,yy] = meshgrid(xx,yy);
                centerMask(xx.^2 + yy.^2 < ip.Results.CircRadius^2) = true;
                
                %Calculate the mean and std of the intensity in the mask
                meanInt = mean(tempNuclImg(centerMask));
                stdInt = std(tempNuclImg(centerMask));
                
                %Calculate the threshold level
                %thLvl = (meanInt - ip.Results.stdLevel * stdInt);
                thLvl = (meanInt - 0.5 * meanInt);
                
                %If the threshold level is too low (below zero), automatically
                %adjust it until it is above the mean intensity level
                if thLvl <= 0
                    factor = 0.5;
                    while thLvl <= meanInt
                        %Try reducing the stdLevel
                        factor = factor - 0.1;
                        thLvl = (meanInt -  factor * stdInt);
                    end
                    warning('DNARepairAnalyzer:ThresholdTooLow',...
                        'Initial threshold level %.2f was too low and was raised to %.2f',...
                        (meanInt - 0.5 * meanInt),thLvl)
                end
            else
                thLvl = ip.Results.ThresholdLevel;
            end
            
            %Generate the mask
            objMask = nuclImg > thLvl;

            if ~ip.Results.CellsCanTouchSides
                objMask = imclearborder(objMask);
            end            
            
            objMask = bwareaopen(objMask,500);
            objMask = imfill(objMask,'holes');
            objMask = imdilate(objMask,strel('disk',3));
            
            %Keep only objects at the center of the image
            detObjs = regionprops(objMask,'Centroid','PixelIdxList');
            
            if numel(detObjs) > 1
                %Find the object closest to the center of the image
                centroids = cat(1,detObjs.Centroid);
                
                imgCenter = size(nuclImg)./2;
                distances = sqrt((centroids(:,1) - imgCenter(1)).^2 + (centroids(:,2) - imgCenter(2)).^2);
                [~,minObjInd] = min(distances);
                
                nuclMask = false(size(nuclImg));
                nuclMask(detObjs(minObjInd).PixelIdxList) = true;
                
            else
                nuclMask = objMask;
            end
            
            %Make sure all holes are filled
            nuclMask = imfill(nuclMask,'holes');
            
            if all(~nuclMask(:))
                error('DNARepairAnalyzer:getNuclMask:CouldNotFindNucleus',...
                    'Error finding nucleus. Please specify a threshold level manually.');
            end

        end
        
        function nuclPosition = getNuclPosition(nuclMask)
            %GETNUCLPOSITION  Get the centroid of the nuclear mask

            objProps = regionprops(nuclMask,'Centroid');
            
            nuclPosition = objProps(closestObjIdx).Centroid;
            
        end
        
        function exportToCSV(filename, dataIn)
            
            %CSV OUTPUT:
            % Filename:
            % Date:
            % Time, ROI -N, ..., ROI 0, ..., ROI N
            
            %Make sure the filename ends in a CSV
            [fpath,fname,fext] = fileparts(filename);
            if ~strcmpi(fext,'.csv')
               filename = fullfile(fpath,[fname,'.csv']);
            end
            
            numRows = numel(dataIn.timestamps);
            
            channelNames = fieldnames(dataIn);
            if ismember('backgroundLvl',channelNames)
                channelNames(strcmpi(channelNames,'backgroundLvl')) = [];
            end
            
            %List of fields to ignore. This function should be done in a
            %better way.
            channelNames(strcmpi(channelNames,'filename')) = [];
            channelNames(strcmpi(channelNames,'timestamps')) = [];
            channelNames(strcmpi(channelNames,'timeUnit')) = [];
            channelNames(strcmpi(channelNames,'phbCorrection')) = [];
            channelNames(strcmpi(channelNames,'PreIrradiationFrame')) = [];
            numChannels = numel(channelNames);
            
            if isfield(dataIn.(channelNames{1}),'ROIAbove')
                numAdditionalROIs = numel(dataIn.(channelNames{1}).ROIAbove);
            else
                numAdditionalROIs = 0;
            end
            
            fid = fopen(filename,'w');
            if fid == -1
                error('Error opening file to write. Check that directory exists.');     
            end
            
            fprintf(fid,'Filename: %s\n',dataIn.filename);
            fprintf(fid,'%s\n',datestr(now));
            
            %Print the headers
            fprintf(fid,',');
            for iC = 1:numChannels
                for ii = 1:(2*numAdditionalROIs) + 1
                    if ii == (numAdditionalROIs + 1)
                        fprintf(fid,'%s',channelNames{iC});
                    else
                        fprintf(fid,',');
                    end
                end
            end
            fprintf(fid,'\n');
            fprintf(fid,'Time (%s)',dataIn.timeUnit);
            for iC = 1:numChannels
                fprintf(fid,',ROI %d',-numAdditionalROIs:numAdditionalROIs);
            end
            fprintf(fid,'\n');
            
            %Print the data            
            for iR = 1:numRows
                
                fprintf(fid,'%g',dataIn.timestamps(iR));
                
                for iC = 1:numChannels
                    for ii = -numAdditionalROIs:numAdditionalROIs
                        
                        if ii < 0
                            fprintf(fid,',%g',dataIn.(channelNames{iC}).ROIAbove(abs(ii)).RawIntensity(iR));
                            
                        elseif ii == 0
                            fprintf(fid,',%g',dataIn.(channelNames{iC}).RawIntensity(iR));
                            
                        elseif ii > 0
                            fprintf(fid,',%g',dataIn.(channelNames{iC}).ROIBelow(abs(ii)).RawIntensity(iR));
                        end     
                    end
                end
                fprintf(fid,'\n');
            end
            
            fclose(fid);
            
            %---Export the normalized data values---
            %Rename the file
            filename = [filename(1:end-4),'_normalized.csv'];
            fid = fopen(filename,'w');
            if fid == -1
                error('Error opening file to write. Check that directory exists.');     
            end
            
            fprintf(fid,'Filename: %s\n',dataIn.filename);
            fprintf(fid,'%s\n',datestr(now));
            
            %Print the headers
            fprintf(fid,',');
            for iC = 1:numChannels
                for ii = 1:(2*numAdditionalROIs) + 1
                    if ii == (numAdditionalROIs + 1)
                        fprintf(fid,'%s',channelNames{iC});
                    else
                        fprintf(fid,',');
                    end
                end
            end
            fprintf(fid,'\n');
            fprintf(fid,'Time (%s)',dataIn.timeUnit);
            for iC = 1:numChannels
                fprintf(fid,',ROI %d',-numAdditionalROIs:numAdditionalROIs);
            end
            fprintf(fid,'\n');
            
            %Print the data            
            for iR = 1:numRows
                
                fprintf(fid,'%g',dataIn.timestamps(iR));
                
                for iC = 1:numChannels
                    for ii = -numAdditionalROIs:numAdditionalROIs
                        
                        if ii < 0
                            fprintf(fid,',%g',dataIn.(channelNames{iC}).ROIAbove(abs(ii)).NormIntensity(iR));
                            
                        elseif ii == 0
                            fprintf(fid,',%g',dataIn.(channelNames{iC}).NormIntensity(iR));
                            
                        elseif ii > 0
                            fprintf(fid,',%g',dataIn.(channelNames{iC}).ROIBelow(abs(ii)).NormIntensity(iR));
                        end     
                    end
                end
                fprintf(fid,'\n');
            end
            
            fclose(fid);
            
        end
        
        function imgOut = annotateImage(imgIn, roiImg, roiLabel)
           
            imgOut = DNARepairAnalyzer.showoverlay(DNARepairAnalyzer.normalizeimg(imgIn),bwperim(roiImg),[1 1 1]);
            
            %Find the middle of the ROI
            txtStrPos = [find(any(roiImg,1),1,'last'),...
                0.5*(find(any(roiImg,2),1,'last') - find(any(roiImg,2),1,'first')) + find(any(roiImg,2),1,'first')];
            imgOut = insertText(imgOut,txtStrPos,roiLabel);
            
        end
        
        function exportInfo(saveDir, fname, ROImask, nuclMask)
           %EXPORTINFO  Export additional data for Eric
           
           %From the ROI, calculate the [left top width heigth]
           rowStart = find(any(ROImask,2),1,'first');
           rowEnd = find(any(ROImask,2),1,'last');
           height = rowEnd - rowStart;
           
           colStart = find(any(ROImask,1),1,'first');
           colEnd = find(any(ROImask,1),1,'last');
           width = colEnd - colStart - 1;
           
           fid = fopen(fullfile(saveDir,[fname,'ROI.txt']),'w');
           fprintf(fid,'%d, %d, %d, %d', colStart, rowStart, width, height);
           fclose(fid);
           
           %From the nuclear mask, get the nuclear edge coordinate and
           %export it as a txt file
           [startI, startJ] = find(nuclMask,1,'first');
           contour = bwtraceboundary(nuclMask,[startI, startJ], 'NE', 8, Inf, 'clockwise');

           fid = fopen(fullfile(saveDir,[fname,'NuclMask.txt']),'w');
           fprintf(fid,'%d, %d\n', contour');
           fclose(fid);
           
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
        
        function varargout = showoverlay(baseimage, mask, color, varargin)
            %SHOWOVERLAY    Plot an overlay mask on an image
            %
            %  SHOWOVERLAY(IMAGE,MASK,COLOR) will plot an overlay specified by a binary
            %  MASK on the IMAGE. The color of the overlay is specified using a three
            %  element vector COLOR.
            %
            %  Example:
            %
            %    mainImg = imread('cameraman')
            %
            %
            %  Downloaded from http://cellmicroscopy.wordpress.com
            
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

    end
    
    methods (Hidden)
        
        function optStruct = props2struct(obj)
            %PROPS2STRUCT  Convert object properties into a struct
            
            props = properties(obj);
            
            for iP = 1:numel(props)
                optStruct.(props{iP}) = obj.(props{iP});
            end
            
        end
        
    end
    
end