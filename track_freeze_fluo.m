file = 'R1500_I10_R1';
recordingtime = 300;
freezetime = 5;
cropFactor = 0.5;
bleachedzone = 30; %radius of bleached zone in um
fps = 0; %frame per second for movie
bit8 = 1; %convert images to 8 bit
crop10X = 0; % do not crop frame
totalTime = 60*60; %total time limit for entire experiment

availableConfigs = cell(mmc.getAvailableConfigs('Channel').toArray());
%%
configTracking = 'Phase_DAPI_10X';
configPhase = 'Phase_100X';
configFluoCFP = 'CFP_100X';
configFluoYFP = 'YFP_100X';
configFluoRFP = 'RFP_100X';

cellstatsFile = '';
backgroundFile = '';

if exist(cellstatsFile,'file')
    load(cellstatsFile);
else
    cellstats = struct('Area',3,'Eccentricity',1,'MajorAxisLength',2.5,'stdArea',2,'stdEccentricity',0.5,'stdMajorAxisLength',1);
end

%%
waitfor(Live(mmc,configTracking));
%%
tic
mmc.setOriginXY(mmc.getXYStageDevice());
sequenceAcquisitionFreeze(mmc, file, recordingtime, configTracking, freezetime, fps, bit8, crop10X);
%%
%process the last frame to get the coordinate of the frozen cells
objects = getFrozenObjects(file, mmc, 1);

imWidth = mmc.getImageWidth();
imHeight = mmc.getImageHeight();
pixelSize = mmc.getPixelSizeUm();

if exist(backgroundFile,'file')
    load(backgroundFile);
    if sum(size(background) == [imHeight,imWidth]) ~= 2
        background = zeros(imWidth,imHeight);
    end
else
    background = zeros(imWidth,imHeight);
end

%%
%configure the microscope for acquisition at 100X
q = questdlg('Switch objective to 100 X. Switch phase ring to Ph3. Turn off the light. Disengage ND filter for DiaLight. Close the epifluorescence lamp diaphragm.', 'Warning','OK','Cancel','OK');
%
%start acquisition of single cell fluorescence
if strcmp(q,'OK')
    
    mmc.setConfig('Channel',configPhase);
    %calculate z-score of object sizes extracted from the last frame at 10X
    z = (objects(:,3)-mean(objects(:,3)))./std(objects(:,3));
    %remove objects of unusual sizes
    objects(abs(z)>3,:)=[];
    
    %calculate the coordinates of each object from the center of the frame
    %to translate it to the motorizes stage coordinates
    coordinates = objects(:,[1 2]);
    [b, ix] = sort((coordinates(:,1)-imWidth/2).^2 + (coordinates(:,2)-imHeight/2).^2);
    
    %move the stage to the first object to help the user refocus the
    %microscope
    x = coordinates(ix(1),1);
    y = coordinates(ix(1),2);
    
    xtarget = round(pixelSize * (imWidth/2  - x));
    ytarget = round(pixelSize * (imHeight/2 - y));
    
    if strcmp(char(mmc.getCameraDevice()),'HamamatsuHam_DCAM')
        xcorrection = 10;
        ycorrection = -12;
    elseif strcmp(char(mmc.getCameraDevice()),'Andor')
        xcorrection = -53;
        ycorrection = -45;
    else
        xcorrection = 0;
        ycorrection = 0;
    end
    
    mmc.setXYPosition(mmc.getXYStageDevice(),xtarget-xcorrection,ytarget-ycorrection);
    
    %start live acquisition for the user to focus on the first object
    waitfor(Live(mmc,configPhase));
    
    xcorrection = xtarget - mmc.getXPosition(mmc.getXYStageDevice());
    ycorrection = ytarget - mmc.getYPosition(mmc.getXYStageDevice());
    
    quest = questdlg('Are you ready to continue?', 'Warning', 'Yes','No','No');
    
    %start automated cell acquisition
    if strcmp(quest,'Yes')
        
        %create cell array to store pictures
        pictures = cell(size(coordinates,1),1);
        
        %get z coordinates
        pos = mmc.getPosition('ZStage');
        newpos = [];
        
        mmc.setConfig('Channel',configPhase);
        mmc.waitForSystem();
        
        % write metadata for fluorescence pictures
        picmetafile = strcat(file,'.pictures.meta.txt');
        fid = fopen(picmetafile,'w');
        fprintf(fid,'%s;%s;%s\n', 'Experiment', 'Date', date);
        fprintf(fid,'%s;%s;%s\n', 'Experiment', 'Time', datestr(rem(now,1)));
        fprintf(fid,'%s;%s;%s\n', 'Experiment', 'DateVector', num2str(clock));
        fprintf(fid,'%s;%s;%s\n', 'Experiment', 'Computer', getenv('COMPUTERNAME'));
        fprintf(fid,'%s;%s;%s\n', 'Experiment', 'User', getenv('USERNAME'));
        metadata = mmc.getLoadedDevices();
        for i=1:metadata.size()
            device = metadata.get(i-1);
            deviceProp = mmc.getDevicePropertyNames(char(device));
            for j = 1:deviceProp.size()
                prop = deviceProp.get(j-1);
                propVal = mmc.getProperty(char(device),char(prop));
                fprintf(fid,'%s;%s;%s\n', char(device), char(prop), char(propVal));
            end
        end
        fprintf(fid,'Camera;CurrentPixelSize_um;%f\n', mmc.getPixelSizeUm());
        fprintf(fid,'Camera;ImageWidth;%f\n', mmc.getImageWidth());
        fprintf(fid,'Camera;ImageHeight;%f\n', mmc.getImageHeight());
        fprintf(fid,'ObjectiveOffset;Xcorrection;%f\n', xcorrection);
        fprintf(fid,'ObjectiveOffset;Ycorrection;%f\n', ycorrection);
        fclose(fid);
        
        currentcoord = coordinates(ix(1),:);
        j=1;
        f=0;
        
        close all;
        drawnow;
        
        while size(coordinates,1) > 0 && toc < totalTime
            
            close all;
            fprintf('%d objects captured, %d objects remaining\n',f,size(coordinates,1));
            
            %get the object closest to the current stage position to
            %minimize stage position
            [b, ix] = sort((coordinates(:,1)-currentcoord(1)).^2 + (coordinates(:,2)-currentcoord(2)).^2);
            currentcoord = coordinates(ix(1),:);
            
            xtarget = round(pixelSize * (imWidth/2  - currentcoord(1))) - xcorrection;
            ytarget = round(pixelSize * (imHeight/2 - currentcoord(2))) - ycorrection;
            
            mmc.setXYPosition(mmc.getXYStageDevice(),xtarget,ytarget);
            %waitfor(Live(mmc,configPhase));
                        
            %start automated focusing and picture acquisition
            pictures{j} = autofocus_cell(mmc,cellstats,cropFactor,configPhase,configFluoCFP,configFluoYFP,configFluoRFP,background);
            pictures{j}.time = toc;
            %delete current object from queue
            coordinates(ix(1),:) = [];
            
            if ~isempty(pictures{j}.Phase)
                %calculate the distance of the remaining objects from the
                %center of the bleached zone
                xtarget_all = round(pixelSize * (imWidth/2  - coordinates(:,1))) - xcorrection;
                ytarget_all = round(pixelSize * (imHeight/2 - coordinates(:,2))) - ycorrection;
                d = sqrt((xtarget_all-pictures{j}.XYcoord(1)).^2 + (ytarget_all-pictures{j}.XYcoord(2)).^2);
                %remove bleached objects from the queue
                coordinates(d<bleachedzone,:)=[];
                
                %add the object z coordinate to the list and calculate the
                %mean z coordinate as a starting pint for the next object
                newpos = [newpos pictures{j}.Zcoord];
                pos = mean(newpos);
                f=f+1;
            end
            j=j+1;
        end
        
        %remove objects with no picture taken and save data
        pictures(cellfun('isempty',pictures))=[];
        for i = 1: numel(pictures)
            if isempty(pictures{i}.Phase)
                pictures{i}.Centroid = [];
            end
        end
        pictures = horzcat(pictures{:});
        pictures(cellfun('isempty',{pictures.Phase}))=[];
        save(strcat(file,'.pictures.mat'),'pictures');
        
    end
    
    %%
    %calculate the statistics of the captured objects
    ind = [];
    meanepi = [];
    for i=1:numel(pictures)
        if ~isempty(pictures(i).CFP)
            ind = [ind i];
            meanepi = [meanepi mean(pictures(i).CFP(:))];
        end
    end
    
    cellstats.Area = mean([pictures(ind).Area]);
    cellstats.stdArea = std([pictures(ind).Area]);
    cellstats.MajorAxisLength = mean([pictures(ind).MajorAxisLength]);
    cellstats.stdMajorAxisLength = std([pictures(ind).MajorAxisLength]);
    cellstats.Eccentricity = mean([pictures(ind).Eccentricity]);
    cellstats.stdEccentricity = std([pictures(ind).Eccentricity]);
    
    save(strcat(file,'.cellstats.mat'),'cellstats');
    
end

close all;