function sequenceAcquisitionFreeze(mmc, file, runtime, config, freezeTime, fps, bit8, crop)
    
    if nargin < 6 || isempty(fps) || fps < 0
        fps = 0;
    end
    if nargin < 7 || isempty(bit8) || bit8 ~= 0
        bit8 = 1;
    end
    if nargin < 8 || isempty(crop) || crop ~= 1
        crop = 0;
    end
    
    %% microscope settings
    import mmcorej.*;
    f_name = getFunctionName();
    computer_name = getenv('COMPUTERNAME');
    camera = mmc.getCameraDevice();
    
    mmc.setConfig('System','Startup');
    metadata = mmc.getLoadedDevices();
    mmc.setConfig('Channel',config);
    mmc.waitForSystem();
    checkTemperature(mmc);
    
    if mmc.getBytesPerPixel == 2
        pixelType = 'uint16';
        realbit = 16;
    else
        pixelType = 'uint8';
        realbit = 8;
    end
    
    if strcmp(computer_name, 'EMONETSCOPE')
        circ_buff_size_mb = 5000;
    elseif strcmp(computer_name, 'ZENSCOPE')
        circ_buff_size_mb = 10000;
    else
        error(['[' f_name ']' ' Do not know memory requirements of computer ' ...
            '''' computer_name, '''.']);
    end
    
    if strcmp('HamamatsuHam_DCAM',camera)
        interval = str2double(mmc.getProperty(camera,'Exposure'));
    end
    %mmc.setAutoShutter(0);
    %% open files to save data
    imgfile = strcat(file, FileInfo.binExt);
    metafile = strcat(file, FileInfo.metaExt);
    
    fid = fopen(imgfile,'w');
    fid1 = fopen(metafile,'w');
    
    fprintf(fid1,'%s;%s;%s\n', 'Experiment', 'Date', date);
    fprintf(fid1,'%s;%s;%s\n', 'Experiment', 'Time', datestr(rem(now,1)));
    fprintf(fid1,'%s;%s;%s\n', 'Experiment', 'DateVector', num2str(clock));
    fprintf(fid1,'%s;%s;%s\n', 'Experiment', 'Computer', computer_name);
    fprintf(fid1,'%s;%s;%s\n', 'Experiment', 'User', getenv('USERNAME'));
    %% start sequence acquisition
    
    mmc.setCircularBufferMemoryFootprint(circ_buff_size_mb);
    mmc.initializeCircularBuffer();
    
    mmc.prepareSequenceAcquisition(camera)
    
    for i=1:metadata.size()
        device = metadata.get(i-1);
        deviceProp = mmc.getDevicePropertyNames(char(device));
        for j = 1:deviceProp.size()
            prop = deviceProp.get(j-1);
            propVal = mmc.getProperty(char(device),char(prop));
            fprintf(fid1,'%s;%s;%s\n', char(device), char(prop), char(propVal));
        end
    end
    
    if strcmp('HamamatsuHam_DCAM',char(camera))
        mmc.setProperty(camera,'TRIGGER SOURCE','INTERNAL');
    end
    
    w=mmc.getImageWidth();
    h=mmc.getImageHeight();
    
    mmc.waitForImageSynchro();
    mmc.startContinuousSequenceAcquisition(0);
    
    % HAVE to set interval AFTER starting sequence acquisition, since the camera
    % resets it when the method is called.
    if strcmp('Andor',char(camera))
        interval = str2double(char(mmc.getProperty('Andor','ActualInterval-ms')));
        realbit = 14;
    end
    
    nImages = ceil(runtime*1000/interval);
    skip_frame = 1;
    
    if fps ~=0
        skip_frame = floor(1000/(fps * interval));
        if skip_frame < 1
            mmc.stopSequenceAcquisition();
            mmc.setShutterOpen(0);
            mmc.waitForSystem();
            fclose(fid);
            fclose(fid1);
            error('The target aquisition frame rate cannot be achieved with this configuration! Use fps = 0 if you want to record at the maximum possible frame rate.');
        end
    end
    
    if bit8 == 1
        fprintf('Converting images to 8 bits\n');
    else
        fprintf('Recording images in 16 bits\n');
    end
    if crop == 1
        fprintf('Frame is cropped by 1/2\n');
    else
        fprintf('Recording full frame\n');
    end
    fprintf('Skipping %d frame(s)\n', skip_frame-1);
    fprintf('Image interval = %.2f ms\n', skip_frame*interval);
    fprintf('Recording at %.2f frames/sec\n', 1000/(skip_frame*interval));
    
    i=0;
    tic1=tic;
    nImages_actual = 0;
    
    
    while i < nImages
        if mmc.getRemainingImageCount()>0
            img = mmc.popNextImage();
            if mod(i, skip_frame) == 0
                img = typecast(img, pixelType);
                if crop == 1
                    img = reshape(img,w,h);
                    img = img((w/2-w/4+1):(w/2+w/4), (h/2-h/4+1):(h/2+h/4));
                    img = img(:);
                end
                if bit8 == 1 && realbit ~= 8
                    img = uint8(double(img) * (2^8 - 1) / (2^realbit - 1));
                    fwrite(fid, img, 'uint8');
                else
                    fwrite(fid, img, pixelType);
                end
                nImages_actual = nImages_actual + 1;
            end
            i = i+1;
            if mod(i,ceil(60000/interval)) == 0
                toc(tic1);
            end
        end
    end
    toc(tic1);
    
    imagesBeforeFreeze = nImages_actual;
    
    %% freeze cells in place with UV
    tic;
    if strcmp(mmc.getShutterDevice(),'TIDiaShutter')
        mmc.setProperty('TIEpiShutter','State', '1');
    elseif strcmp(mmc.getShutterDevice(),'Arduino-Shutter')
        mmc.setProperty('Spectra','White_Enable',1);
        mmc.setProperty('Spectra','White_Level',100);
        mmc.setProperty('Spectra','State',1);
    end
    
    while toc < freezeTime,
        if mmc.getRemainingImageCount()>0
            img = mmc.popNextImage();
            if mod(i, skip_frame) == 0
                img = typecast(img, pixelType);
                if crop == 1
                    img = reshape(img,w,h);
                    img = img((w/2-w/4+1):(w/2+w/4), (h/2-h/4+1):(h/2+h/4));
                    img = img(:);
                end
                if bit8 == 1 && realbit ~= 8
                    img = uint8(double(img) * (2^8 - 1) / (2^realbit - 1));
                    fwrite(fid, img, 'uint8');
                else
                    fwrite(fid, img, pixelType);
                end
                nImages_actual = nImages_actual + 1;
            end
            i = i+1;
        end
    end
    
    
    if strcmp(mmc.getShutterDevice(),'TIDiaShutter')
        mmc.setProperty('TIEpiShutter','State', '0');
    elseif strcmp(mmc.getShutterDevice(),'Arduino-Shutter')
        mmc.setProperty('Spectra','State',0);
    end
    
    tic;
    while toc < 5,
        if mmc.getRemainingImageCount()>0
            img = mmc.popNextImage();
            if mod(i, skip_frame) == 0
                img = typecast(img, pixelType);
                if crop == 1
                    img = reshape(img,w,h);
                    img = img((w/2-w/4+1):(w/2+w/4), (h/2-h/4+1):(h/2+h/4));
                    img = img(:);
                end
                if bit8 == 1 && realbit ~= 8
                    img = uint8(double(img) * (2^8 - 1) / (2^realbit - 1));
                    fwrite(fid, img, 'uint8');
                else
                    fwrite(fid, img, pixelType);
                end
                nImages_actual = nImages_actual + 1;
            end
            i = i+1;
        end
    end
    mmc.stopSequenceAcquisition();
    toc(tic1);
    mmc.waitForImageSynchro();
    
    while mmc.getRemainingImageCount()>0
        img = mmc.popNextImage();
        if mod(i, skip_frame) == 0
            img = typecast(img, pixelType);
            if crop == 1
                img = reshape(img,w,h);
                img = img((w/2-w/4+1):(w/2+w/4), (h/2-h/4+1):(h/2+h/4));
                img = img(:);
            end
            if bit8 == 1 && realbit ~= 8
                img = uint8(double(img) * (2^8 - 1) / (2^realbit - 1));
                fwrite(fid, img, 'uint8');
            else
                fwrite(fid, img, pixelType);
            end
            nImages_actual = nImages_actual + 1;
        end
        i = i+1;
    end
    
    mmc.setShutterOpen(0);
    fprintf('Images recorded: %d\n', nImages_actual);
    
    %% close files
    fprintf(fid1,'SequenceAcquisition;ActualIntervalBurst-ms;%f\n', skip_frame * interval);
    fprintf(fid1,'SequenceAcquisition;NumberImages;%d\n', nImages_actual);
    fprintf(fid1,'SequenceAcquisition;NImagesBeforeFreeze;%d\n', imagesBeforeFreeze);
    
    fprintf(fid1,'XYStage;XPosition_um;%d\n', mmc.getXPosition('XYStage'));
    fprintf(fid1,'XYStage;YPosition_um;%d\n', mmc.getYPosition('XYStage'));
    fprintf(fid1,'ZStage;ZPosition_um;%d\n', mmc.getPosition('ZStage'));
    
    fprintf(fid1,'Camera;CurrentPixelSize_um;%f\n', mmc.getPixelSizeUm());
    if bit8 == 1
        fprintf(fid1,'SequenceAcquisition;FilePixelType;*%s\n', 'uint8');
    else
        fprintf(fid1,'SequenceAcquisition;FilePixelType;*%s\n', pixelType);
    end
    if crop == 0
        fprintf(fid1,'Camera;ImageWidth;%f\n', mmc.getImageWidth());
        fprintf(fid1,'Camera;ImageHeight;%f\n', mmc.getImageHeight());
    else
        fprintf(fid1,'Camera;ImageWidth;%f\n', mmc.getImageWidth()/2);
        fprintf(fid1,'Camera;ImageHeight;%f\n', mmc.getImageHeight()/2);
    end
      
    fclose(fid1);
    fclose(fid);    
    
    fid2 = fopen(strcat('last_',imgfile), 'w');
    if bit8 == 1 && realbit ~= 8
        fwrite(fid, img, 'uint8');
    else
        fwrite(fid, img, pixelType);
    end
    fclose(fid2);
end