function data = autofocus_cell(mmc,cellstats,cropFactor,configPhase,configFluoCFP,configFluoYFP,configFluoRFP,background)

mmc.setConfig('Channel',configPhase);
mmc.waitForSystem();

imWidth = mmc.getImageWidth();
imHeight = mmc.getImageHeight();
pixelSize = mmc.getPixelSizeUm();
zPosition = mmc.getPosition('ZStage');

if mmc.getBytesPerPixel == 2
    pixelType = 'uint16';
else
    pixelType = 'uint8';
end

data = struct('Area',[],'Eccentricity',[],'MajorAxisLength',[],'MinIntensity',[],'Phase',[],'YFP',[],'CFP',[],'RFP',[],'Zcoord',[],'XYcoord',[],'time',[],'Centroid',[]);

%% take z-stack pictures with 150X phase

nstacks = 7; %must be odd

imstack = zeros(2*round(imWidth*0.5*cropFactor),2*round(imHeight*0.5*cropFactor),nstacks,'double');

mmc.setAutoShutter(0);
mmc.setShutterOpen(1);

steps = 3;

for i=1:nstacks
    mmc.setPosition('ZStage',zPosition+(steps*(i-(nstacks+1)/2)));
    mmc.waitForImageSynchro();
    mmc.snapImage();
    imgtmp=mmc.getImage;
    imgtmp = typecast(imgtmp, pixelType);
    img=rotateFrame(imgtmp,mmc);
    imstack(:,:,i)= double(img((1+imWidth/2-round(imWidth*0.5*cropFactor)):(imWidth/2+round(imWidth*0.5*cropFactor)),(1+imHeight/2-round(imHeight*0.5*cropFactor)):(imHeight/2+round(imHeight*0.5*cropFactor)))) - background((1+imWidth/2-round(imWidth*0.5*cropFactor)):(imWidth/2+round(imWidth*0.5*cropFactor)),(1+imHeight/2-round(imHeight*0.5*cropFactor)):(imHeight/2+round(imHeight*0.5*cropFactor)));
end

imstack  = imstack - min(imstack(:));
imstack = imstack./max(imstack(:));

mmc.setShutterOpen(0);

%% detect dark object in all z

best = -1;
boxp = [1 1 2 2];
zp = zPosition;

%  figure,
for i=1:nstacks
    
    imnorm = imstack(:,:,i);

    meanint = mean(imnorm(:));
    bwim = ~im2bw(imnorm,meanint);

    imef = rangefilt(imnorm);

    imbw = im2bw(imef,0.1);
    imbw = imdilate(imbw, [1 1 1;1 1 1; 1 1 1]);
    imbw = imfill(imbw,'holes');

    
    imbw = imbw & bwim;

    imbw_filled = medfilt2(imbw,[5 5]);
    
    objects=regionprops(imbw_filled,'Area','Centroid','Eccentricity','BoundingBox','MajorAxisLength');
    objects(([objects.Area].*pixelSize^2)<0.5)=[];

    % pick most likely to be E. coli cell
    if numel(objects)>0
        cent = vertcat(objects.Centroid);
        dcent = ((size(imstack,1)/2-cent(:,1)')./(imWidth/100)).^2 + ((size(imstack,2)/2-cent(:,2)')./(imHeight/100)).^2;
        dsize = ((cellstats.Area - [objects.Area].*pixelSize^2)/cellstats.stdArea).^2;
        decc = ((cellstats.Eccentricity - [objects.Eccentricity])/cellstats.stdEccentricity).^2;
        daxis = ((cellstats.MajorAxisLength - [objects.MajorAxisLength]*pixelSize)/cellstats.stdMajorAxisLength).^2;
        idx = sub2ind(size(imnorm), round(cent(:,2)), round(cent(:,1)));
        dintensity = imnorm(idx)';
        scores = dsize + decc + daxis + dcent + dintensity;
        [s, ind] = min(scores);
    else
        continue;
    end
    
    if ~isempty(ind)

        box = [round(objects(ind).Centroid(1)-5/pixelSize) round(objects(ind).Centroid(2)-5/pixelSize) round(10/pixelSize) round(10/pixelSize)];
        object = imcrop(imnorm, box);
        glcmb = graycomatrix(object,'NumLevels',32);
        statsb = graycoprops(glcmb,'Contrast');
        if (statsb.Contrast/s)>best
            boxp = box;
            best=statsb.Contrast/s;
            zp = zPosition+(steps*(i-(nstacks+1)/2));
            data.Centroid = objects(ind).Centroid;
        end
        
    end
    
end


if best<0
    return;
end

%% center stage on selected object

xtarget = mmc.getXPosition(mmc.getXYStageDevice()) + round(pixelSize*(size(imstack,1)/2 - data.Centroid(1)));
ytarget = mmc.getYPosition(mmc.getXYStageDevice()) + round(pixelSize*(size(imstack,2)/2 - data.Centroid(2)));

mmc.setXYPosition('XYStage',xtarget,ytarget);
%%
boxp = [round(imWidth/2-5/pixelSize), round(imHeight/2-5/pixelSize), round(10/pixelSize), round(10/pixelSize)];

mmc.setPosition('ZStage',zp);
mmc.setShutterOpen(1);
mmc.waitForImageSynchro();

mmc.snapImage;
imgtmp = mmc.getImage;
imgtmp = typecast(imgtmp, pixelType);
img = double(rotateFrame(imgtmp,mmc));

imgp = imcrop(img, boxp);
mingray = min(imgp(:)) - (max(imgp(:))-min(imgp(:)))/2;
maxgray = max(imgp(:)) + (max(imgp(:))-min(imgp(:)))/2;

glcmb = graycomatrix(imgp,'GrayLimits',[mingray maxgray],'NumLevels',32);
statsb = graycoprops(glcmb,'Contrast');
bw = edge(imgp,'log');
bw = imfill(bw,'holes');
prop = regionprops(bw,'Area','Centroid');
[ma, mi] = max([prop.Area]);
mint = sign(-imgp(round(prop(mi).Centroid(2)),round(prop(mi).Centroid(1))) + mean(imgp(:)));
contrast = mint * statsb.Contrast;

if ma*pixelSize^2 > 1
    hasarea = 1;
    contarea = contrast;
else
    hasarea = 0;
    contarea = 0;
end

subplot(2,3,1), imshow(imgp,[],'InitialMagnification',100);
xlabel([hasarea statsb.Contrast]);
hold off;

drawnow;

steps = steps / 3;
ozp = zp;
areazp = zp;

while abs(steps) > 0.3
    j=0;
    mmc.setPosition('ZStage',zp+steps);
    mmc.waitForImageSynchro();
    mmc.snapImage;
    imgtmp = mmc.getImage;
    imgtmp = typecast(imgtmp, pixelType);
    img=double(rotateFrame(imgtmp,mmc));
    object = imcrop(img, boxp);
    glcmb = graycomatrix(object,'GrayLimits',[mingray maxgray],'NumLevels',32);
    statsb = graycoprops(glcmb,'Contrast');
    bw = edge(object,'log');
    bw = imfill(bw,'holes');
    
    prop = regionprops(bw,'Area','Centroid');
    [ma, mi] = max([prop.Area]);
    mint = sign(-object(round(prop(mi).Centroid(2)),round(prop(mi).Centroid(1))) + mean(object(:)));
    
    if (mint * statsb.Contrast) > (contrast + abs(contrast/50)) %&& ma*pixelSize^2 > 1
        contrast = mint * statsb.Contrast;
        ozp = zp;
        nzp = zp+steps;
        j=1;
        if ma*pixelSize^2 > 1
            areazp = nzp;
            hasarea = 1;
            contarea = contrast;
        end
    end
    
    if (zp-steps) ~= ozp
        steps = -steps;
        mmc.setPosition('ZStage',zp+steps);
        mmc.waitForImageSynchro();
        mmc.snapImage;
        imgtmp=mmc.getImage;
        imgtmp = typecast(imgtmp, pixelType);
        img=double(rotateFrame(imgtmp,mmc));
        object = imcrop(img, boxp);
        glcmb = graycomatrix(object,'GrayLimits',[mingray maxgray],'NumLevels',32);
        statsb = graycoprops(glcmb,'Contrast');
        bw = imfill(bw,'holes');
        prop = regionprops(bw,'Area','Centroid');
        [ma, mi] = max([prop.Area]);
        mint = sign(-object(round(prop(mi).Centroid(2)),round(prop(mi).Centroid(1))) + mean(object(:)));
        
        if (mint * statsb.Contrast)>(contrast + abs(contrast/50)) %&& ma*pixelSize^2 > 1
            contrast = mint * statsb.Contrast;
            ozp = zp;
            nzp = zp+steps;
            j=1;
            if ma*pixelSize^2 > 1
                areazp = nzp;
                hasarea = 1;
                contarea = contrast;
            end
        end
    end
    if j == 0
        steps = steps/3;
    else
        zp = nzp;
    end
end

if hasarea && contarea > (contrast - abs(contrast/50));
    zp = areazp;
end

mmc.setPosition('ZStage',zp);
mmc.waitForImageSynchro();
mmc.snapImage;
imgtmp=mmc.getImage;
imgtmp = typecast(imgtmp, pixelType);
img=double(rotateFrame(imgtmp,mmc));
imgp = imcrop(img, boxp);
glcmb = graycomatrix(imgp,'GrayLimits',[mingray maxgray],'NumLevels',32);
statsb = graycoprops(glcmb,'Contrast');

mmc.setShutterOpen(0);

data.Zcoord = zp;
data.Phase = imgp;



%% check if cell is good

bw = edge(imgp,'log');
bw = imfill(bw,'holes');

prop = regionprops(bw, imgp, 'Area','Eccentricity','Centroid','MajorAxisLength','MinIntensity');

[ma, mi] = max([prop.Area]);

data.Area = ma*pixelSize^2;
data.Eccentricity = prop(mi).Eccentricity;
data.MajorAxisLength = prop(mi).MajorAxisLength*pixelSize;
data.MinIntensity = prop(mi).MinIntensity;

data.XYcoord = [xtarget+pixelSize*(prop(mi).Centroid(1)-size(bw,1)/2); ytarget+pixelSize*(prop(mi).Centroid(2)-size(bw,2)/2)];

subplot(2,3,2), imshow(imgp,[],'InitialMagnification',100);
xlabel([contrast statsb.Contrast]);


subplot(2,3,3), imshow(bw,[],'InitialMagnification',100);
xlabel([data.Area data.MajorAxisLength]);


drawnow;

if data.Area > 1 && data.Area < 8 && data.MajorAxisLength > 1 && data.MajorAxisLength < 10 && data.Eccentricity > 0.8 && (data.MinIntensity - mean(imgp(:)))<0
    % take fluorescence image
    
    
    %% take CFP
    
    mmc.setConfig('Channel',configFluoCFP);
    mmc.waitForImageSynchro();
    
    mmc.snapImage;
    imgtmp=mmc.getImage;
    imgtmp = typecast(imgtmp, pixelType);
    imgf=imcrop(double(rotateFrame(imgtmp,mmc)), boxp);
    

    subplot(2,3,4), imshow(imgf,[],'InitialMagnification',100);
    xlabel('CFP');

    data.CFP = imgf;
    
    %% take YFP
    mmc.setConfig('Channel',configFluoYFP);
    mmc.waitForImageSynchro();
    
    mmc.snapImage;
    imgtmp=mmc.getImage;
    imgtmp = typecast(imgtmp, pixelType);
    
    imgf=imcrop(double(rotateFrame(imgtmp,mmc)), boxp);
    

    subplot(2,3,5), imshow(imgf,[],'InitialMagnification',100);
    xlabel('YFP');

    data.YFP = imgf;
    
    %% take RFP
    
    mmc.setConfig('Channel',configFluoRFP);
    mmc.waitForImageSynchro();
    
    mmc.snapImage;
    imgtmp=mmc.getImage;
    imgtmp = typecast(imgtmp, pixelType);
    
    imgf=imcrop(double(rotateFrame(imgtmp,mmc)), boxp);

    subplot(2,3,6), imshow(imgf,[],'InitialMagnification',100);
    xlabel('RFP');
    data.RFP = imgf;
    drawnow;
end

end