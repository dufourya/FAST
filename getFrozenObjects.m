function last_objects = getFrozenObjects(file, mmc,show_im)
%% read frame timestamp from file metadata
imgfile = strcat(file,'.bin');
metafile = strcat(file,'.meta.txt');

fid = fopen(metafile,'r');
metadata = textscan(fid,'%[^;];%[^;];%[^\n]');
fclose(fid);

for p = 1:size(metadata{:,2},1)
    if strcmp(metadata{1,2}(p,1),'NImagesBeforeFreeze')
        numImagesbeforeFreeze = str2double(metadata{1,3}(p,1));
    elseif strcmp(metadata{1,2}(p,1),'ImageWidth')
        imWidth = str2double(metadata{1,3}(p,1));
    elseif strcmp(metadata{1,2}(p,1),'ImageHeight')
        imHeight = str2double(metadata{1,3}(p,1));
    end
end

%% calculate image background to be substracted from images

background = zeros(imWidth,imHeight,'double');
j= floor(numImagesbeforeFreeze/100);

fid = fopen(imgfile,'r');
for i=1:100,
    fseek(fid,i*j*imWidth*imHeight*2,'bof');
    frame = rotateFrame(fread(fid,imWidth*imHeight,'*int16'),mmc);
    background = background + double(frame);
end
background = background./100;
meanbg = mean(reshape(background,1,[]));
fclose(fid);

fid = fopen(strcat('last_',file,'.bin'));
last_frame = fread(fid,'*int16');
fclose(fid);

last_frame = rotateFrame(last_frame,mmc);
last_objects = process_frame(last_frame,background, meanbg);

if show_im==1
    figure, imshow(double(last_frame)-background,[]);
    hold on
    scatter(last_objects(:,1),last_objects(:,2));
    hold off
%     figure;
%     imshow(background,[]);
end

return;