function [tracks, good_pics, leftover_pics] = analyze_freezing_assay(file, do_plots)
%% load data file

fileInfo = FileInfo(file);

load(strcat(fileInfo.name, RadFileInfo.swimtrackerExt),'tracks');
metadataTracks = Metadata(fileInfo.meta_file);
metadataTracks.read();
metadataPics = Metadata(strcat(fileInfo.name, FileInfo.picsMetaExt));
metadataPics.read();
load(strcat(fileInfo.name, FileInfo.picturesExt),'pictures');

pixelSizeTracks = metadataTracks.getPixelSize();
pixelSizePics = metadataPics.getPixelSize();

imageHeight = metadataTracks.getImageHeightPixel();
imageWidth = metadataTracks.getImageWidthPixel();

xcorrection = metadataPics.getXCorrection();
ycorrection = metadataPics.getYCorrection();
%%
tracks = likelihood_tumblebias(tracks, do_plots);
%% find trajectories that finish at the freezing time
end_times = cellfun(@max,{tracks.time},'UniformOutput', false);
t = [end_times{:}] == max([end_times{:}]);
tracks = tracks(t);
traj_coordinates = zeros(numel(tracks),2);

for i = 1:numel(tracks)
    traj_coordinates(i,:) = [tracks(i).x(end), tracks(i).y(end)];
end
%
coordinates = cat(2,pictures.XYcoord);
coordinates = coordinates';
coordinates(:,1) = (imageWidth/2*pixelSizeTracks) - coordinates(:,1) - xcorrection;
coordinates(:,2) = (imageHeight/2*pixelSizeTracks) - coordinates(:,2) - ycorrection;

%% match frozen objects with trajectories
dx = repmat(coordinates(:,1),1,size(traj_coordinates,1))-repmat(traj_coordinates(:,1),1,size(coordinates,1))';
dy = repmat(coordinates(:,2),1,size(traj_coordinates,1))-repmat(traj_coordinates(:,2),1,size(coordinates,1))';

D = dx.^2 + dy.^2;

[minx,  idx] = min(D,[],1);
[miny,  idy] = min(D,[],2);

ex = find(idx(idy') == 1:length(idy));
ey = find(idy(idx)' == 1:length(idx));

%mean(sqrt([miny(ex)' minx(ey)]))
%%

if do_plots
    figure,
    scatter(coordinates(:,1), coordinates(:,2),'o');
    axis equal
    hold on
    scatter(traj_coordinates(:,1), traj_coordinates(:,2),'+');
    hold off
    
    figure,
    hist(sqrt([miny(ex)' minx(ey)]),50)
    
    figure,
    hold on
    scatter(traj_coordinates(idy(ex),1),traj_coordinates(idy(ex),2));
    scatter(coordinates(idx(ey),1),coordinates(idx(ey),2),'+','red');
    hold off
    
    figure,
    imshow(zeros(ceil(imageWidth*pixelSizeTracks),ceil(imageHeight*pixelSizeTracks)));
    hold on;
    r = rand(length(ey),3);
    
    for i = 1:length(ey)
        plot(tracks(ey(i)).x,tracks(ey(i)).y,'color',r(i,:));
        scatter(coordinates(idx(ey(i)),1),coordinates(idx(ey(i)),2),'MarkerFaceColor',r(i,:),'MarkerEdgeColor',r(i,:));
    end
    hold off;
    figure, fi = 0;
end


%% analyze object fluorescence levels and sizes from 100x pictures

for i = 1:size(pictures,2)
    
    bw = edge(pictures(i).Phase,'log');
    bw = imfill(bw,'holes');
    ang = regionprops(bw, 'Orientation','Eccentricity','MajorAxisLength','Area','Centroid');
    [~, mi] = max([ang.Area]);
    bw = bwselect(bw,ang(mi).Centroid(1),ang(mi).Centroid(2));
    r = imrotate(bw,-ang(mi).Orientation);
    volume = sum((sum(r,1)./2.*pixelSizePics).^2.*pi.*pixelSizePics);
    bw1 = imdilate(bw,strel('disk',5,4));
    axislength = ang(mi).MajorAxisLength * pixelSizePics;
    
    if isempty(pictures(i).CFP) || volume == 0
        YFP = NaN;
        CFP = NaN;
        RFP = NaN;
        bgYFP = NaN;
        bgCFP = NaN;
        bgRFP = NaN;
        volume = NaN;
        axislength = NaN;
    else
        bgYFP = mean(pictures(i).YFP(~bw1));
        bgCFP = mean(pictures(i).CFP(~bw1));
        bgRFP = mean(pictures(i).RFP(~bw1));
        YFP = (sum(pictures(i).YFP(bw1)) - (sum(sum(bw1)) * bgYFP));
        CFP = (sum(pictures(i).CFP(bw1)) - (sum(sum(bw1)) * bgCFP));
        RFP = (sum(pictures(i).RFP(bw1)) - (sum(sum(bw1)) * bgRFP));
        
        if do_plots
            boundary1 = bwboundaries(bw1);
            boundary = bwboundaries(bw);
            if numel(boundary) > 0
                subplot(6,8,mod(fi,48)+1),
                fi = fi +1;
                imshow(pictures(i).CFP',[]);
                b = boundary{1};
                b1 = boundary1{1};
                hold on
                plot(b(:,1),b(:,2),'color','green')
                plot(b1(:,1),b1(:,2),'color','red')
                hold off
            end
            title(volume);
        end
    end
    
    if sum(idx(ey)==i)>0
        tracks(idy(i)).sumYFP = YFP;
        tracks(idy(i)).sumRFP = RFP;
        tracks(idy(i)).sumCFP = CFP;
        tracks(idy(i)).bgYFP = bgYFP;
        tracks(idy(i)).bgRFP = bgRFP;
        tracks(idy(i)).bgCFP = bgCFP;        
        tracks(idy(i)).axislength = axislength;
        tracks(idy(i)).Phase = pictures(i).Phase;
        tracks(idy(i)).volume = volume;
        tracks(idy(i)).pictime = pictures(i).time;
        tracks(idy(i)).match_dist = sqrt(miny(i));
    end
    
end
        
leftover_pics = pictures(setdiff(1:numel(pictures),idx(ey)));
good_pics =  pictures(idx(ey));
%%

if do_plots
    temp = tracks(~cellfun('isempty',{tracks.sumYFP}));
    tempfluo = temp(~isnan([temp.sumYFP]));
    tempfluo = tempfluo([tempfluo.sumYFP]./[tempfluo.sumRFP]>0);
    tempfluo = tempfluo([tempfluo.trajtime]>30);
    fratio = log2([tempfluo.sumYFP]./[tempfluo.sumRFP]);
    fratio = (fratio - min(fratio))/max(fratio - min(fratio));
    
    figure,
    imshow(zeros(ceil(imageWidth*pixelSizeTracks),ceil(imageHeight*pixelSizeTracks)));
    hold on;
    for i = 1:numel(fratio)
        plot(tempfluo(i).x,tempfluo(i).y,'color',[1-fratio(i) fratio(i) 0],'linewidth',2);
        scatter(tempfluo(i).x(end),tempfluo(i).y(end),'MarkerFaceColor',[1-fratio(i) fratio(i) 0],'MarkerEdgeColor',[1-fratio(i) fratio(i) 0]);
    end
    hold off;
    
    figure,
    scatter(log10([tempfluo.fit_runtime]) , fratio(~cellfun('isempty',{tempfluo.fit_runtime})));
    xlabel('Log10 run length (s)');
    ylabel('Log2 YFP/RFP');
end
end