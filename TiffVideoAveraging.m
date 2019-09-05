
%% Load tiff videos
redfile_path = '/Volumes/AARON FILES/Two photon/2PM 004/AW055/2P red 20190331/Silence/Tiff videos/20190331AW055_tifStacks_AW055_20_plane1_RED.tif';
redVid = zeros(512,512,6000);
greenfile_path = '/Volumes/AARON FILES/Two photon/2PM 004/AW055/2P red 20190331/Silence/Tiff videos/20190331AW055_tifStacks_3_AW055_2P_plane1_GREEN.tif';
greenVid = zeros(512,512,6000);
for frame = 1:6000
    frame
    redVid(:,:,frame) = imread(redfile_path,frame);
    greenVid(:,:,frame) = imread(greenfile_path,frame);
end

%% Compare average pixel green-red intensity
redAvg = mean(redVid,3);
greenAvg = mean(greenVid,3);

figure;
scatter(reshape(greenAvg,[1 512*512]),reshape(redAvg,[1 512*512]));
xlabel('Avg green intensity')
ylabel('Avg red intensity')

%% Compare time-correlation of green-red signals
stats = zeros(1000,2);
xList = randi([1 512],1,1000);
yList = randi([1 512],1,1000);
for i = 1:1000
    i
    x = xList(i);
    y = yList(i);
    [rho pval] = corr(squeeze(greenVid(x,y,:)),squeeze(redVid(x,y,:)));
    stats(i,1) = rho;
    stats(i,2) = pval;
end
figure;
histogram((stats(:,1)));