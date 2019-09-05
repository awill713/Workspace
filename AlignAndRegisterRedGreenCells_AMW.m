
clear
%% Load processedData file
[file folder] = uigetfile('/Users/Aaron/Documents/');
if file==0
    return;
end
load([folder '/' file]);


xOffset = spatialInfo.xrange(1)-1;
yOffset = spatialInfo.yrange(1)-1;

tempImage = zeros(1,size(spatialInfo.im,1)*size(spatialInfo.im,2));
roiLocations = cell(1,length(spatialInfo.ROIs));
for n = 1:length(spatialInfo.ROIs)
    tempImage(spatialInfo.ipix{1,n})=1;
    
    roiRow = mod(spatialInfo.ipix{1,n},size(spatialInfo.im,1));
    roiRow(find(roiRow==0)) = size(spatialInfo.im,1);
    roiColumn = floor(spatialInfo.ipix{1,n}/size(spatialInfo.im,1))+1;

    roiRow = roiRow + yOffset;
    roiColumn = roiColumn + xOffset;
    roiLocations{1,n} = [roiRow roiColumn];
    
end
tempImage = reshape(tempImage,size(spatialInfo.im,1),size(spatialInfo.im,2));

image = zeros(512,512);
image(spatialInfo.yrange,spatialInfo.xrange) = tempImage;
figure;imagesc(image);
title('Final image');

testImage = zeros(512,512);
% for n = 1:length(spatialInfo.ROIs)
for n = 1:length(roiLocations)
    for p = 1:length(roiLocations{1,n})
        testImage(roiLocations{1,n}(p,1),roiLocations{1,n}(p,2)) = 1;
    end
end
figure;imagesc(testImage);
title('ROI-wise');

%% Choose folder with green and red videos
recFolder = uigetdir('/Users/Aaron/Documents/');
greenFile = dir([recFolder '/*GREEN1.tif']);
redFile = dir([recFolder '/*RED1.tif']);
% greenFile(contains({greenFile.name},'._'))=[];
% redFile(contains({redFile.name},'._'))=[];

totalNeurons = size(calcium.rawTraces,1);
grTraces = zeros(totalNeurons,2,2000);
for f = 1:2000
    f
    gFrame = imread([recFolder '/' greenFile.name],f);
    rFrame = imread([recFolder '/' redFile.name],f);
    for n = 1:size(calcium.rawTraces,1)
        gPixels = zeros(1,length(roiLocations{1,n}));
        rPixels = zeros(1,length(roiLocations{1,n}));
        
        for p = 1:length(roiLocations{1,n})
            gPixels(p) = gFrame(roiLocations{1,n}(p,1),roiLocations{1,n}(p,2));
            rPixels(p) = rFrame(roiLocations{1,n}(p,1),roiLocations{1,n}(p,2));
        end
        
        grTraces(n,1,f) = mean(gPixels);
        grTraces(n,2,f) = mean(rPixels);
    end
end

grSlopeIntercept = zeros(totalNeurons,2);
for n = 1:totalNeurons
    p = polyfit(grTraces(n,1,:),grTraces(n,2,:),1);
    grSlopeIntercept(n,1) = p(1);
    grSlopeIntercept(n,2) = p(2);
end
ff1 = figure;histogram(grSlopeIntercept(:,1));
xlabel('Slope (dred / dgreen)');
ylabel('Number of neurons');
title('Neuron-wise slope');
ff2 = figure;histogram(grSlopeIntercept(:,2));
xlabel('Baseline red intensity');
ylabel('Number of neurons');
title('Neuron-wise red intensity');
saveas(ff1, fullfile(recFolder,['Neuron-wise red green slope.fig']));
saveas(ff2, fullfile(recFolder,['Neuron-wise red green y-intercept.fig']));

grMean = zeros(totalNeurons,3);
for n = 1:totalNeurons
    gMean = mean(grTraces(n,1,:));
    rMean = mean(grTraces(n,2,:));
    grMean(n,1:2) = [gMean rMean];
    grMean(n,3) = rMean / gMean;
end
figure;histogram(grMean(:,3))

figure;scatter(grSlopeIntercept(:,2),grMean(:,3));
%     f
%     greenVid(:,:,f) = imread([recFolder '/' greenVidFile.name],f);
%     redVid(:,:,f) = imread([recFolder '/' redVidFile.name],f);
% end

%% Choose folder with green and red videos
recFolder = uigetdir('/Users/Aaron/Documents/');

greenAvgFile = dir([recFolder '/*GREENavg.tif']);
greenAvgFile(contains({greenAvgFile.name},'._'))=[];
greenAvg = imread([recFolder '/' greenAvgFile.name]);

redAvgFile = dir([recFolder '/*REDavg.tif']);
redAvgFile(contains({redAvgFile.name},'._'))=[];
redAvg = imread([recFolder '/' redAvgFile.name]);

totalNeurons = length(calcium.n);
fluorIntensity = zeros(totalNeurons,4); %green, red, r:g ratio, z-scored ratio
for n = 1:totalNeurons
    pixels = length(roiLocations{1,n});
    gPixels = zeros(1,pixels);
    rPixels = zeros(1,pixels);
    for p = 1:pixels
    	gPixels(p) = greenAvg(roiLocations{1,n}(p,1),roiLocations{1,n}(p,2));
        rPixels(p) = redAvg(roiLocations{1,n}(p,1),roiLocations{1,n}(p,2));
    end
    fluorIntensity(n,1) = mean(gPixels);
    fluorIntensity(n,2) = mean(rPixels);
    fluorIntensity(n,3) = mean(rPixels)/mean(gPixels);
end
fluorIntensity(:,4) = zscore(fluorIntensity(:,3));
f1 = figure;scatter(fluorIntensity(:,1),fluorIntensity(:,2));
f2 = figure;histogram(fluorIntensity(:,3));
f3 = figure;histogram(fluorIntensity(:,4));

saveas(f1, fullfile(recFolder,['Neuron-wise red green scatter plot.fig']))
saveas(f2, fullfile(recFolder,['Neuron-wise red green ratio histogram.fig']))
saveas(f3, fullfile(recFolder,['Neuron-wise red green z-score histogram.fig']))
rawZ = zscore(fluorIntensity(:,2));
figure;histogram(rawZ)
figure;histogram(zscore(fluorIntensity(:,2)) - zscore(fluorIntensity(:,1)))
figure;histogram(fluorIntensity(:,2) - fluorIntensity(:,1));

thresholdPercentage = 0.1;
supraThresholdCount = round(thresholdPercentage * totalNeurons);
redNeurons = zeros(1,supraThresholdCount);
tempFluor = fluorIntensity(:,3);
for k = 1:supraThresholdCount
    [value index] = max(tempFluor);
    redNeurons(k) = index;
    tempFluor(index) = [];
end
calcium.redNeurons = redNeurons;
save([folder file],'calcium','-append');

% 
% totalFrames = length(calcium.rawTraces);
% greenVid = zeros(512,512,totalFrames);
% greenVidFile = dir([recFolder '/*GREENvid.tif']);
% greenVidFile(contains({greenVidFile.name},'._'))=[];
% redVid = zeros(512,512,totalFrames);
% redVidFile = dir([recFolder '/*REDvid.tif']);
% redVidFile(contains({redVidFile.name},'._'))=[];
% for f = 2940:totalFrames
%     f
%     greenVid(:,:,f) = imread([recFolder '/' greenVidFile.name],f);
%     redVid(:,:,f) = imread([recFolder '/' redVidFile.name],f);
% end

% 
% redFluor = zeros(totalNeurons,totalFrames);
% greenFluor = zeros(totalNeurons,totalFrames);
% for frame = 1:totalFrames
%     frame
%     green = imread([recFolder '/' greenVidFile.name],frame);
%     red = imread([recFolder '/' redVidFile.name],frame);
%     
%     for n = 1:totalNeurons
%         for p = 1:length(roiLocations{1,n})
%             greenFluor(n,frame) = size(green([roiLocations{1,n}(:,1),roiLocations{1,n}(:,2)]));
%             redFluor(n,frame) = mean(red(roiLocations{1,n}(:,1),roiLocations{1,n}(:,2)));
%     end
% end
