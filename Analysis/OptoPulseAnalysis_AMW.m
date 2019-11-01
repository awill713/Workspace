
clear;

experiment = 'EP001';
mouseID = 'AW093';
session = 'Session1';
date = '20191029';
% stimPath = 'D:\Code\Aaron\StimInfo\20190623_optoLaserPulse_1and5ms_400k_005_stimInfo';
stimPath = 'C:\Users\Aaron\Documents\MATLAB\Workspace\Analysis\EphyStimInfo\20190623_optoLaserPulse_1and5ms_400k_005_stimInfo';
analysisWindow = [-50 200]; %ms relative to stimulus onset
quantWindow = [5 15]; %ms after stimulus onset
frBinWidth = 10; %ms
baselineWindow = [-10 0];
rasterColors = [0 0.447 0.741;... %blue
    0.850 0.325 0.098];    %orange

% dataFolder = fullfile('D:\KiloSort\',mouseID,session,date,'SpikeMat');
% dataFolder = 'D:\Electrophysiology\EP001\AW091\20191021';
dataFolder = fullfile('D:\Electrophysiology\',experiment,mouseID,date,'SpikeMat');
dataFiles = dir(fullfile(dataFolder,'*optoLaserPulse*'));

% newDir = fullfile('D:\KiloSort\',mouseID,session,date,'OptoLaserResponses');
newDir = fullfile('D:\Electrophysiology\',experiment,mouseID,date,'OptoLaserResponses');
figureDir = fullfile(newDir,'Figures');
if ~exist(newDir)
    mkdir(newDir);
end
if ~exist(figureDir)
    mkdir(figureDir);
end

load(stimPath);
optoPulseDuration = stimInfo.optoPulseDuration;

totalUnits = length(dataFiles);

analysisEdges = analysisWindow*1000;
quantificationEdges = quantWindow*1000; %temporal resolution of cheetah acquisition is microseconds
binCount = (analysisWindow(2)-analysisWindow(1))/frBinWidth; %NOT a sliding window, should adjust to be sliding in the future
binScalar = 1000/frBinWidth;
binEdges = analysisWindow(1):frBinWidth:analysisWindow(2);
frScalar = 1000/(quantWindow(2)-quantWindow(1));

optoResponsiveUnits = [];
optoResponsiveSingleUnits = [];
optoResponsiveMultiUnits = [];
singleUnits = 0;
multiUnits = 0;
latencies = [];

for n = 1:totalUnits
    nData = load(fullfile(dataFiles(n).folder,dataFiles(n).name));
    eventTimes = nData.Events;
    spikeTimes = nData.SpikeData(1,:);
    
    optoRaster = cell(length(optoPulseDuration),3); %spike times (1), trial number (y-axis) (2), color (3)
    optoResponse = zeros(length(optoPulseDuration),5); %duration (1),  mean baseline (2), baseline std error (3), mean opto resp (4), opto std error (5),
    optoFR = zeros(length(optoPulseDuration),binCount);
    
    tempY = 0;
    peakBin = zeros(0,2);
    for d = 1:length(optoPulseDuration)
        duration = optoPulseDuration(d);
        durEvents = find(stimInfo.order==d);
        color = rasterColors(d,:);
        
        frTrain = zeros(length(durEvents),binCount);
        rasterX = [];
        rasterY = [];
        
        prePost = zeros(length(durEvents),2);
        for e = 1:length(durEvents)
            tempY = tempY+1;
            
            time = eventTimes(durEvents(e));
            alignedSpikes = spikeTimes-time;
            trialSpikes = alignedSpikes(find(alignedSpikes>analysisEdges(1) & alignedSpikes<analysisEdges(end)));
            rasterX = [rasterX trialSpikes/1000];
            rasterY = [rasterY repmat(tempY,[1 length(trialSpikes)])];
            
            frTrain(e,:) = histcounts(trialSpikes/1000,binEdges);
            prePost(e,1) = histcounts(trialSpikes/1000,[baselineWindow(1) baselineWindow(2)]);
%             prePost(e,2) = histcounts(trialSpikes/1000,[quantWindow(1) quantWindow(2)]);
%             prePost(e,2) = max(frTrain(e,(analysisWindow(1)/frBinWidth)*sign(analysisWindow(1))+1:end));
        end
        
        meanFrTrain = mean(frTrain,1);
        optoFR(d,:) = meanFrTrain*frScalar;
        
        peakBin(1,d) = min(find(meanFrTrain==max(meanFrTrain))); %min just in case there are multiple bins where it reaches an equal peak, choosing the earlier one
        prePost(:,2) = frTrain(:,peakBin(1,d));
        
        optoRaster{d,1} = rasterX;
        optoRaster{d,2} = rasterY;
        optoRaster{d,3} = color;
        optoResponse(d,1) = duration;
        optoResponse(d,2) = mean(prePost(:,1))*frScalar;
        optoResponse(d,3) = std(prePost(:,1))*frScalar/sqrt(length(prePost(:,1)));
        optoResponse(d,4) = mean(prePost(:,2))*frScalar;
        optoResponse(d,5) = std(prePost(:,2))*frScalar/sqrt(length(prePost(:,2)));
        
        
        [hh pp] = ttest2(prePost(:,1),prePost(:,2),'alpha',0.05/totalUnits/2);
        h(d) = hh; h(isnan(h))=0;
        p(d) = pp; p(isnan(p))=1;
    end
    
    unitData(n).raster = optoRaster;
    unitData(n).optoResponse = optoResponse;
    unitData(n).optoFR = optoFR;
    unitData(n).type = nData.CellInfo(6);
    if h(1) || h(2)
        optoResponsiveUnits = [optoResponsiveUnits nData.CellInfo(4)];
        if nData.CellInfo(6)==1
            optoResponsiveSingleUnits = [optoResponsiveSingleUnits nData.CellInfo(4)];
        elseif nData.CellInfo(6)==2
            optoResponsiveMultiUnits = [optoResponsiveMultiUnits nData.CellInfo(4)];
        end
        sigPeak = peakBin(h==1); %in the case of responsive to one duration and not the other, only take the responsive one
        lat = (min(sigPeak)-1)*frBinWidth + analysisWindow(1);
        latencies = [latencies; nData.CellInfo(4) lat];
    end
    if nData.CellInfo(6)==1
        singleUnits = singleUnits+1;
    elseif nData.CellInfo(6)==2
        multiUnits = multiUnits+1;
    end
    
    f1 = figure;
    subplot(1,3,1);
    for r = 1:size(optoRaster,1)
        scatter(optoRaster{r,1},optoRaster{r,2},[],optoRaster{r,3},'.');
        hold on;
    end
    xlim([analysisWindow(1) analysisWindow(end)]);
    xlabel('Time (ms)');
    ylabel('Trial type');
    
    subplot(1,3,2); hold on;
    for d = 1:size(optoFR,1)
        plot(binEdges(2:end)-frBinWidth/2,optoFR(d,:),'Color',optoRaster{d,3});
    end
    xlabel('Time (ms)');
    ylabel('Firing rate (Hz)');
    legend(strsplit(num2str(optoPulseDuration)));
    
    subplot(1,3,3);
    data = [optoResponse(:,2) optoResponse(:,4)];
    bar(data);
    xlist = [1:1:size(optoResponse,2)];
    xticks(xlist);
    xticklabels(optoPulseDuration)
    xlabel('Opto pulse duration (ms)');
    ylabel('Firing rate (Hz)');
    title(['p(1)= ' num2str(p(1)) ', p(2) = ' num2str(p(2))]);
    
    suptitle(['Unit ' num2str(nData.CellInfo(4)) ' (n = ' num2str(n) '), unitType = ' num2str(nData.CellInfo(6)) ', opto = ' num2str(h)]);
    
    saveas(f1,fullfile(figureDir,['Unit ' num2str(nData.CellInfo(4)) ' opto response plots.fig']));
    saveas(f1,fullfile(figureDir,['Unit ' num2str(nData.CellInfo(4)) ' opto response plots.jpg']));
    close(f1);
end

f2 = figure;
histogram(latencies(:,2));
xlabel('Latency to response (ms)');
ylabel('Number of units');
title('Latency of response to opto stimulation');
saveas(f2,fullfile(newDir,'ResponseLatencyHistogram.fig'));


analysisParams.mouseID = mouseID;
analysisParams.session = session;
analysisParams.date = date;
analysisParams.stimPath = stimPath;
analysisParams.analysisWindow = analysisWindow; %ms
analysisParams.quantificationWindow = quantWindow;
analysisParams.frBinWidth = frBinWidth;

responsiveUnits.optoResponsiveUnits = optoResponsiveUnits;
responsiveUnits.optoResponsiveSingleUnits = optoResponsiveSingleUnits;
responsiveUnits.optoResponsiveMultiUnits = optoResponsiveMultiUnits;
responsiveUnits.totalUnits = totalUnits;
responsiveUnits.singleUnits = singleUnits;
responsiveUnits.multiUnits = multiUnits;
save(fullfile(newDir,'OptoResponseData.mat'),'unitData','analysisParams','responsiveUnits');
