
clear;

experiment = 'EP005';
mouseID = 'AW122';
session = 'Session1';
date = '20200228-1';
stimPath = fullfile('E:\Electrophysiology\',experiment,mouseID,date,'StimInfo',[date '_' mouseID '_LEDFlashes_stimInfo']);

analysisWindow = [-50 100]; %ms relative to stimulus onset
quantWindow = [0 100]; %ms relative to stimulus onset
baselineWindow = [-40 0]; %ms relative to stimulus onset

frBinWidth = 10; %ms

dataFolder = fullfile('E:\Electrophysiology\',experiment,mouseID,date,'SpikeMat');
dataFiles = dir(fullfile(dataFolder,'*LEDflashes*'));

% newDir = fullfile('D:\KiloSort\',mouseID,session,folder,'OptoNoiseResponses');
newDir = fullfile('E:\Electrophysiology\',experiment,mouseID,date,'LED flashes data');
figureDir = fullfile(newDir,'Figures');
if ~exist(newDir)
    mkdir(newDir);
end
if ~exist(figureDir)
    mkdir(figureDir);
end

load(stimPath);
repeats = stimInfo.repeats;

totalUnits = length(dataFiles);

binCount = analysisWindow(2)-analysisWindow(1)-frBinWidth+1; %sliding window
frScalar = 1000/frBinWidth;
quantScalar = 1000/(quantWindow(2)-quantWindow(1));
quantFirstBin = quantWindow(1)-analysisWindow(1)-frBinWidth+1;
quantLastBin = quantWindow(2)-analysisWindow(1)-frBinWidth+1;
baselineFirstBin = baselineWindow(1) - analysisWindow(1)-frBinWidth+1;
baselineLastBin = baselineWindow(2) - analysisWindow(1)-frBinWidth+1;
baselineScalar = 1000/(baselineWindow(2) - baselineWindow(1));
binEdges = analysisWindow(1)+frBinWidth:1:analysisWindow(2);
stimBin = -analysisWindow(1)-frBinWidth + 1;

lightResponsiveUnit = [];
singleUnits = [];
multiUnits = [];
nonNoiseUnits = 0;
responseStats = zeros(0,6); %neuronNumber, neuronN, amplitude, Fano factor, latency, time to peak

rasterColor = [0 0 0];

for n = 1:totalUnits
    n
    nData = load(fullfile(dataFiles(n).folder,dataFiles(n).name));
    
    eventTimes = nData.Events;
    spikeTimes = nData.SpikeData(1,:);
    
    spikeRaster = cell(1,3); %spike times (1), trial number (y-axis) (2), color (3)
    meanResponse = zeros(1,5); %index (1),  mean baseline (2), baseline std (3), mean response (4), response std (5)
    frTrain = zeros(1,binCount);
    frTrainTrials = zeros(repeats,binCount);
    quantTrials = zeros(repeats,2);
    
    spikeFR = zeros(repeats,binCount);
    rasterX = [];
    rasterY = [];
    
    prePost = zeros(repeats,2);
    tempY = 0;
    
    for e = 1:length(eventTimes)
        tempY = tempY+1;
        
        time = eventTimes(e);
        alignedSpikes = (spikeTimes-time)/1000; %temporal resolution of cheetah is in microseconds, converting to milliseconds
        trialSpikes = alignedSpikes(find(alignedSpikes>analysisWindow(1) & alignedSpikes<analysisWindow(end)));
        rasterX = [rasterX trialSpikes];
        rasterY = [rasterY repmat(tempY,[1 length(trialSpikes)])];
        
        for b = 1:binCount %sliding window
            timeMin = analysisWindow(1)+b-1;
            timeMax = timeMin+frBinWidth;
            spikeFR(e,b) = length(find(alignedSpikes>timeMin & alignedSpikes<timeMax));
        end
        prePost(e,1) = histcounts(trialSpikes,[baselineWindow(1) baselineWindow(2)]);
        prePost(e,2) = histcounts(trialSpikes,[quantWindow(1) quantWindow(2)]);
    end
    
    spikeRaster{1,1} = rasterX;
    spikeRaster{1,2} = rasterY;
    spikeRaster{1,3} = rasterColor;
    meanResponse(1,1) = 1;
    meanResponse(1,2) = mean(prePost(:,1))*baselineScalar;
    meanResponse(1,3) = std(prePost(:,1))*baselineScalar/sqrt(repeats);
    meanResponse(1,4) = mean(prePost(:,2))*quantScalar;
    meanResponse(1,5) = std(prePost(:,2))*quantScalar/sqrt(repeats);
    
    frTrainTrials = spikeFR*frScalar;
    frTrain(1,:) = mean(spikeFR,1)*frScalar;
    
    unitData(n).raster = spikeRaster;
    unitData(n).meanResponse = meanResponse;
    unitData(n).frTrain = frTrain;
    unitData(n).frTrainTrials = frTrainTrials;
    unitData(n).type = nData.CellInfo(6);
    unitData(n).neuronNumber = nData.CellInfo(4);
    neuronNumber = nData.CellInfo(4);
    
    if nData.CellInfo(6)==1
        singleUnits = [singleUnits neuronNumber];
        nonNoiseUnits = nonNoiseUnits+1;
    elseif nData.CellInfo(6)==2
        multiUnits = [multiUnits neuronNumber];
        nonNoiseUnits = nonNoiseUnits + 1;
    end
    
    [h p] = ttest2(prePost(:,1)*baselineScalar,prePost(:,2)*quantScalar);h(isnan(h))=0;p(isnan(p))=1;
    if h && nData.CellInfo(6)~=0
        lightResponsiveUnit = [lightResponsiveUnit neuronNumber];
        
        amp = mean(prePost(:,2))*quantScalar - mean(prePost(:,1))*baselineScalar;
        
        fanoFactor = std(prePost(:,2)*quantScalar - prePost(:,1)*baselineScalar)^2 / mean(prePost(:,2)*quantScalar - prePost(:,1)*baselineScalar);
        
        tempTrain = sign(amp)*frTrain(1,:);
        
        baseMean = mean(tempTrain(1,baselineFirstBin:baselineLastBin));
        baseSTD = std(tempTrain(1,baselineFirstBin:baselineLastBin));
        latency = find(tempTrain(1,stimBin:quantLastBin)>(baseMean+2*baseSTD),1);
        
        timeToPeak = find(tempTrain(1,stimBin:quantLastBin)==max(tempTrain(1,stimBin:quantLastBin)),1);
        
        responseStats = [responseStats; neuronNumber n amp fanoFactor latency timeToPeak];
    end
    
    
    f1 = figure;
%     set(f1,'Position',[150 80 1250 700]);
    subplot(1,3,1);
    scatter(spikeRaster{1,1},spikeRaster{1,2},[],spikeRaster{1,3},'.');
    xlim([analysisWindow(1) analysisWindow(end)]);
    xlabel('Time (ms)');
    ylabel('Trial number');
    ylim([0 repeats]);
    
    subplot(1,3,2);
    for r = 1:repeats
        plot(binEdges,frTrainTrials(r,:),'Color',[0 0 0 0.15]);
        hold on;
    end
    plot(binEdges,frTrain,'Color',[0 0 0 1],'LineWidth',2);
    xlim([analysisWindow(1) analysisWindow(end)]);
    xlabel('Time (ms)');
    ylabel('Firing rate (Hz)');
    
    
    subplot(1,3,3);hold on;
    bar([1 2],[meanResponse(1,2) meanResponse(1,4)]);
    hold on;
    errorbar([1 2],[meanResponse(1,2) meanResponse(1,4)],[meanResponse(1,3) meanResponse(1,5)],'Color',[0 0 0],'LineStyle','none');
    xticks([1 2]);
    xticklabels({'Pre','Post'});
    ylabel('Firing rate (Hz)');
    
    suptitle(['Unit ' num2str(nData.CellInfo(4)) ' (n = ' num2str(n) '), unitType = ' num2str(nData.CellInfo(6)) ', p = ' num2str(p) ]);
    
    
            saveas(f1,fullfile(figureDir,['Unit ' num2str(nData.CellInfo(4)) ' response plots.fig']));
            saveas(f1,fullfile(figureDir,['Unit ' num2str(nData.CellInfo(4)) ' response plots.jpg']));
        close(f1);
end

f2 = figure;
set(f2,'Position',[400 200 750 500]);
subplot(1,4,1);
histogram(responseStats(:,3));
title('Response amplitude');
xlabel('\Delta firing rate (Hz)');
subplot(1,4,2);
histogram(responseStats(:,4));
title('Fano factor');
xlabel('Fano factor (Hz)');
subplot(1,4,3);
histogram(responseStats(:,5));
xlabel('Time (ms)');
title('Latency to response');
subplot(1,4,4);
histogram(responseStats(:,6));
xlabel('Time (ms)');
title('Latency to peak FR');
suptitle(['Light responsive units (' num2str(size(responseStats,1)) ' out of ' num2str(nonNoiseUnits) ' units)']);
saveas(f2,fullfile(newDir,['Light responsive units data.fig']));



analysisParams.mouseID = mouseID;
analysisParams.session = session;
analysisParams.date = date;
analysisParams.stimPath = stimPath;
analysisParams.analysisWindow = analysisWindow;
analysisParams.frBinWidth = frBinWidth;
analysisParams.quantWindow = quantWindow;
analysisParams.baselineWindow = baselineWindow;

responsiveUnits.lightResponsiveUnit = lightResponsiveUnit;
responsiveUnits.multiUnits = multiUnits;
responsiveUnits.singleUnits = singleUnits;
responsiveUnits.responseStats = responseStats;


save(fullfile(newDir,'LEDFlashesData.mat'),'unitData','analysisParams','responsiveUnits');
responsiveUnits


