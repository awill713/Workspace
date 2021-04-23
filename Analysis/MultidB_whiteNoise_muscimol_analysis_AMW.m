
clear;

experiment = 'EP009';
mouseID = 'AW156';
session = 'Session2';
date = '20201121';
stimPath = fullfile('D:\Electrophysiology\',experiment);
stimFile = '20200826_noise_100rep_multidB_400k_005_stimInfo';

analysisWindow = [-100 400]; %ms relative to stimulus onset
quantWindow = [0 100]; %ms relative to stimulus onset
baselineWindow = [-90 0]; %ms relative to stimulus onset

frBinWidth = 10; %ms

dataFolder = fullfile('D:\Electrophysiology\',experiment,mouseID,date,'SpikeMat');
dataFiles{1} = dir(fullfile(dataFolder,'*ephys_whiteNoise_multidB_pre-injection*'));
dataFiles{2} = dir(fullfile(dataFolder,'*ephys_whiteNoise_multidB_15min_post-injection*'));
dataFiles{3} = dir(fullfile(dataFolder,'*ephys_whiteNoise_multidB_30min_post-injection*'));
dataFiles{4} = dir(fullfile(dataFolder,'*ephys_whiteNoise_multidB_45min_post-injection*'));
dataFiles{5} = dir(fullfile(dataFolder,'*ephys_whiteNoise_multidB_60min_post-injection*'));
times = {'Pre inj','15 min post','30 min post','45 min post','60 min post'};

% newDir = fullfile('D:\KiloSort\',mouseID,session,folder,'OptoNoiseResponses');
newDir = fullfile('D:\Electrophysiology\',experiment,mouseID,date,'WhiteNoiseMultidB');
figureDir = fullfile(newDir,'Figures');
if ~exist(newDir)
    mkdir(newDir);
end
if ~exist(figureDir)
    mkdir(figureDir);
end

load(fullfile(stimPath,mouseID,date,stimFile));
uniqueEvents = size(stimInfo.index,1);
indices = stimInfo.index(:,1);
repeats = stimInfo.repeats;

totalUnits = length(dataFiles{1});

binCount = analysisWindow(2)-analysisWindow(1)-frBinWidth+1; %sliding window
frScalar = 1000/frBinWidth;
quantScalar = 1000/(quantWindow(2)-quantWindow(1));
quantFirstBin = quantWindow(1)-analysisWindow(1)-frBinWidth+1;
quantLastBin = quantWindow(2)-analysisWindow(1)-frBinWidth+1;
baselineFirstBin = baselineWindow(1) - analysisWindow(1)-frBinWidth+1;
baselineLastBin = baselineWindow(2) - analysisWindow(1)-frBinWidth+1;
baselineScalar = 1000/(baselineWindow(2)-baselineWindow(1));
binEdges = analysisWindow(1)+frBinWidth:1:analysisWindow(2);

soundResponsiveUnits = [];
singleUnits = [];
multiUnits = [];

rasterColorMap = 'copper';
map = colormap(rasterColorMap);close;
rasterColors = map(round(linspace(1,length(map),uniqueEvents)),:);
% rasterColors = [0 0 0; 0 1 0; 0 1 1];

timePoints = length(dataFiles);

quantColorMap = 'copper';
qmap = colormap(quantColorMap);close;
quantColors = qmap(round(linspace(1,length(qmap),timePoints)),:);


for n = 1:totalUnits
    n
    nData = load(fullfile(dataFiles{1}(n).folder,dataFiles{1}(n).name));
    
    eventTimes = nData.Events;
    spikeTimes = nData.SpikeData(1,:);
    
    spikeRaster = cell(timePoints,uniqueEvents,3); %spike times (1), trial number (y-axis) (2), color (3)
    meanResponse = zeros(timePoints,uniqueEvents,5); %index (1),  mean baseline (2), baseline std (3), mean response (4), response std (5)
    trialResponse = zeros(timePoints,uniqueEvents,repeats);
    baselineTrials = zeros(timePoints,uniqueEvents,repeats);
    frTrain = zeros(timePoints,uniqueEvents,binCount);
    frTrainTrials = zeros(timePoints,uniqueEvents,repeats,binCount);
    quantTrials = zeros(timePoints,uniqueEvents,repeats);
    
    %% Pre injection data
    for u = 1:uniqueEvents
        eventID = indices(u);
        eventsOfInterest = find(stimInfo.order==eventID);
        
        color = rasterColors(u,:);
        
        spikeFR = zeros(length(eventsOfInterest),binCount);
        rasterX = [];
        rasterY = [];
        
        prePost = zeros(length(eventsOfInterest),2);
        tempY = mod(u-1,uniqueEvents)*repeats;
        for e = 1:length(eventsOfInterest)
            tempY = tempY+1;
            
            time = eventTimes(eventsOfInterest(e));
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
            quantTrials(1,u,e) = histcounts(trialSpikes,[quantWindow(1) quantWindow(2)])*quantScalar;
        end
        
        spikeRaster{1,u,1} = rasterX;
        spikeRaster{1,u,2} = rasterY;
        spikeRaster{1,u,3} = color;
        meanResponse(1,u,1) = eventID;
        meanResponse(1,u,2) = mean(prePost(:,1))*baselineScalar;
        meanResponse(1,u,3) = std(prePost(:,1))*baselineScalar/sqrt(repeats);
        meanResponse(1,u,4) = mean(prePost(:,2))*quantScalar;
        meanResponse(1,u,5) = std(prePost(:,2))*quantScalar/sqrt(repeats);
        trialResponse(1,u,:) = prePost(:,2)*quantScalar;
        baselineTrials(1,u,:) = prePost(:,1)*baselineScalar;
        
        frTrainTrials(1,u,:,:) = spikeFR; %JUST NUMBER OF SPIKES IN BIN, NOT FIRING RATE
        frTrain(1,u,:) = mean(spikeFR,1)*frScalar;
    end
    
    neuronNumber = nData.CellInfo(4);
    anovaMat = [];
    for un = 1:uniqueEvents
        col = squeeze(trialResponse(1,un,:));
        anovaMat = [anovaMat col];
    end
    pSound = anova1(anovaMat,[],'off');
    if pSound<0.05
        soundResponsiveUnits = [soundResponsiveUnits neuronNumber];
    end
    
    %% Post injection data
    for tp = 2:timePoints
        
    nData = load(fullfile(dataFiles{tp}(n).folder,dataFiles{tp}(n).name));
    
    eventTimes = nData.Events;
    spikeTimes = nData.SpikeData(1,:);
    
    for u = 1:uniqueEvents
        eventID = indices(u);
        eventsOfInterest = find(stimInfo.order==eventID);
        
        color = rasterColors(u,:);
        
        spikeFR = zeros(length(eventsOfInterest),binCount);
        rasterX = [];
        rasterY = [];
        
        prePost = zeros(length(eventsOfInterest),2);
        tempY = mod(u-1,uniqueEvents)*repeats;
        for e = 1:length(eventsOfInterest)
            tempY = tempY+1;
            
            time = eventTimes(eventsOfInterest(e));
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
            quantTrials(tp,u,e) = histcounts(trialSpikes,[quantWindow(1) quantWindow(2)])*quantScalar;
        end
        
        spikeRaster{tp,u,1} = rasterX;
        spikeRaster{tp,u,2} = rasterY;
        spikeRaster{tp,u,3} = color;
        meanResponse(tp,u,1) = eventID;
        meanResponse(tp,u,2) = mean(prePost(:,1))*baselineScalar;
        meanResponse(tp,u,3) = std(prePost(:,1))*baselineScalar/sqrt(repeats);
        meanResponse(tp,u,4) = mean(prePost(:,2))*quantScalar;
        meanResponse(tp,u,5) = std(prePost(:,2))*quantScalar/sqrt(repeats);
        trialResponse(tp,u,:) = prePost(:,2)*quantScalar;
        baselineTrials(tp,u,:) = prePost(:,1)*baselineScalar;
        
        frTrainTrials(tp,u,:,:) = spikeFR; %JUST NUMBER OF SPIKES IN BIN, NOT FIRING RATE
        frTrain(tp,u,:) = mean(spikeFR,1)*frScalar;
    end
    end
    
    %% Save and plot
    unitData(n).raster = spikeRaster;
    unitData(n).meanResponse = meanResponse;
    unitData(n).frTrain = frTrain;
    unitData(n).frTrainSTD = squeeze(std(frTrainTrials,[],2));
    unitData(n).frTrainTrials = frTrainTrials; %JUST NUMBER OF SPIKES IN BIN, NOT FIRING RATE
    unitData(n).trialResponse = trialResponse;
    unitData(n).type = nData.CellInfo(6);
    unitData(n).neuronNumber = nData.CellInfo(4);
    
    
    if nData.CellInfo(6)==1
        singleUnits = [singleUnits neuronNumber];
    elseif nData.CellInfo(6)==2
        multiUnits = [multiUnits neuronNumber];
    end
    
    
    
    f1 = figure;
    set(f1,'Position',[150 150 1250 500]);
    for tp = 1:timePoints
        subplot(1,timePoints,tp);
    for r = 1:size(spikeRaster,2)
        scatter(spikeRaster{tp,r,1},spikeRaster{tp,r,2},[],spikeRaster{tp,r,3},'.');
        hold on;
    end
    xlim([analysisWindow(1) analysisWindow(end)]);
    xlabel('Time (ms)');
    ylabel('Trial type');
    title(times{tp});
    end
    suptitle({['Unit ' num2str(nData.CellInfo(4)) ' (n = ' num2str(n)...
        '), unitType = ' num2str(nData.CellInfo(6))] ['Sound responsive = ' num2str(pSound)]});
    
    saveas(f1,fullfile(figureDir,['Unit ' num2str(nData.CellInfo(4)) ' rasters.fig']));
    saveas(f1,fullfile(figureDir,['Unit ' num2str(nData.CellInfo(4)) ' rasters.jpg']));
    close(f1);
    
    f2 = figure;
    set(f2,'Position',[150 150 1000 400]);
    subplot(1,2,1);
    hold on;
%     plot(binEdges,squeeze(frTrain(1,uniqueEvents,:)),'Color',rasterColors(uniqueEvents,:));
%     plot(binEdges,squeeze(frTrain(2,uniqueEvents,:)),'Color',rasterColors(uniqueEvents,:));
    plot(binEdges,squeeze(frTrain(1,uniqueEvents,:)),'Color',[0 0 0]);
    plot(binEdges,squeeze(frTrain(timePoints,uniqueEvents,:)),'Color',[0 0 1]);
    xlabel('Time (ms)');
    ylabel('Firing rate (Hz)');
    legend({times{1},times{timePoints}});
    
    subplot(1,2,2);
    hold on;
    for tp = 1:timePoints
        xx = stimInfo.intensity;
        yy = meanResponse(tp,:,4);
        err = meanResponse(tp,:,5);
        errr = errorbar(xx,yy,err);
        errr.Color = quantColors(tp,:);
        errr.LineWidth = 2;
    end
    legend(times);
    xlabel('Sound intensity');
    ylabel('Firing rate (Hz)');
    
    
    suptitle({['Unit ' num2str(nData.CellInfo(4)) ' (n = ' num2str(n)...
        '), unitType = ' num2str(nData.CellInfo(6))] ['Sound responsive = ' num2str(pSound)]});
    
    saveas(f2,fullfile(figureDir,['Unit ' num2str(nData.CellInfo(4)) ' PSTH and quant.fig']));
    saveas(f2,fullfile(figureDir,['Unit ' num2str(nData.CellInfo(4)) ' PSTH and quant.jpg']));
    close(f2);
    
end


analysisParams.mouseID = mouseID;
analysisParams.session = session;
analysisParams.date = date;
analysisParams.stimPath = stimPath;
analysisParams.analysisWindow = analysisWindow;
analysisParams.frBinWidth = frBinWidth;
analysisParams.quantWindow = quantWindow;
analysisParams.baselineWindow = baselineWindow;

responsiveUnits.soundResponsiveUnits = intersect(soundResponsiveUnits,[singleUnits multiUnits]);
responsiveUnits.multiUnits = multiUnits;
responsiveUnits.singleUnits = singleUnits;


save(fullfile(newDir,'WhiteNoiseMultidBMuscimol_Data.mat'),'unitData','analysisParams','responsiveUnits');
responsiveUnits


meanQuants = cell(1,timePoints);
for nn = 1:totalUnits
    neuronNumber = unitData(nn).neuronNumber;
    
    if ismember(neuronNumber,responsiveUnits.soundResponsiveUnits)
        for tp = 1:timePoints
            mq = squeeze(unitData(nn).meanResponse(tp,:,4));
            meanQuants{tp} = [meanQuants{tp}; mq];
        end
    end
end
for tp = 1:timePoints
    meanMean(tp,:) = mean(meanQuants{tp});
    stdMean(tp,:) = std(meanQuants{tp})./sqrt(size(meanQuants{tp},1));
end
f3 = figure;
for tp = 1:timePoints
    xx = stimInfo.intensity;
    yy = meanMean(tp,:);
    erx = stdMean(tp,:);
    erplot = errorbar(xx,yy,erx);
    erplot.Color = quantColors(tp,:);
    erplot.LineWidth = 2;
    hold on;
end
legend(times);
xlabel('Sound intensity');
ylabel('Firing rate (Hz)');