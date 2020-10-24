
clear;

experiment = 'EP007';
mouseID = 'AW142';
session = 'Session1';
date = '20200912';
stimPath = fullfile('D:\Electrophysiology\',experiment);
stimFile = '20200826_noise_100rep_multidB_400k_005_stimInfo';

analysisWindow = [-100 400]; %ms relative to stimulus onset
quantWindow = [0 100]; %ms relative to stimulus onset
baselineWindow = [-90 0]; %ms relative to stimulus onset

frBinWidth = 10; %ms

dataFolder = fullfile('D:\Electrophysiology\',experiment,mouseID,date,'SpikeMat');
dataFilesPre = dir(fullfile(dataFolder,'*ephys_whiteNoise_multidB_pre-injection*'));
dataFilesPost = dir(fullfile(dataFolder,'*ephys_whiteNoise_multidB_post-injection*'));

% newDir = fullfile('D:\KiloSort\',mouseID,session,folder,'OptoNoiseResponses');
newDir = fullfile('D:\Electrophysiology\',experiment,mouseID,date,'WhiteNoiseMultidB');
figureDir = fullfile(newDir,'Figures');
if ~exist(newDir)
    mkdir(newDir);
end
if ~exist(figureDir)
    mkdir(figureDir);
end

load(fullfile(stimPath,stimFile));
uniqueEvents = size(stimInfo.index,1);
indices = stimInfo.index(:,1);
repeats = stimInfo.repeats;

totalUnits = length(dataFilesPre);

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

for n = 1:totalUnits
    n
    nData = load(fullfile(dataFilesPre(n).folder,dataFilesPre(n).name));
    
    eventTimes = nData.Events;
    spikeTimes = nData.SpikeData(1,:);
    
    spikeRaster = cell(2,uniqueEvents,3); %spike times (1), trial number (y-axis) (2), color (3)
    meanResponse = zeros(2,uniqueEvents,5); %index (1),  mean baseline (2), baseline std (3), mean response (4), response std (5)
    trialResponse = zeros(2,uniqueEvents,repeats);
    baselineTrials = zeros(2,uniqueEvents,repeats);
    frTrain = zeros(2,uniqueEvents,binCount);
    frTrainTrials = zeros(2,uniqueEvents,repeats,binCount);
    quantTrials = zeros(2,uniqueEvents,repeats);
    
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
    nData = load(fullfile(dataFilesPre(n).folder,dataFilesPost(n).name));
    
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
            quantTrials(2,u,e) = histcounts(trialSpikes,[quantWindow(1) quantWindow(2)])*quantScalar;
        end
        
        spikeRaster{2,u,1} = rasterX;
        spikeRaster{2,u,2} = rasterY;
        spikeRaster{2,u,3} = color;
        meanResponse(2,u,1) = eventID;
        meanResponse(2,u,2) = mean(prePost(:,1))*baselineScalar;
        meanResponse(2,u,3) = std(prePost(:,1))*baselineScalar/sqrt(repeats);
        meanResponse(2,u,4) = mean(prePost(:,2))*quantScalar;
        meanResponse(2,u,5) = std(prePost(:,2))*quantScalar/sqrt(repeats);
        trialResponse(2,u,:) = prePost(:,2)*quantScalar;
        baselineTrials(2,u,:) = prePost(:,1)*baselineScalar;
        
        frTrainTrials(2,u,:,:) = spikeFR; %JUST NUMBER OF SPIKES IN BIN, NOT FIRING RATE
        frTrain(2,u,:) = mean(spikeFR,1)*frScalar;
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
    set(f1,'Position',[150 80 1250 700]);
    subplot(1,4,1);
    for r = 1:size(spikeRaster,2)
        scatter(spikeRaster{1,r,1},spikeRaster{1,r,2},[],spikeRaster{1,r,3},'.');
        hold on;
    end
    xlim([analysisWindow(1) analysisWindow(end)]);
    xlabel('Time (ms)');
    ylabel('Trial type');
    title('Pre-injection');
    
    subplot(1,4,2);
    for r = 1:size(spikeRaster,2)
        scatter(spikeRaster{2,r,1},spikeRaster{2,r,2},[],spikeRaster{2,r,3},'.');
        hold on;
    end
    xlim([analysisWindow(1) analysisWindow(end)]);
    xlabel('Time (ms)');
    ylabel('Trial type');
    title('Post-injection');
    
    
    subplot(1,4,3);
    hold on;
%     plot(binEdges,squeeze(frTrain(1,uniqueEvents,:)),'Color',rasterColors(uniqueEvents,:));
%     plot(binEdges,squeeze(frTrain(2,uniqueEvents,:)),'Color',rasterColors(uniqueEvents,:));
    plot(binEdges,squeeze(frTrain(1,uniqueEvents,:)),'Color',[0 0 0]);
    plot(binEdges,squeeze(frTrain(2,uniqueEvents,:)),'Color',[0 0 1]);
    xlabel('Time (ms)');
    ylabel('Firing rate (Hz)');
    legend({'Pre-inj','Post-inj'});
    
    subplot(1,4,4);
    hold on;
    xx = stimInfo.intensity;
    yy = meanResponse(1,:,4);
    err = meanResponse(1,:,5);
    errr = errorbar(xx,yy,err);
    errr.Color = [0 0 0];
    yyPost = meanResponse(2,:,4);
    errPost = meanResponse(2,:,5);
    errrPost = errorbar(xx,yyPost,errPost);
    errrPost.Color = [0 0 1];
    xlabel('Sound intensity');
    ylabel('Firing rate (Hz)');
    
    
    suptitle({['Unit ' num2str(nData.CellInfo(4)) ' (n = ' num2str(n)...
        '), unitType = ' num2str(nData.CellInfo(6))] ['Sound responsive = ' num2str(pSound)]});
    
    saveas(f1,fullfile(figureDir,['Unit ' num2str(nData.CellInfo(4)) ' response rasters and PSTH.fig']));
    saveas(f1,fullfile(figureDir,['Unit ' num2str(nData.CellInfo(4)) ' response rasters and PSTH.jpg']));
%     close(f1);
    
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


save(fullfile(newDir,'WhiteNoiseMultidBSaline_Data.mat'),'unitData','analysisParams','responsiveUnits');
responsiveUnits

