
clear;

experiment = 'EP011';
mouseID = 'AW178';
session = 'Session3';
date = '20210327-1';
stimPath = fullfile('D:\Electrophysiology\',experiment,mouseID,date,'StimInfo');
stimFile = '20210219_multidBNoise_opto_50rep_400k_005_stimInfo';

analysisWindow = [-100 400]; %ms relative to stimulus onset
quantWindow = [0 100]; %ms relative to stimulus onset
baselineWindow = [-90 0]; %ms relative to stimulus onset

frBinWidth = 10; %ms

dataFolder = fullfile('D:\Electrophysiology\',experiment,mouseID,date,'SpikeMat');
dataFiles = dir(fullfile(dataFolder,'*ephys_multidBNoise_laser*'));

newDir = fullfile('D:\Electrophysiology\',experiment,mouseID,date,'WhiteNoiseMultidB_laser');
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
uniqueSounds = length(stimInfo.intensity);

totalUnits = length(dataFiles);

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
laserResponsiveUnits = [];
singleUnits = [];
multiUnits = [];

rasterColorMap = 'copper';
map = colormap(rasterColorMap);close;
rasterColors = map(round(linspace(1,length(map),uniqueSounds)),:);
% rasterColors = [0 0 0; 0 1 0; 0 1 1];

for n = 1:totalUnits
    n
    nData = load(fullfile(dataFiles(n).folder,dataFiles(n).name));
    
    eventTimes = nData.Events;
    spikeTimes = nData.SpikeData(1,:);
    
    spikeRaster = cell(uniqueEvents,3); %spike times (1), trial number (y-axis) (2), color (3)
    meanResponse = zeros(uniqueEvents,5); %index (1),  mean baseline (2), baseline std (3), mean response (4), response std (5)
    trialResponse = zeros(uniqueEvents,repeats);
    baselineTrials = zeros(uniqueEvents,repeats);
    frTrain = zeros(uniqueEvents,binCount);
    frTrainTrials = zeros(uniqueEvents,repeats,binCount);
    quantTrials = zeros(uniqueEvents,repeats);
    
    for u = 1:uniqueEvents
        eventID = indices(u);
        eventsOfInterest = find(stimInfo.order==eventID);
        
        color = rasterColors(mod(u-1,uniqueSounds)+1,:);
        
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
            quantTrials(u,e) = histcounts(trialSpikes,[quantWindow(1) quantWindow(2)])*quantScalar;
        end
        
        spikeRaster{u,1} = rasterX;
        spikeRaster{u,2} = rasterY;
        spikeRaster{u,3} = color;
        meanResponse(u,1) = eventID;
        meanResponse(u,2) = mean(prePost(:,1))*baselineScalar;
        meanResponse(u,3) = std(prePost(:,1))*baselineScalar/sqrt(repeats);
        meanResponse(u,4) = mean(prePost(:,2))*quantScalar;
        meanResponse(u,5) = std(prePost(:,2))*quantScalar/sqrt(repeats);
        trialResponse(u,:) = prePost(:,2)*quantScalar;
        baselineTrials(u,:) = prePost(:,1)*baselineScalar;
        
        frTrainTrials(u,:,:) = spikeFR; %JUST NUMBER OF SPIKES IN BIN, NOT FIRING RATE
        frTrain(u,:) = mean(spikeFR,1)*frScalar;
    end
    
    unitData(n).raster = spikeRaster;
    unitData(n).meanResponse = meanResponse;
    unitData(n).frTrain = frTrain;
    unitData(n).frTrainSTD = squeeze(std(frTrainTrials,[],2));
    unitData(n).frTrainTrials = frTrainTrials; %JUST NUMBER OF SPIKES IN BIN, NOT FIRING RATE
    unitData(n).trialResponse = trialResponse;
    unitData(n).type = nData.CellInfo(6);
    unitData(n).neuronNumber = nData.CellInfo(4);
    neuronNumber = nData.CellInfo(4);
    
    if nData.CellInfo(6)==1
        singleUnits = [singleUnits neuronNumber];
    elseif nData.CellInfo(6)==2
        multiUnits = [multiUnits neuronNumber];
    end
    
    anovaMat = [];
    for un = 1:uniqueSounds
        col = trialResponse(un,:)';
        laserCol = trialResponse(un+uniqueSounds,:)';
        anovaMat = [anovaMat [col; laserCol]];
    end
    pValues = anova2(anovaMat,repeats,'off');
    if pValues(1)<0.05
        soundResponsiveUnits = [soundResponsiveUnits neuronNumber];
    end
    if pValues(2)<0.05
        laserResponsiveUnits = [laserResponsiveUnits neuronNumber];
    end
    
    f1 = figure;
    set(f1,'Position',[150 80 1250 700]);
    subplot(1,4,1);
    for r = 1:uniqueSounds
        scatter(spikeRaster{r,1},spikeRaster{r,2},[],spikeRaster{r,3},'.');
        hold on;
    end
    xlim([analysisWindow(1) analysisWindow(end)]);
    xlabel('Time (ms)');
    ylabel('Trial type');
    title('Sound only');
    
    subplot(1,4,2);
    for r = 1:uniqueSounds
        scatter(spikeRaster{r+uniqueSounds,1},spikeRaster{r+uniqueSounds,2},[],spikeRaster{r+uniqueSounds,3},'.');
        hold on;
    end
    xlim([analysisWindow(1) analysisWindow(end)]);
    xlabel('Time (ms)');
    ylabel('Trial type');
    title('Sound with laser');
    
    
    subplot(1,4,3);
    hold on;
    baselineTrain = frTrain(1,:);
    laserTrain = frTrain(1+uniqueSounds,:);
    soundTrain = frTrain(uniqueSounds,:);
    soundLaserTrain = frTrain(uniqueEvents,:);
    plot(binEdges,baselineTrain,':','Color',rasterColors(1,:));
    plot(binEdges,laserTrain,':','Color',[0 1 0]);
    plot(binEdges,soundTrain,'Color',rasterColors(1,:));
    plot(binEdges,soundLaserTrain,'Color',[0 1 0]);
    xlabel('Time (ms)');
    ylabel('Firing rate (Hz)');
    legend('Baseline','Laser only','Noise','Noise+laser');
    
    subplot(1,4,4);
    hold on;
    xx = stimInfo.intensity;
    yy = meanResponse(1:uniqueSounds,4);
    err = meanResponse(1:uniqueSounds,5);
    errorbar(xx,yy,err,'Color',[0 0 0]);
    lyy = meanResponse((1:uniqueSounds)+uniqueSounds,4);
    lerr = meanResponse((1:uniqueSounds)+uniqueSounds,5);
    errorbar(xx,lyy,lerr,'Color',[0 1 0]);
    xlabel('Sound intensity');
    ylabel('Firing rate (Hz)');
    legend({'Sound only','Sound with laser'});
    
    
    suptitle({['Unit ' num2str(nData.CellInfo(4)) ' (n = ' num2str(n)...
        '), unitType = ' num2str(nData.CellInfo(6))] ['Sound responsive = ' num2str(pValues(1))...
        ', Laser responsive = ' num2str(pValues(2))]});
    
    saveas(f1,fullfile(figureDir,['Unit ' num2str(nData.CellInfo(4)) ' response rasters and PSTH.fig']));
    saveas(f1,fullfile(figureDir,['Unit ' num2str(nData.CellInfo(4)) ' response rasters and PSTH.jpg']));
    close(f1);
    
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
responsiveUnits.laserResponsiveUnits = intersect(laserResponsiveUnits,[singleUnits multiUnits]);
responsiveUnits.multiUnits = multiUnits;
responsiveUnits.singleUnits = singleUnits;


save(fullfile(newDir,'WhiteNoiseMultidB_Data.mat'),'unitData','analysisParams','responsiveUnits');
responsiveUnits

