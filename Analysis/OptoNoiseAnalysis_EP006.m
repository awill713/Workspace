
clear;

experiment = 'EP006';
mouseID = 'AW131';
session = 'Session1';
date = '20200731';
stimPath = fullfile('D:\Electrophysiology\',experiment,mouseID,date,'StimInfo');
stimFile = '20200720_noise_opto_100rep_70dB_400k_005_stimInfo';
% stimFile = '20200720_noise_opto_100rep_90dB_400k_005_stimInfo';
% stimFile = '20200731_noise_opto_100rep_50dB_400k_005_stimInfo';

analysisWindow = [-100 400]; %ms relative to stimulus onset
quantWindow = [0 100]; %ms relative to stimulus onset
baselineWindow = [-90 0]; %ms relative to stimulus onset

frBinWidth = 10; %ms

dataFolder = fullfile('D:\Electrophysiology\',experiment,mouseID,date,'SpikeMat');
dataFiles = dir(fullfile(dataFolder,'*ephys_whiteNoise_70dB_laser_100repeats*'));

% newDir = fullfile('D:\KiloSort\',mouseID,session,folder,'OptoNoiseResponses');
newDir = fullfile('D:\Electrophysiology\',experiment,mouseID,date,'OptoNoise - 70dB');
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
laserResponsiveUnits = []; %units whose laser is different from baseline
interactUnits = []; %units whose sound+laser is different from sound
singleUnits = [];
multiUnits = [];

% rasterColorMap = 'hsv';
% map = colormap(rasterColorMap);close;
% rasterColors = map(round(linspace(1,length(map),length(orientations))),:);
rasterColors = [0 0 0; 0 1 0; 0 1 1];

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
    
    [hSound pSound] = ttest2(baselineTrials(1,:),trialResponse(1,:));hSound(isnan(hSound))=0;
    if hSound
        soundResponsiveUnits = [soundResponsiveUnits neuronNumber];
    end
    
    [hLaser pLaser] = ttest2(baselineTrials(2,:),trialResponse(2,:));hLaser(isnan(hLaser))=0;
    if hLaser
        laserResponsiveUnits = [laserResponsiveUnits neuronNumber];
    end
    
    [hInteract pInteract] = ttest2(trialResponse(3,:),trialResponse(1,:));hInteract(isnan(hInteract))=0;
    if hInteract
        interactUnits = [interactUnits neuronNumber];
    end
    
    
    f1 = figure;
    set(f1,'Position',[150 80 1250 700]);
    subplot(1,3,1);
    for r = 1:size(spikeRaster,1)
        scatter(spikeRaster{r,1},spikeRaster{r,2},[],spikeRaster{r,3},'.');
        hold on;
    end
    xlim([analysisWindow(1) analysisWindow(end)]);
    xlabel('Time (ms)');
    ylabel('Trial type');
    
    
    subplot(1,3,2);
    hold on;
    for frt = 1:uniqueEvents
        plot(binEdges,frTrain(frt,:),'Color',rasterColors(frt,:));
    end
    xlabel('Time (ms)');
    ylabel('Firing rate (Hz)');
    legend({'Noise','Laser','Noise+Laser'});
    
    subplot(1,3,3);
    hold on;
    data = [mean(meanResponse(:,2)) meanResponse(:,4)'];
    error = [mean(meanResponse(:,3)) meanResponse(:,5)'];
    bar(1:4,data);
    er = errorbar(1:4,data,error);
    er.Color = [0 0 0];
    er.LineStyle = 'none';
    xticks(1:4);
    xticklabels({'Base','Noise', 'Laser', 'Noise+Laser'});
    xtickangle(45);
    ylabel('Firing rate (Hz)');
    
    
    suptitle({['Unit ' num2str(nData.CellInfo(4)) ' (n = ' num2str(n)...
        '), unitType = ' num2str(nData.CellInfo(6))] ['Sound responsive = ' num2str(pSound)...
        ', Laser responsive = ' num2str(pLaser) ', Laser mods sound = ' num2str(pInteract)]});
    
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
responsiveUnits.interactUnits = intersect(interactUnits,[singleUnits multiUnits]);
responsiveUnits.multiUnits = multiUnits;
responsiveUnits.singleUnits = singleUnits;




save(fullfile(newDir,'OptoNoise70dB_Data.mat'),'unitData','analysisParams','responsiveUnits');
responsiveUnits

