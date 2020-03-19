
clear;

experiment = 'EP004';
mouseID = 'AW115';
session = 'Session2';
date = '20200222';
stimPath = fullfile('E:\Electrophysiology\',experiment,mouseID,date,'StimInfo',[date '_' mouseID '_AVstaticMultiAmpTones_stimInfo']);

analysisWindow = [-50 200];
quantWindow = [0 100]; %ms relative to stimulus onset
baselineWindow = [-40 0]; %ms relative to stimulus onset
rasterHeatMap = 'jet';

frBinWidth = 10; %sliding window

% dataFolder = fullfile('D:\KiloSort\',mouseID,session,date,'SpikeMat');
dataFolder = fullfile('E:\Electrophysiology\',experiment,mouseID,date,'SpikeMat');
dataFiles = dir(fullfile(dataFolder,'*AV_static_multiAmpTones*'));

% newDir = fullfile('D:\KiloSort\',mouseID,session,date,'OptoMultiAmpToneResponses');
newDir = fullfile('E:\Electrophysiology\',experiment,mouseID,date,'AVStaticMultiAmpTones');
figureDir = fullfile(newDir,'Figures');
if ~exist(newDir)
    mkdir(newDir);
end
if ~exist(figureDir)
    mkdir(figureDir);
end

load(stimPath);
uniqueEvents = size(stimInfo.index,1);
indices = stimInfo.index(:,1);
repeats = stimInfo.repeats;
frequencies = stimInfo.frequencies;
intensities = stimInfo.intensities;
uniqueFreq = length(frequencies);
uniqueIntensities= length(intensities);

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

frequencySelectiveUnits = [];
soundResponsiveUnits = [];
lightResponsiveUnits = [];
interactUnits = [];
singleUnits = [];
multiUnits = [];

map = colormap(rasterHeatMap);close;


for n = 1:totalUnits
    
    n
    nData = load(fullfile(dataFiles(n).folder,dataFiles(n).name));
    
    eventTimes = nData.Events;
    spikeTimes = nData.SpikeData(1,:);
    
    spikeRaster = cell(uniqueEvents,3); %spike times (1), trial number (y-axis) (2), color (3)
    meanResponse = zeros(uniqueEvents,5); %index (1),  mean baseline (2), baseline std (3), mean response (4), response std (5)
    frTrain = zeros(uniqueEvents,binCount);
    frTrainTrials = zeros(uniqueEvents,repeats,binCount);
    tuningMat = zeros(uniqueFreq,uniqueIntensities,2);
    
    for u = 1:uniqueEvents
        eventID = indices(u);
        eventsOfInterest = find(stimInfo.order==eventID);
        
        freq = stimInfo.index(eventID,2);
        freqIdx = find(frequencies==freq);
        color = map(round(freqIdx/length(frequencies)*length(map)),:);
        
        intensity = stimInfo.index(eventID,3);
        intenseIdx = find(intensities==intensity);
        
        lightStatus = stimInfo.index(eventID,4);
        lightIdx = lightStatus+1;
        
        
        spikeFR = zeros(length(eventsOfInterest),binCount);
        rasterX = [];
        rasterY = [];
        
        prePost = zeros(length(eventsOfInterest),2);
        tempY = (intenseIdx-1)*uniqueFreq*repeats + (freqIdx-1)*repeats;
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
        end
        
        spikeRaster{u,1} = rasterX;
        spikeRaster{u,2} = rasterY;
        spikeRaster{u,3} = color;
        meanResponse(u,1) = eventID;
        meanResponse(u,2) = mean(prePost(:,1))*baselineScalar;
        meanResponse(u,3) = std(prePost(:,1))*baselineScalar/sqrt(repeats);
        meanResponse(u,4) = mean(prePost(:,2))*quantScalar;
        meanResponse(u,5) = std(prePost(:,2))*quantScalar/sqrt(repeats);
        tuningMat(freqIdx,intenseIdx,lightIdx) = mean(prePost(:,2))*quantScalar;
        
        frTrainTrials(u,:,:) = spikeFR*frScalar;
        frTrain(u,:) = mean(spikeFR,1)*frScalar;
    end
    
    noLightVector = reshape(tuningMat(:,:,1),[uniqueFreq*uniqueIntensities 1]);
    yesLightVector = reshape(tuningMat(:,:,2),[uniqueFreq*uniqueIntensities 1]);
    fit = polyfit(noLightVector,yesLightVector,1);
    
    unitData(n).raster = spikeRaster;
    unitData(n).meanResponse = meanResponse;
    unitData(n).tuningMat = tuningMat;
    unitData(n).frTrain = frTrain;
    unitData(n).frTrainTrials = frTrainTrials;
    unitData(n).lightEffectRegression = fit;
    unitData(n).type = nData.CellInfo(6);
    unitData(n).neuronNumber = nData.CellInfo(4);
    neuronNumber = nData.CellInfo(4);
    
    if nData.CellInfo(6)==1
        singleUnits = [singleUnits neuronNumber];
    elseif nData.CellInfo(6)==2
        multiUnits = [multiUnits neuronNumber];
    end
    
    
    p = anovan(meanResponse(:,4),{stimInfo.index(:,2),stimInfo.index(:,3),stimInfo.index(:,4)},'model','interaction','varnames',{'freq','amp','light'},'display','off');
    if p(1)<0.05 || p(4)<0.05
        frequencySelectiveUnits = [frequencySelectiveUnits neuronNumber];
    end
    if p(2)<0.05
        soundResponsiveUnits = [soundResponsiveUnits neuronNumber];
    end
    if p(3)<0.05
        lightResponsiveUnits = [lightResponsiveUnits neuronNumber];
    end
    if p(5)<0.05 || p(6)<0.05
        interactUnits = [interactUnits neuronNumber];
    end
    
        
    
    
    low = min(min(min(tuningMat)));
    high = max(max(max(tuningMat)));
    
    f1 = figure;
    set(f1,'Position',[100 100 1200 600])
    subplot(2,4,1); colormap 'jet';
    for r = 1:uniqueFreq*uniqueIntensities
        scatter(spikeRaster{r,1},spikeRaster{r,2},[],spikeRaster{r,3},'.');
        hold on;
    end
    xlim([analysisWindow(1) analysisWindow(end)]);
    xlabel('Time (ms)');
    ylabel('Trial type (freq and intensity)');
    yticks([]);
    title('Sound')
    
    subplot(2,4,2); colormap 'jet';
    for r = uniqueFreq*uniqueIntensities+1:uniqueEvents
        scatter(spikeRaster{r,1},spikeRaster{r,2},[],spikeRaster{r,3},'.');
        hold on;
    end
    xlim([analysisWindow(1) analysisWindow(end)]);
    xlabel('Time (ms)');
    ylabel('Trial type (freq and intensity)');
    yticks([]);
    title('Sound/Light')
    
    subplot(2,4,5);
    imagesc(flip(tuningMat(:,:,1)'));
    xlist = [1 5:5:uniqueFreq];
    xticks(xlist);
    xticklabels(round(frequencies(xlist)/1000,1))
    xlabel('Frequency (kHz)');
    ylist = [1 round(uniqueIntensities/2) uniqueIntensities];
    yticks(ylist);
    yticklabels(flip(intensities(ylist)));
    ylabel('SPL (dB)');
    cb = colorbar;
    caxis([low high]);
    ylabel(cb,'Firing rate (Hz)');
    title('Sound');
    
    subplot(2,4,6);
    imagesc(flip(tuningMat(:,:,2)'));
    xlist = [1 5:5:uniqueFreq];
    xticks(xlist);
    xticklabels(round(frequencies(xlist)/1000,1))
    xlabel('Frequency (kHz)');
    ylist = [1 round(uniqueIntensities/2) uniqueIntensities];
    yticks(ylist);
    yticklabels(flip(intensities(ylist)));
    ylabel('SPL (dB)');
    cb = colorbar;
    caxis([low high]);
    ylabel(cb,'Firing rate (Hz)');
    title('Sound/Light');
    
    subplot(2,4,7);
    diffMat = flip(tuningMat(:,:,2)') - flip(tuningMat(:,:,1)');
    extreme = max(max(abs(diffMat)));
    imagesc(diffMat);
    xlist = [1 5:5:uniqueFreq];
    xticks(xlist);
    xticklabels(round(frequencies(xlist)/1000,1))
    xlabel('Frequency (kHz)');
    ylist = [1 round(uniqueIntensities/2) uniqueIntensities];
    yticks(ylist);
    yticklabels(flip(intensities(ylist)));
    ylabel('SPL (dB)');
    cb = colorbar;
    caxis([-extreme extreme]);
    ylabel(cb,'Firing rate (Hz)');
    title('Sound/Light - Sound');
    
    subplot(2,4,8);
    scatter(noLightVector,yesLightVector);
    unity = refline(1,0); unity.Color = 'k'; unity.LineStyle = '--';
    fitLine = refline(fit(1),fit(2)); fitLine.Color = 'b';
    xlabel('Tone evoked FR, sound (Hz)');
    ylabel('Tone evoked FR, sound/light (Hz)');
    title(['Slope = ' num2str(fit(1)) ', y-int = ' num2str(fit(2))]);
    
    subplot(2,4,3);
    tcScalar = 20/max(tuningMat(:));
    for spl = 1:size(tuningMat,2)
        semilogx(frequencies/1000,tuningMat(:,spl,1)*tcScalar + (spl-1)*20,'Color',[0 0 0]);
        hold on;
        semilogx(frequencies/1000,tuningMat(:,spl,2)*tcScalar + (spl-1)*20,'Color',[0 0 1]);
        hold on;
        line([frequencies(1)/1000 frequencies(end)/1000],[(spl-1)*20 (spl-1)*20],'LineStyle','--');
    end
    xlabel('Frequency (kHz)');
    ylabel('\Delta FR (Hz)');
    xlist = [1 5:5:uniqueFreq];
    xticks(round(frequencies(xlist)/1000,1));
    yticks([]);
    title('Tuning curves');
%     legend(strsplit(num2str(intensities)));
    
    subplot(2,4,4);
    tcScalar = 20/max(tuningMat(:));
    for spl = 1:size(tuningMat,2)
        semilogx(frequencies/1000,(tuningMat(:,spl,2)-tuningMat(:,spl,1))*tcScalar + (spl-1)*20,'Color',[0 0 0]);
        hold on;
        line([frequencies(1)/1000 frequencies(end)/1000],[(spl-1)*20 (spl-1)*20],'LineStyle','--');
    end
    xlabel('Frequency (kHz)');
    ylabel('\Delta FR (Hz)');
    xlist = [1 5:5:uniqueFreq];
    xticks(round(frequencies(xlist)/1000,1));
    yticks([]);
    title('\Delta tuning curves');
%     legend(strsplit(num2str(intensities)));


    suptitle({['Unit ' num2str(nData.CellInfo(4)) ' (n = ' num2str(n)...
        '), unitType = ' num2str(nData.CellInfo(6))] ['Sound responsive = ' num2str(p(2))...
        ', Light responsive = ' num2str(p(3)) ', Frequency selective = ' num2str(p(1))...
        ', Light-modulated = ' num2str(p(5)) ', ' num2str(p(6))]});
    
    saveas(f1,fullfile(figureDir,['Unit ' num2str(nData.CellInfo(4)) ' response plots.fig']));
    saveas(f1,fullfile(figureDir,['Unit ' num2str(nData.CellInfo(4)) ' response plots.jpg']));
%     close(f1);
    
end


analysisParams.mouseID = mouseID;
analysisParams.session = session;
analysisParams.folder = date;
analysisParams.stimPath = stimPath;
analysisParams.analysisWindow = analysisWindow; %ms
analysisParams.quantificationWindow = quantWindow;
analysisParams.frBinWidth = frBinWidth;
analysisParams.baselineWindow = baselineWindow;

responsiveUnits.soundResponsiveUnits = soundResponsiveUnits;
responsiveUnits.lightResponsiveUnits = lightResponsiveUnits;
responsiveUnits.interactUnits = interactUnits;
responsiveUnits.frequencySelectiveUnits = frequencySelectiveUnits;
responsiveUnits.multiUnits = multiUnits;
responsiveUnits.singleUnits = singleUnits;

    
save(fullfile(newDir,'AVStaticMultiAmpTonesData.mat'),'unitData','analysisParams','responsiveUnits');


freqToOctave = ceil(log(2)/log(frequencies(2)/frequencies(1)));
sideBands = cell(uniqueIntensities,2*freqToOctave+1,2);
diffSideBands = cell(uniqueIntensities,2*freqToOctave+1,2);
bestFrequencies = [];
tcStacks = cell(uniqueIntensities,2);
neuronCount = 0;
for n = 1:totalUnits
    if unitData(n).type~=0 && ismember(unitData(n).neuronNumber,responsiveUnits.soundResponsiveUnits)
        neuronCount = neuronCount+1;
        
        tunMat = unitData(n).tuningMat;
        maxVal = max(tunMat(:));
        minVal= min(tunMat(:));
        [val ind] = max(tunMat(:,uniqueIntensities,1));
        
        for i = 1:uniqueIntensities
            tCurveAud = tunMat(:,i,1);
            tCurveAudVis = tunMat(:,i,2);
            
            tCurveAudVis = (tCurveAudVis-minVal)/(maxVal - minVal);
            tCurveAud = (tCurveAud-minVal)/(maxVal - minVal);
            
%             [val ind] = max(tCurveAudVis);
            indexRange = (ind-freqToOctave):ind+freqToOctave;
            startInd = max(1,ind-freqToOctave);
            startBin = find(indexRange == startInd);
            endInd = min(uniqueFreq,ind+freqToOctave);
            endBin = find(indexRange==endInd);
            bin = startInd-1;
            for b = startBin:endBin
                bin = bin+1;
                sideBands{i,b,1} = [sideBands{i,b,1} tCurveAud(bin)];
                sideBands{i,b,2} = [sideBands{i,b,2} tCurveAudVis(bin)];
            end
            if i == uniqueIntensities
                bestFrequencies = [bestFrequencies ind];
            end
            
            tcStacks{i,1} = [tcStacks{i,1}; tCurveAud'];
            tcStacks{i,2} = [tcStacks{i,2}; tCurveAudVis'];
                
        end
    end
end
sbCurves = zeros(uniqueIntensities,size(sideBands,2),2);
for i = 1:uniqueIntensities
    tempCurveA = zeros(1,size(sideBands,2));
    tempCurveAV = zeros(1,size(sideBands,2));
    
    for s = 1:length(tempCurveA)
        tempCurveA(s) = mean(sideBands{i,s,1});
        tempCurveAV(s) = mean(sideBands{i,s,2});
    end
    
    sbCurves(i,:,1) = tempCurveA;
    sbCurves(i,:,2) = tempCurveAV;
end
f2 = figure;
set(f2,'Position',[25 200 1500 400])
for i = 1:uniqueIntensities
    subplot(1,uniqueIntensities,i);hold on;
    plot((1:size(sbCurves,2))-freqToOctave-1,sbCurves(i,:,1),'Color',[0 0 0],'LineWidth',2);
    plot((1:size(sbCurves,2))-freqToOctave-1,sbCurves(i,:,2),'Color',[0 0 1],'LineWidth',2);
    title(['SPL = ' num2str(intensities(i)) ' dB']);
    xticks([-(log(2)/log(frequencies(2)/frequencies(1))) 0 (log(2)/log(frequencies(2)/frequencies(1)))]);
    xticklabels({'-1','0','1'});
    xlabel('Octaves from BF');
    ylabel('Normalized FR');
    ylim([0 1]);
end
suptitle(['Changes in individual neuron tuning curves (' num2str(neuronCount) ' neurons)']);  
saveas(f2,fullfile(newDir,'Neuron-wise tuning curve changes.fig'));

f3 = figure;hold on;plot(mean(tcStacks{1,1}),'Color',[0 0 0]);plot(mean(tcStacks{1,2}),'Color',[0 0 1]);
set(f3,'Position',[25 200 1500 400])
for i = 1:uniqueIntensities
    subplot(1,uniqueIntensities,i);
    semilogx(frequencies/1000,mean(tcStacks{i,1}),'Color',[0 0 0],'LineWidth',2);
    hold on;
    semilogx(frequencies/1000,mean(tcStacks{i,2}),'Color',[0 0 1],'LineWidth',2);
    xlabel('Frequency (kHz)');
    ylabel('\Delta FR (Hz)');
    xlist = [1 5:5:uniqueFreq];
    xticks(round(frequencies(xlist)/1000,1));
    title(['SPL = ' num2str(intensities(i)) ' dB']);
    ylim([0 1]);
end
suptitle(['Population tuning curve (' num2str(neuronCount) ' neurons)']);  
saveas(f3,fullfile(newDir,'Population tuning curves.fig'));

f4 = figure;histogram(bestFrequencies,0.5:1:uniqueFreq+1);
ticks = xticks;ticks = ticks(2:end);
xticks([1 round(uniqueFreq/4):round(uniqueFreq/4):uniqueFreq]);
xticklabels(round(frequencies(xticks)/1000,1));
xlabel('Frequency (kHz)');
title('Histogram of best frequencies');
saveas(f4,fullfile(newDir,'Best frequency histogram.fig'));
