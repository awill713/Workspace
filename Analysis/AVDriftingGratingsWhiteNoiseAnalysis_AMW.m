
clear;

experiment = 'EP004';
mouseID = 'AW118';
session = 'Session1';
date = '20200201';
stimPath = fullfile('E:\Electrophysiology\',experiment,mouseID,date,'StimInfo',[date '_' mouseID '_AVdriftingGratingsWhiteNoise_stimInfo']);

analysisWindow = [-100 1500]; %ms relative to stimulus onset
quantWindow = [50 100]; %ms relative to stimulus onset
baselineWindow = [-20 0]; %ms relative to stimulus onset

frBinWidth = 10; %ms

dataFolder = fullfile('E:\Electrophysiology\',experiment,mouseID,date,'SpikeMat');
dataFiles = dir(fullfile(dataFolder,'*AV_driftingGratings_whiteNoise*'));

% newDir = fullfile('D:\KiloSort\',mouseID,session,folder,'OptoNoiseResponses');
newDir = fullfile('E:\Electrophysiology\',experiment,mouseID,date,'AVDriftingGratingsWhiteNoise');
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
orientations = stimInfo.orientations;

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

lightResponsiveUnits = [];
orientationSelectiveUnits = [];
soundResponsiveUnits = [];
soundResponsiveFromANOVA = [];
interactUnits = [];
singleUnits = [];
multiUnits = [];

rasterColorMap = 'hsv';
map = colormap(rasterColorMap);close;
rasterColors = map(round(linspace(1,length(map),length(orientations))),:);

unNormalizedVis = [];
unNormalizedVisAud = [];
unNormalizedVis_Orientation = [];
unNormalizedVis_NoOrientation = [];
unNormalizedVisAud_Orientation = [];
unNormalizedVisAud_NoOrientation = [];
normalizedVis = [];
normalizedVisAud = [];
normalizedVis_Orientation = [];
normalizedVis_NoOrientation = [];
normalizedVisAud_Orientation = [];
normalizedVisAud_NoOrientation = [];

for n = 1:totalUnits
    n
    nData = load(fullfile(dataFiles(n).folder,dataFiles(n).name));
    
    eventTimes = nData.Events;
    spikeTimes = nData.SpikeData(1,:);
    
    spikeRaster = cell(uniqueEvents,3); %spike times (1), trial number (y-axis) (2), color (3)
    meanResponse = zeros(uniqueEvents,5); %index (1),  mean baseline (2), baseline std (3), mean response (4), response std (5)
    frTrain = zeros(uniqueEvents,binCount);
    frTrainTrials = zeros(uniqueEvents,repeats,binCount);
    quantTrials = zeros(uniqueEvents,repeats);
    
    for u = 1:uniqueEvents
        eventID = indices(u);
        eventsOfInterest = find(stimInfo.order==eventID);
        
        if u == uniqueEvents
            color = [0 0 0];
        else
            color = rasterColors(mod(u-1,size(rasterColors,1))+1,:);
        end
        
        spikeFR = zeros(length(eventsOfInterest),binCount);
        rasterX = [];
        rasterY = [];
        
        prePost = zeros(length(eventsOfInterest),2);
        tempY = mod(u-1,length(orientations))*repeats;
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
        
        frTrainTrials(u,:,:) = spikeFR*frScalar;
        frTrain(u,:) = mean(spikeFR,1)*frScalar;
        
        if u==uniqueEvents
            [hSound, pSound] = ttest2(prePost(:,1),prePost(:,2)); hSound(isnan(hSound)) = 0; pSound(isnan(pSound))= 1;
            if hSound
                soundResponsiveUnits = [soundResponsiveUnits nData.CellInfo(4)];
            end
        end
    end
    
    unitData(n).raster = spikeRaster;
    unitData(n).meanResponse = meanResponse;
    unitData(n).frTrain = frTrain;
    unitData(n).frTrainTrials = frTrainTrials;
    unitData(n).type = nData.CellInfo(6);
    unitData(n).neuronNumber = nData.CellInfo(4);
    neuronNumber = nData.CellInfo(4);
    
    if nData.CellInfo(6)==1
        singleUnits = [singleUnits neuronNumber];
    elseif nData.CellInfo(6)==2
        multiUnits = [multiUnits neuronNumber];
    end
    
    statsPre = mean(frTrainTrials(1:length(orientations),:,baselineFirstBin:baselineLastBin),3);
    statsPre = statsPre(:);
    statsPost = mean(frTrainTrials(1:length(orientations),:,quantFirstBin:quantLastBin),3);
    statsPost = statsPost(:);
    [hLight pLight] = ttest2(statsPre,statsPost); hLight(isnan(hLight)) = 0; pLight(isnan(pLight))= 1;
    if hLight
        lightResponsiveUnits = [lightResponsiveUnits neuronNumber];
    end
    
    anovaMatrix = zeros(repeats*2,length(orientations));
    for i = 1:uniqueEvents-1
        row = floor((i-1)/length(orientations))*repeats+1;
        column = mod(i-1,length(orientations))+1;
        anovaMatrix(row:row+repeats-1,column) = quantTrials(i,:);
    end
    [pOrientation, table, stats] = anova2(anovaMatrix,repeats,'off'); pOrientation(isnan(pOrientation)) = 0;
    if pOrientation(2)<0.05
        soundResponsiveFromANOVA = [soundResponsiveFromANOVA neuronNumber];
    end
    if pOrientation(1)<0.05
        orientationSelectiveUnits = [orientationSelectiveUnits neuronNumber];
    end
    if pOrientation(3)<0.05
        interactUnits = [interactUnits neuronNumber];
    end
    
    f1 = figure;
    set(f1,'Position',[150 80 1250 700]);
    subplot(1,4,1);
    for r = 1:length(orientations)
        scatter(spikeRaster{r,1},spikeRaster{r,2},[],spikeRaster{r,3},'.');
        hold on;
    end
    xlim([analysisWindow(1) analysisWindow(end)]);
    xlabel('Time (ms)');
    ylabel('Trial type (orientation)');
    ylim([-repeats repeats*length(orientations)]);
    yticks([]);
    title('Visual');
    legend(cellstr(num2str(orientations')));
    %         c = colorbar;
    %         set(c,'YTick',[0 1]);
    %         set(c,'YTickLabel',[intensities(1,1) intensities(1,end)])
    %         ylabel(c,'Intensity (dB)');
    
    subplot(1,4,2);
    for r = length(orientations)+1 : length(orientations)*2
        scatter(spikeRaster{r,1},spikeRaster{r,2},[],spikeRaster{r,3},'.');
        hold on;
    end
    scatter(spikeRaster{uniqueEvents,1},spikeRaster{uniqueEvents,2}-repeats,[],spikeRaster{uniqueEvents,3},'.');
    xlim([analysisWindow(1) analysisWindow(end)]);
    xlabel('Time (ms)');
    ylabel('Trial type (orientation)');
    ylim([-repeats repeats*length(orientations)]);
    yticks([]);
    title('Visual/Audio');
    %         c = colorbar;
    %         colormap(rasterColorMap);
    %         set(c,'YTick',[0 1]);
    %         set(c,'YTickLabel',[intensities(1,1) intensities(1,end)])
    %         ylabel(c,'Intensity (dB)');
    
    subplot(1,4,3);hold on;
    psthScalar = 20/max(max(frTrain));
    yMin = 0;
    for i = 1:length(orientations)
        trace = (frTrain(i,:)-mean(frTrain(i,baselineFirstBin:baselineLastBin)))*psthScalar + 20*(i-1);
        plot(binEdges,trace,'Color',spikeRaster{i,3});
        yMin = min([yMin min(trace)]);
        
        trace = (frTrain(i+length(orientations),:)-mean(frTrain(i+length(orientations),baselineFirstBin:baselineLastBin)))*psthScalar + 20*(i-1);
        plot(binEdges,trace,'Color',spikeRaster{i,3},'LineStyle','--');
        yMin = min([yMin min(trace)]);
    end
    trace = (frTrain(uniqueEvents,:)-mean(frTrain(uniqueEvents,baselineFirstBin:baselineLastBin)))*psthScalar - 20;
    plot(binEdges,trace,'Color',spikeRaster{uniqueEvents,3});
    yMin = min([yMin min(trace)]);
    
    yLimits = ylim();
    ylim([yMin-5  yLimits(2)]);
    yticks([]);
    xlim([analysisWindow]);
    xlabel('Time (ms)');
    
    subplot(1,4,4);
    polarplot(2*pi*orientations./360,meanResponse(1:length(orientations),4),'Color',[0 0 0],'LineWidth',2);
    hold on;
    polarplot(2*pi*orientations./360,meanResponse(length(orientations)+1:2*length(orientations),4),'Color',[0 0 1],'LineWidth',2);
    %         plot(meanResponse(1:length(intensities),4),'Color',[0 0 0],'LineWidth',2);
    %         plot(meanResponse(length(intensities)+1:2*(length(intensities)),4),'Color',[0 0 1],'LineWidth',2);
    %         errorbar(meanResponse(1:length(intensities),4),meanResponse(1:length(intensities),5),'Color',[0 0 0],'LineWidth',2);
    %         errorbar(meanResponse(length(intensities)+1:2*(length(intensities)),4),mean(length(intensities)+1:2*length(intensities),5),'Color',[0 0 1],'LineWidth',2);
    %         xticks(1:length(intensities));
    %         xticklabels(intensities)
    %         xlabel('Sound intensity (dB)');
    %         legend({'Vis','Vis/Aud'})
    
    suptitle({['Unit ' num2str(nData.CellInfo(4)) ' (n = ' num2str(n)...
        '), unitType = ' num2str(nData.CellInfo(6))] ['Light responsive = ' num2str(pLight)...
        ', Sound responsive = ' num2str(pSound) ', Orientation selective = ' num2str(pOrientation(1))...
        ', Sound-modulated = ' num2str(pOrientation(2)) ', ' num2str(pOrientation(3))]});
    
    
            saveas(f1,fullfile(figureDir,['Unit ' num2str(nData.CellInfo(4)) ' response plots.fig']));
            saveas(f1,fullfile(figureDir,['Unit ' num2str(nData.CellInfo(4)) ' response plots.jpg']));
            close(f1);
    
    if nData.CellInfo(6)~=0 && hLight
        responseVis = meanResponse(1:length(orientations),4);
        responseVisAud = meanResponse(length(orientations)+1:2*length(orientations),4);
        
        [peak ind] = max(responseVis);
        targetIndex = ceil(length(orientations)/2) + 1;
        
        realignedVis = circshift(responseVis,targetIndex - ind);
        realignedVis = [realignedVis; realignedVis(1)];
        realignedVisAud = circshift(responseVisAud,targetIndex - ind);
        realignedVisAud = [realignedVisAud; realignedVisAud(1)];
        
        unNormalizedVis = [unNormalizedVis; realignedVis'];
        unNormalizedVisAud = [unNormalizedVisAud; realignedVisAud'];
        
        respMin = min(realignedVis); respMax = max(realignedVis);
        normVis = (realignedVis - respMin) / (respMax - respMin);
        normVisAud = (realignedVisAud - respMin) / (respMax - respMin);
        
        normalizedVis = [normalizedVis; normVis'];
        normalizedVisAud = [normalizedVisAud; normVisAud'];
    end
end

for n = 1:totalNeurons
    if unitData(n).type~=0 && ismember(unitData(n).neuronNumber,responsiveUnits.lightResponsiveUnits)
        responseVis = unitData(n).meanResponse(1:length(orientations),4);
        responseVisAud = unitData(n).meanResponse(length(orientations)+1:2*length(orientations),4);
        
        [peak ind] = max(responseVis);
        targetIndex = ceil(length(orientations)/2) + 1;
        
        realignedVis = circshift(responseVis,targetIndex - ind);
        realignedVis = [realignedVis; realignedVis(1)];
        realignedVisAud = circshift(responseVisAud,targetIndex - ind);
        realignedVisAud = [realignedVisAud; realignedVisAud(1)];
        
        unNormalizedVis = [unNormalizedVis; realignedVis'];
        unNormalizedVisAud = [unNormalizedVisAud; realignedVisAud'];
        
        respMin = min(realignedVis); respMax = max(realignedVis);
        normVis = (realignedVis - respMin) / (respMax - respMin);
        normVisAud = (realignedVisAud - respMin) / (respMax - respMin);
        
        normalizedVis = [normalizedVis; normVis'];
        normalizedVisAud = [normalizedVisAud; normVisAud'];
        
        if ismember(unitData(n).neuronNumber,responsiveUnits.orientationSelectiveUnits)
            unNormalizedVis_Orientation = [unNormalizedVis_Orientation; realignedVis'];
            unNormalizedVisAud_Orientation = [unNormalizedVisAud_Orientation; realignedVisAud'];
            normalizedVis_Orientation = [normalizedVis_Orientation; normVis'];
            normalizedVisAud_Orientation = [normalizedVisAud_Orientation; normVisAud'];
        else
            unNormalizedVis_NoOrientation = [unNormalizedVis_NoOrientation; realignedVis'];
            unNormalizedVisAud_NoOrientation = [unNormalizedVisAud_NoOrientation; realignedVisAud'];
            normalizedVis_NoOrientation = [normalizedVis_NoOrientation; normVis'];
            normalizedVisAud_NoOrientation = [normalizedVisAud_NoOrientation; normVisAud'];
        end
    end
end

f2 = figure;
set(f2,'Position',[150 80 1250 700]);
subplot(2,3,1);hold on;
plot(linspace(-180,180,length(orientations)+1),mean(unNormalizedVis),'Color',[0 0 0],'LineWidth',2);
plot(linspace(-180,180,length(orientations)+1),mean(unNormalizedVisAud),'Color',[0 0 1],'LineWidth',2);
ylabel('Firing rate (Hz)');
xlabel('\Delta orientation (degrees)');
xticks(-180:90:180);
title('All light responsive (unnormalized)');

subplot(2,3,2);hold on;
plot(linspace(-180,180,length(orientations)+1),mean(unNormalizedVis_Orientation),'Color',[0 0 0],'LineWidth',2);
plot(linspace(-180,180,length(orientations)+1),mean(unNormalizedVisAud_Orientation),'Color',[0 0 1],'LineWidth',2);
ylabel('Firing rate (Hz)');
xlabel('\Delta orientation (degrees)');
xticks(-180:90:180);
title('Orientation selective (unnormalized)');

subplot(2,3,3);hold on;
plot(linspace(-180,180,length(orientations)+1),mean(unNormalizedVis_NoOrientation),'Color',[0 0 0],'LineWidth',2);
plot(linspace(-180,180,length(orientations)+1),mean(unNormalizedVisAud_NoOrientation),'Color',[0 0 1],'LineWidth',2);
ylabel('Firing rate (Hz)');
xlabel('\Delta orientation (degrees)');
xticks(-180:90:180);
title('Orientation UNselective (unnormalized)');

subplot(2,3,4);hold on;
plot(linspace(-180,180,length(orientations)+1),mean(normalizedVis),'Color',[0 0 0],'LineWidth',2);
plot(linspace(-180,180,length(orientations)+1),mean(normalizedVisAud),'Color',[0 0 1],'LineWidth',2);
title('All light responsive (normalized)');
ylabel('Normalized firing rate');
xlabel('\Delta orientation (degrees)');
xticks(-180:90:180);

subplot(2,3,5);hold on;
plot(linspace(-180,180,length(orientations)+1),mean(normalizedVis_Orientation),'Color',[0 0 0],'LineWidth',2);
plot(linspace(-180,180,length(orientations)+1),mean(normalizedVisAud_Orientation),'Color',[0 0 1],'LineWidth',2);
ylabel('Normalized firing rate');
xlabel('\Delta orientation (degrees)');
xticks(-180:90:180);
title('Orientation selective (normalized)');

subplot(2,3,6);hold on;
plot(linspace(-180,180,length(orientations)+1),mean(normalizedVis_NoOrientation),'Color',[0 0 0],'LineWidth',2);
plot(linspace(-180,180,length(orientations)+1),mean(normalizedVisAud_NoOrientation),'Color',[0 0 1],'LineWidth',2);
ylabel('Normalized firing rate');
xlabel('\Delta orientation (degrees)');
xticks(-180:90:180);
title('Orientation UNselective (normalized)');
legend({'Vis','Vis/Aud'});

suptitle({'Population orientation preference, +/- white noise'});
saveas(f2,fullfile(newDir,'Population orientation preference'));


analysisParams.mouseID = mouseID;
analysisParams.session = session;
analysisParams.date = date;
analysisParams.stimPath = stimPath;
analysisParams.analysisWindow = analysisWindow;
analysisParams.frBinWidth = frBinWidth;
analysisParams.quantWindow = quantWindow;
analysisParams.baselineWindow = baselineWindow;

responsiveUnits.soundResponsiveUnits = soundResponsiveUnits;
responsiveUnits.lightResponsiveUnits = lightResponsiveUnits;
responsiveUnits.interactUnits = interactUnits;
responsiveUnits.orientationSelectiveUnits = orientationSelectiveUnits;
responsiveUnits.soundResponsiveFromANOVA = soundResponsiveFromANOVA;


save(fullfile(newDir,'AVStaticMultiAmpNoiseData.mat'),'unitData','analysisParams','responsiveUnits');
responsiveUnits


