
clear;

experiment = 'EP001';
mouseID = 'AW093';
session = 'Session2';
date = '20191104-2';
% stimPath = 'D:\Code\Aaron\StimInfo\20190621_noise_optoOffset_50rep_70dB_400k_005_stimInfo';
stimPath = 'C:\Users\Aaron\Documents\MATLAB\Workspace\Analysis\EphyStimInfo\20191103_noise_optoOffset_25msPulse_50rep_70dB_400k_005_stimInfo';
analysisWindow = [-50 200]; %ms relative to stimulus onset
quantWindow = [20 40]; %ms after stimulus onset
frBinWidth = 10; %ms
baselineWindow = [-10 0];

% dataFolder = fullfile('D:\KiloSort\',mouseID,session,folder,'SpikeMat');
dataFolder = fullfile('E:\Electrophysiology\',experiment,mouseID,date,'SpikeMat');
dataFiles = dir(fullfile(dataFolder,'*noise_optoOffset_25ms*'));

% newDir = fullfile('D:\KiloSort\',mouseID,session,folder,'OptoNoiseResponses');
newDir = fullfile('E:\Electrophysiology\',experiment,mouseID,date,'OptoNoiseResponses_25ms');
figureDir = fullfile(newDir,'Figures');
if ~exist(newDir)
    mkdir(newDir);
end
if ~exist(figureDir)
    mkdir(figureDir);
end

load(stimPath);
optoOffsets = stimInfo.optoOffsets;
uniqueEvents = size(stimInfo.index,1);
indices = stimInfo.index(:,1);

totalUnits = length(dataFiles);

binCount = analysisWindow(2)-analysisWindow(1)-frBinWidth+1; %sliding window
binScalar = 1000/frBinWidth;
binEdges = analysisWindow(1)+frBinWidth:1:analysisWindow(2);
frScalar = 1000/(quantWindow(2)-quantWindow(1));

soundResponsiveUnits = [];
optoResponsiveUnits = [];
interactUnits = [];
singleUnits = [];
multiUnits = [];
bestDuration = zeros(1,totalUnits);

rasterColorMap = 'cool';
map = colormap(rasterColorMap);close;
optoColors = map(round(linspace(1,length(map),uniqueEvents-2)),:);
rasterColors = [0 0 0; 0 0 1; optoColors];

for n = 1:totalUnits
    n
    nData = load(fullfile(dataFiles(n).folder,dataFiles(n).name));
    eventTimes = nData.Events;
    spikeTimes = nData.SpikeData(1,:);
    
    spikeRaster = cell(uniqueEvents,3); %spike times (1), trial number (y-axis) (2), color (3)
    meanResponse = zeros(uniqueEvents,5); %index (1),  mean baseline (2), baseline std error (3), mean response (4), response std error (5),
    trialResponses = zeros(uniqueEvents,stimInfo.repeats);
    frTrain = zeros(uniqueEvents,binCount);
    
    
    
    tempY = 0;
    for u = 1:uniqueEvents
        eventID = indices(u);
        eventsOfInterest = find(stimInfo.order==eventID);
        color = rasterColors(u,:);
        
        spikeFR = zeros(length(eventsOfInterest),binCount);
        rasterX = [];
        rasterY = [];
        
        prePost = zeros(length(eventsOfInterest),2);
        eventBaseline = min([0 stimInfo.index(eventID,4)]);
        for e = 1:length(eventsOfInterest)
            tempY = tempY+1;
            
            time = eventTimes(eventsOfInterest(e));
            alignedSpikes = (spikeTimes-time)/1000; %temporal resolution of cheetah is in microseconds, converting to milliseconds
            trialSpikes = alignedSpikes(find(alignedSpikes>analysisWindow(1) & alignedSpikes<analysisWindow(end)));
            rasterX = [rasterX trialSpikes];
            rasterY = [rasterY repmat(tempY,[1 length(trialSpikes)])];
            
            for b = 1:binCount
                timeMin = analysisWindow(1)+b-1;
                timeMax = timeMin+frBinWidth;
                spikeFR(e,b) = length(find(alignedSpikes>timeMin & alignedSpikes<timeMax));
            end
            prePost(e,1) = histcounts(trialSpikes,[baselineWindow(1) baselineWindow(2)]+eventBaseline);
            prePost(e,2) = histcounts(trialSpikes,[quantWindow(1) quantWindow(2)]);
        end
        
        spikeRaster{u,1} = rasterX;
        spikeRaster{u,2} = rasterY;
        spikeRaster{u,3} = color;
        meanResponse(u,1) = eventID;
        meanResponse(u,2) = mean(prePost(:,1))*frScalar;
        meanResponse(u,3) = std(prePost(:,1))*frScalar/sqrt(length(prePost(:,1)));
        meanResponse(u,4) = mean(prePost(:,2))*frScalar;
        meanResponse(u,5) = std(prePost(:,2))*frScalar/sqrt(length(prePost(:,2)));
        frTrain(u,:) = mean(spikeFR,1)*frScalar;
        
        [hh pp] = ttest2(prePost(:,1),prePost(:,2));
        h(u) = hh; h(isnan(h))=0;
        p(u) = pp; p(isnan(p))=1;
        
        trialResponses(u,:) = prePost(:,2);
        
    end
    
    unitData(n).raster = spikeRaster;
    unitData(n).meanResponse = meanResponse;
    unitData(n).frTrain = frTrain;
    unitData(n).type = nData.CellInfo(6);
    unitData(n).neuronNumber = nData.CellInfo(4);
    
    [val ind] = max(meanResponse(3:end,4));
    bestDuration(n) = ind;
    
    if h(1)
        soundResponsiveUnits = [soundResponsiveUnits nData.CellInfo(4)];
%         if nData.CellInfo(6)==1
%             soundResponsiveSingleUnits = [soundResponsiveSingleUnits nData.CellInfo(4)];
%         elseif nData.CellInfo(6)==2
%             soundResponsiveMultiUnits = [soundResponsiveMultiUnits nData.CellInfo(4)];
%         end
    end
    if h(2)
        optoResponsiveUnits = [optoResponsiveUnits nData.CellInfo(4)];
%         if nData.CellInfo(6)==1
%             optoResponsiveSingleUnits = [optoResponsiveSingleUnits nData.CellInfo(4)];
%         elseif nData.CellInfo(6)==2
%             optoResponsiveMultiUnits = [optoResponsiveMultiUnits nData.CellInfo(4)];
%         end
    end
    hInteract = ttest2(trialResponses(1,:),trialResponses(ind+2,:)); hInteract(isnan(hInteract))=0;
    if hInteract
        interactUnits = [interactUnits nData.CellInfo(4)];
%         if nData.CellInfo(6)==1
%             interactSingleUnits = [interactSingleUnits nData.CellInfo(6)];
%         elseif nData.CellInfo(6)==2
%             interactMultiUnits = [interactMultiUnits nData.CellInfo(4)];
%         end
    end
    if nData.CellInfo(6)==1
        singleUnits = [singleUnits nData.CellInfo(4)];
    elseif nData.CellInfo(6)==2
        multiUnits = [multiUnits nData.CellInfo(4)];
    end
    
    f1 = figure;
    subplot(1,3,1);
    for r = 1:size(spikeRaster,1)
        scatter(spikeRaster{r,1},spikeRaster{r,2},[],spikeRaster{r,3},'.');
        hold on;
    end
    xlim([analysisWindow(1) analysisWindow(end)]);
    xlabel('Time (ms)');
    ylabel('Trial type');
    
    subplot(1,3,2); hold on;
    for u = 1:size(frTrain,1)
        if u== 1 || u==2 || u==6
            plot(binEdges,frTrain(u,:),'Color',spikeRaster{u,3});
        end
    end
    xlabel('Times (ms)');
    ylabel('Firing rate (Hz)');
    
    subplot(1,3,3);
    data = [meanResponse(:,2) meanResponse(:,4)];
    bar(data);
    xlist = 1:1:size(meanResponse,1);
    xticks(xlist);
    xticklabels(['Sound' 'Opto' strsplit(num2str(optoOffsets))]);
    xlabel('Stimulus condition');
    ylabel('Firing rate (Hz)');
    %     title(['p(1)= ' num2str(p(1)) ', p(2) = ' num2str(p(2))]);
    
    suptitle(['Unit ' num2str(nData.CellInfo(4)) ' (n = ' num2str(n) '), unitType = ' num2str(nData.CellInfo(6))]);
    
    
    
    saveas(f1,fullfile(figureDir,['Unit ' num2str(nData.CellInfo(4)) ' opto response plots.fig']));
    saveas(f1,fullfile(figureDir,['Unit ' num2str(nData.CellInfo(4)) ' opto response plots.jpg']));
    close(f1);
end

f2 = figure;
subplot(1,3,1);hold on;
histogram(bestDuration(singleUnits));
histogram(bestDuration(multiUnits));
xticks(1:7);
xticklabels(optoOffsets);
xlabel('Laser offset relative to sound (ms)');
ylabel('Number of units');
legend('Single','Multi');
title('Laser offset with max FR');

subplot(1,3,2);
sound = zeros(1,length(interactUnits));
both = zeros(1,length(interactUnits));
delta = zeros(1,length(interactUnits));
list = zeros(1,length(interactUnits));
for i = 1:length(interactUnits)
    unitNumber = interactUnits(i);
    unitID = find([unitData.neuronNumber]==unitNumber);
    list(i) = unitID;
    sound(i) = unitData(unitID).meanResponse(1,4);
    both(i) = unitData(unitID).meanResponse(6,4);
    delta(i) = unitData(unitID).meanResponse(6,4) - unitData(unitID).meanResponse(1,4);
end
bar([mean(sound) mean(both)]);
xticklabels({'Sound','Sound+opto'});
ylabel('Firing rate (Hz)');
title('Comparison +/- opto');

subplot(1,3,3)
histogram(delta);
title('Difference +/- opto');
xlabel('\Delta FR (Hz)');
ylabel('Number of units');

saveas(f2,fullfile(newDir,'Summary figure'));



analysisParams.mouseID = mouseID;
analysisParams.session = session;
analysisParams.date = date;
analysisParams.stimPath = stimPath;
analysisParams.analysisWindow = analysisWindow; %ms
analysisParams.quantificationWindow = quantWindow;
analysisParams.frBinWidth = frBinWidth;

unitClassification.soundResponsiveUnits = soundResponsiveUnits;
unitClassification.optoResponsiveUnits = optoResponsiveUnits;
unitClassification.interactUnits = interactUnits;
unitClassification.totalUnits = totalUnits;
unitClassification.singleUnits = singleUnits;
unitClassification.multiUnits = multiUnits;
save(fullfile(newDir,'OptoNoiseResponseData.mat'),'unitData','analysisParams','unitClassification');

optoResponsiveUnits
interactUnits