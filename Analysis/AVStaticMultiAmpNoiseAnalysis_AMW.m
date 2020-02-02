
clear;

experiment = 'EP003';
mouseID = 'AW113';
session = 'Session2';
date = '20200120';
stimPath = fullfile('E:\Electrophysiology\',experiment,mouseID,date,'StimInfo',[date '_' mouseID '_AVstaticMultiAmpNoise_stimInfo']);

analysisWindow = [-50 200]; %ms relative to stimulus onset
quantWindow = [10 40]; %ms relative to stimulus onset
baselineWindow = [-10 0]; %ms relative to stimulus onset

frBinWidth = 10; %ms

dataFolder = fullfile('E:\Electrophysiology\',experiment,mouseID,date,'SpikeMat');
dataFiles = dir(fullfile(dataFolder,'*AV_static_multiAmpNoise*'));

% newDir = fullfile('D:\KiloSort\',mouseID,session,folder,'OptoNoiseResponses');
newDir = fullfile('E:\Electrophysiology\',experiment,mouseID,date,'AVStaticMultiAmpNoise');
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
intensities = stimInfo.intensities;

totalUnits = length(dataFiles);

binCount = analysisWindow(2)-analysisWindow(1)-frBinWidth+1; %sliding window
frScalar = 1000/frBinWidth;
quantScalar = 1000/(quantWindow(2)-quantWindow(1));
quantFirstBin = quantWindow(1)-analysisWindow(1)-frBinWidth+1;
quantLastBin = quantWindow(2)-analysisWindow(1)-frBinWidth+1;
baselineFirstBin = baselineWindow(1) - analysisWindow(1)-frBinWidth+1;
baselineLastBin = baselineWindow(2) - analysisWindow(1)-frBinWidth+1;
binEdges = analysisWindow(1)+frBinWidth:1:analysisWindow(2);

soundResponsiveUnit = [];
lightResponsiveUnit = [];
interactUnit = [];
singleUnits = [];
multiUnits = [];

rasterColorMap = 'cool';
map = colormap(rasterColorMap);close;
optoColors = map(round(linspace(1,length(map),length(intensities)-1)),:);
rasterColors = [0 0 0; optoColors];

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
            
            color = rasterColors(mod(u-1,size(rasterColors,1))+1,:);
            
            spikeFR = zeros(length(eventsOfInterest),binCount);
            rasterX = [];
            rasterY = [];
            
            prePost = zeros(length(eventsOfInterest),2);
            tempY = mod(u-1,length(intensities))*repeats;
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
            meanResponse(u,2) = mean(prePost(:,1))*quantScalar;
            meanResponse(u,3) = std(prePost(:,1))*quantScalar/sqrt(repeats);
            meanResponse(u,4) = mean(prePost(:,2))*quantScalar;
            meanResponse(u,5) = std(prePost(:,2))*quantScalar/sqrt(repeats);
            
            frTrainTrials(u,:,:) = spikeFR*frScalar;
            frTrain(u,:) = mean(spikeFR,1)*frScalar;
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
        
        
        anovaMatrix = zeros(repeats*2,length(intensities));
        for i = 1:uniqueEvents
            row = floor((i-1)/length(intensities))*repeats+1;
            column = mod(i-1,length(intensities))+1;
            anovaMatrix(row:row+repeats-1,column) = quantTrials(i,:);
        end
        [p, table, stats] = anova2(anovaMatrix,repeats,'off');
        if p(1)<0.05
            soundResponsiveUnit = [soundResponsiveUnit neuronNumber];
        end
        if p(2)<0.05
            lightResponsiveUnit = [lightResponsiveUnit neuronNumber];
        end
        if p(3)<0.05
            interactUnit = [interactUnit neuronNumber];
        end
        
        f1 = figure;
        set(f1,'Position',[150 80 1250 700]);
        subplot(1,4,1);
        for r = 1:length(intensities)
            scatter(spikeRaster{r,1},spikeRaster{r,2},[],spikeRaster{r,3},'.');
            hold on;
        end
        xlim([analysisWindow(1) analysisWindow(end)]);
        xlabel('Time (ms)');
        ylabel('Trial type');
        ylim([0 repeats*length(intensities)]);
        title('Sound');
%         c = colorbar;
%         set(c,'YTick',[0 1]);
%         set(c,'YTickLabel',[intensities(1,1) intensities(1,end)])
%         ylabel(c,'Intensity (dB)');
        
        subplot(1,4,2);
        for r = length(intensities)+1 : length(intensities)*2
            scatter(spikeRaster{r,1},spikeRaster{r,2},[],spikeRaster{r,3},'.');
            hold on;
        end
        xlim([analysisWindow(1) analysisWindow(end)]);
        xlabel('Time (ms)');
        ylabel('Trial type');
        ylim([0 repeats*length(intensities)]);
        title('Sound/Light');
%         c = colorbar;
%         colormap(rasterColorMap);
%         set(c,'YTick',[0 1]);
%         set(c,'YTickLabel',[intensities(1,1) intensities(1,end)])
%         ylabel(c,'Intensity (dB)');
        
        subplot(1,4,3);hold on;
        psthScalar = 20/max(max(frTrain));
        for i = 1:length(intensities)
            trace = (frTrain(i,:)-mean(frTrain(i,baselineFirstBin:baselineLastBin)))*psthScalar + 20*(i-1);
            plot(binEdges,trace,'Color',spikeRaster{i,3});
            
            trace = (frTrain(i+length(intensities),:)-mean(frTrain(i+length(intensities),baselineFirstBin:baselineLastBin)))*psthScalar + 20*(i-1);
            plot(binEdges,trace,'Color',spikeRaster{i,3},'LineStyle','--');
        end
        yticks([]);
        xlim([analysisWindow]);
        xlabel('Time (ms)');
        
        subplot(1,4,4);hold on;
%         plot(meanResponse(1:length(intensities),4),'Color',[0 0 0],'LineWidth',2);
%         plot(meanResponse(length(intensities)+1:2*(length(intensities)),4),'Color',[0 0 1],'LineWidth',2);
        errorbar(meanResponse(1:length(intensities),4),meanResponse(1:length(intensities),5),'Color',[0 0 0],'LineWidth',2);
        errorbar(meanResponse(length(intensities)+1:2*(length(intensities)),4),mean(length(intensities)+1:2*length(intensities),5),'Color',[0 0 1],'LineWidth',2);
        xticks(1:length(intensities));
        xticklabels(intensities)
        xlabel('Sound intensity (dB)');
        legend({'Sound','Sound/Light'})
        
        suptitle(['Unit ' num2str(nData.CellInfo(4)) ' (n = ' num2str(n) '), unitType = ' num2str(nData.CellInfo(6)) ', p = ' num2str(p) ]);

        
        saveas(f1,fullfile(figureDir,['Unit ' num2str(nData.CellInfo(4)) ' response plots.fig']));
        saveas(f1,fullfile(figureDir,['Unit ' num2str(nData.CellInfo(4)) ' response plots.jpg']));
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

responsiveUnits.soundResponsiveUnits = soundResponsiveUnit;
responsiveUnits.lightResponsiveUnit = lightResponsiveUnit;
responsiveUnits.interactUnit = interactUnit;


save(fullfile(newDir,'AVStaticMultiAmpNoiseData.mat'),'unitData','analysisParams','responsiveUnits');
responsiveUnits


