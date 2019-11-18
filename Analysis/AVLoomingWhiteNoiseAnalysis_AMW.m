
clear;

experiment = 'EP002';
mouseID = 'AW102';
session = 'Session1';
date = '20191117';
% stimPath = 'D:\Code\Aaron\StimInfo\20190621_noise_optoOffset_50rep_70dB_400k_005_exptInfo';
% stimPath = 'C:\Users\Aaron\Documents\MATLAB\Workspace\Analysis\EphyStimInfo\20190621_noise_optoOffset_50rep_70dB_400k_005_exptInfo';
stimPath = fullfile('E:\Electrophysiology\',experiment,mouseID,date,'StimInfo',[date '_' mouseID '_AVloomingWhiteNoise_stimInfo']);
analysisWindow = [-50 2000]; %ms relative to stimulus onset
quantWindow = [5 15]; %ms after stimulus onset
frBinWidth = 10; %ms
baselineWindow = [-10 0];

% dataFolder = fullfile('D:\KiloSort\',mouseID,session,folder,'SpikeMat');
dataFolder = fullfile('E:\Electrophysiology\',experiment,mouseID,date,'SpikeMat');
% dataFiles = dir(fullfile(dataFolder,'*noise_optoOffset*'));
dataFiles = dir(fullfile(dataFolder,'*AV_looming_whiteNoise*'));

% newDir = fullfile('D:\KiloSort\',mouseID,session,folder,'OptoNoiseResponses');
newDir = fullfile('E:\Electrophysiology\',experiment,mouseID,date,'AVLoomingWhiteNoiseResponses');
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

totalUnits = length(dataFiles);

binCount = analysisWindow(2)-analysisWindow(1)-frBinWidth+1; %sliding window
frScalar = 1000/frBinWidth;
% frScalarOLD = 1000/(quantWindow(2)-quantWindow(1));

soundResponsiveUnits = [];
lightResponsiveUnits = [];
audLoomResponsiveUnits = [];
visLoomResponsiveUnits = [];
interactUnits = [];

% rasterColorMap = 'cool';
% map = colormap(rasterColorMap);close;
% optoColors = map(round(linspace(1,length(map),uniqueEvents-2)),:);
% rasterColors = [0 0 0; 0 0 1; optoColors];

for n = 1:totalUnits
    n
    nData = load(fullfile(dataFiles(n).folder,dataFiles(n).name));
    
    if nData.CellInfo(6)==1
        eventTimes = nData.Events;
        spikeTimes = nData.SpikeData(1,:);
        
        spikeRaster = cell(uniqueEvents,2); %spike times (1), trial number (y-axis) (2)
        meanResponse = zeros(uniqueEvents,5); %index (1),  mean baseline (2), baseline std error (3), mean response (4), response std error (5),
        trialResponses = zeros(uniqueEvents,stimInfo.repeats);
        frTrain = zeros(uniqueEvents,binCount);
        
        for u = 1:uniqueEvents
            eventID = indices(u);
            eventsOfInterest = find(stimInfo.order==eventID);
            %         color = rasterColors(u,:);
            
            spikeFR = zeros(length(eventsOfInterest),binCount);
            rasterX = [];
            rasterY = [];
            
            prePost = zeros(length(eventsOfInterest),2);
            tempY = 0;
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
                %             prePost(e,2) = histcounts(trialSpikes,[quantWindow(1) quantWindow(2)]);
                prePost(e,2) = max(spikeFR(e,analysisWindow(1)*sign(analysisWindow(1))-frBinWidth+2:end)); %starts when the front end of the first bin is 1ms, not 0ms, that's why +2. Write out an example to convince yourself
            end
            
            spikeRaster{u,1} = rasterX;
            spikeRaster{u,2} = rasterY;
            meanResponse(u,1) = eventID;
            meanResponse(u,2) = mean(prePost(:,1))*frScalar;
            meanResponse(u,3) = std(prePost(:,1))*frScalar/sqrt(length(prePost(:,1)));
            meanResponse(u,4) = mean(prePost(:,2))*frScalar;
            meanResponse(u,5) = std(prePost(:,2))*frScalar/sqrt(length(prePost(:,2)));
            frTrain(u,:) = mean(spikeFR,1)*frScalar;
            meanResponse(u,4) = max(frTrain(u,analysisWindow(1)*sign(analysisWindow(1))-frBinWidth+2:end)); %this is the max of the mean, which should be more appropriate, as opposed to mean of max
        end
        
        unitData(n).raster = spikeRaster;
        unitData(n).meanResponse = meanResponse;
        unitData(n).frTrain = frTrain;
        unitData(n).type = nData.CellInfo(6);
        unitData(n).neuronNumber = nData.CellInfo(4);
        
        %     if h(1)
        %         soundResponsiveUnits = [soundResponsiveUnits nData.CellInfo(4)];
        % %         if nData.CellInfo(6)==1
        % %             soundResponsiveSingleUnits = [soundResponsiveSingleUnits nData.CellInfo(4)];
        % %         elseif nData.CellInfo(6)==2
        % %             soundResponsiveMultiUnits = [soundResponsiveMultiUnits nData.CellInfo(4)];
        % %         end
        %     end
        %     if h(2)
        %         optoResponsiveUnits = [optoResponsiveUnits nData.CellInfo(4)];
        % %         if nData.CellInfo(6)==1
        % %             optoResponsiveSingleUnits = [optoResponsiveSingleUnits nData.CellInfo(4)];
        % %         elseif nData.CellInfo(6)==2
        % %             optoResponsiveMultiUnits = [optoResponsiveMultiUnits nData.CellInfo(4)];
        % %         end
        %     end
        %     hInteract = ttest2(trialResponses(1,:),trialResponses(ind+2,:)); hInteract(isnan(hInteract))=0;
        %     if hInteract
        %         interactUnits = [interactUnits nData.CellInfo(4)];
        % %         if nData.CellInfo(6)==1
        % %             interactSingleUnits = [interactSingleUnits nData.CellInfo(6)];
        % %         elseif nData.CellInfo(6)==2
        % %             interactMultiUnits = [interactMultiUnits nData.CellInfo(4)];
        % %         end
        %     end
        %     if nData.CellInfo(6)==1
        %         singleUnits = [singleUnits nData.CellInfo(4)];
        %     elseif nData.CellInfo(6)==2
        %         multiUnits = [multiUnits nData.CellInfo(4)];
        %     end
        %
        f1 = figure;
        set(f1,'Position',[500 200 600 500])
        for sp = 1:uniqueEvents
            subplot(4,4,sp);
            scatter(spikeRaster{sp,1},spikeRaster{sp,2},'.');
            xlim([analysisWindow(1) analysisWindow(end)]);
        end
        suptitle(['Unit ' num2str(nData.CellInfo(4)) ' (n = ' num2str(n) '), Raster Plot, unitType = ' num2str(nData.CellInfo(6))]);
        
        f2 = figure;
        set(f2,'Position',[350 100 700 600]);
        modType = {'Sound'; 'Light'};
        direcType = {'app';'rec';'stat';'none'};
        for modality1=1:4
            for modality2 = 1:4
                subplot(2,4,modality1);hold on;
                audInd = find(stimInfo.index(:,2)==modality1);
                visInd = find(stimInfo.index(:,3)==modality2);
                ind = intersect(audInd,visInd);
                
                plot((analysisWindow(1)+frBinWidth):analysisWindow(2),frTrain(ind,:));
                xlim([analysisWindow(1) analysisWindow(end)]);
                title([modType{1} ' ' direcType{modality1}]);
                
                
                subplot(2,4,modality1+4);hold on;
                audInd = find(stimInfo.index(:,2)==modality2);
                visInd = find(stimInfo.index(:,3)==modality1);
                ind = intersect(audInd,visInd);
                
                plot((analysisWindow(1)+frBinWidth):analysisWindow(2),frTrain(ind,:));
                xlim([analysisWindow(1) analysisWindow(end)]);
                title([modType{2} ' ' direcType{modality1}]);
            end
            
        end
        legend('Approach','Recede','Static','None');
        suptitle(['Unit ' num2str(nData.CellInfo(4)) ' (n = ' num2str(n) '), Firing Rate, unitType = ' num2str(nData.CellInfo(6))]);
        
        f3 = figure;
        matrixified = reshape(meanResponse(:,4),[4 4]);
        imagesc(matrixified');
        xticks(1:4)
        xticklabels({'App','Rec','Static','None'});
        xlabel('Light');
        yticks(1:4)
        yticklabels({'App','Rec','Static','None'});
        ylabel('Sound');
        colorbar
        title(['Unit ' num2str(nData.CellInfo(4)) ' (n = ' num2str(n) '), Response Heat Map, unitType = ' num2str(nData.CellInfo(6))]);
        
        saveas(f1,fullfile(figureDir,['Unit ' num2str(nData.CellInfo(4)) ' Raster plots.fig']));
        saveas(f1,fullfile(figureDir,['Unit ' num2str(nData.CellInfo(4)) ' Raster plots.jpg']));
        saveas(f2,fullfile(figureDir,['Unit ' num2str(nData.CellInfo(4)) ' PSTH.fig']));
        saveas(f2,fullfile(figureDir,['Unit ' num2str(nData.CellInfo(4)) ' PSTH.jpg']));
        saveas(f3,fullfile(figureDir,['Unit ' num2str(nData.CellInfo(4)) ' Response heat map.fig']));
        saveas(f3,fullfile(figureDir,['Unit ' num2str(nData.CellInfo(4)) ' Response heat map.jpg']));
        close(f1);
        close(f2);
        close(f3);
    end
end

