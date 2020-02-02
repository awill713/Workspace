
clear;

experiment = 'EP003';
mouseID = 'AW112';
session = 'Session2';
date = '20200120';
% stimPath = 'D:\Code\Aaron\StimInfo\20190621_noise_optoOffset_50rep_70dB_400k_005_exptInfo';
% stimPath = 'C:\Users\Aaron\Documents\MATLAB\Workspace\Analysis\EphyStimInfo\20190621_noise_optoOffset_50rep_70dB_400k_005_exptInfo';
stimPath = fullfile('E:\Electrophysiology\',experiment,mouseID,date,'StimInfo',[date '_' mouseID '_AVloomingWhiteNoise_stimInfo']);
analysisWindow = [-50 2000]; %ms relative to stimulus onset
peakQuantWindow = [0 100]; %ms after stimulus onset
sustQuantWindow = [100 25]; % [ms after stimulus onset, ms after stimulus offset]
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
repeats = stimInfo.repeats;

totalUnits = length(dataFiles);

binCount = analysisWindow(2)-analysisWindow(1)-frBinWidth+1; %sliding window
frScalar = 1000/frBinWidth;
peakFirstBin = peakQuantWindow(1)-analysisWindow(1)-frBinWidth+1;
peakLastBin = peakQuantWindow(2)-analysisWindow(1)-frBinWidth+1;
quantFirstBin = sustQuantWindow(1)-analysisWindow(1)-frBinWidth+1;
quantLastBin = sustQuantWindow(2)+(1000*stimInfo.stimDuration)-analysisWindow(1)-frBinWidth+1;

peakSoundResponsiveUnit = [];
peakLightResponsiveUnit = [];
peakSoundTunedUnit = [];
peakLightTunedUnit = [];
peakInteractUnit = [];
sustSoundResponsiveUnit = [];
sustLightResponsiveUnit = [];
sustSoundTunedUnit = [];
sustLightTunedUnit = [];
sustInteractUnit = [];

% rasterColorMap = 'cool';
% map = colormap(rasterColorMap);close;
% optoColors = map(round(linspace(1,length(map),uniqueEvents-2)),:);
% rasterColors = [0 0 0; 0 0 1; optoColors];

for n = 1:totalUnits
    n
    nData = load(fullfile(dataFiles(n).folder,dataFiles(n).name));
    
%     if nData.CellInfo(6)==1
        eventTimes = nData.Events;
        spikeTimes = nData.SpikeData(1,:);
        
        spikeRaster = cell(uniqueEvents,2); %spike times (1), trial number (y-axis) (2)
%         meanResponse = zeros(uniqueEvents,5); %index (1),  mean baseline (2), baseline std error (3), mean response (4), response std error (5),
        meanSustResponse = zeros(uniqueEvents,3);%index (1), mean response (2), response std (3)
        meanPeakResponse = zeros(uniqueEvents,3);%index (1),  mean response (2), response std (3)
        frTrain = zeros(uniqueEvents,binCount);
        
        peakTrialResponses = zeros(uniqueEvents,repeats);
        sustTrialResponses = zeros(uniqueEvents,repeats);
        
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
                sustTrialResponses(u,e) = mean(spikeFR(e,quantFirstBin:quantLastBin))*frScalar;
                peakTrialResponses(u,e) = mean(spikeFR(e,peakFirstBin:peakLastBin))*frScalar;
            end
            
            spikeRaster{u,1} = rasterX;
            spikeRaster{u,2} = rasterY;
            meanSustResponse(u,1) = eventID;
            meanSustResponse(u,2) = mean(sustTrialResponses(u,:));
            meanSustResponse(u,3) = std(sustTrialResponses(u,:));
            meanPeakResponse(u,1) = eventID;
            meanPeakResponse(u,2) = mean(peakTrialResponses(u,:));
            meanPeakResponse(u,3) = std(sustTrialResponses(u,:));
            
            frTrain(u,:) = mean(spikeFR,1)*frScalar;
            
%             meanResponse(u,1) = eventID;
%             meanResponse(u,2) = mean(prePost(:,1))*frScalar;
%             meanResponse(u,3) = std(prePost(:,1))*frScalar/sqrt(length(prePost(:,1)));
%             meanResponse(u,4) = mean(prePost(:,2))*frScalar;
%             meanResponse(u,5) = std(prePost(:,2))*frScalar/sqrt(length(prePost(:,2)));
%             meanResponse(u,4) = max(frTrain(u,analysisWindow(1)*sign(analysisWindow(1))-frBinWidth+2:end)); %this is the max of the mean, which should be more appropriate, as opposed to mean of max
        end
        
        unitData(n).raster = spikeRaster;
        unitData(n).sustResponse = meanSustResponse;
        unitData(n).peakResponse = meanPeakResponse;
        unitData(n).frTrain = frTrain;
        unitData(n).type = nData.CellInfo(6);
        unitData(n).neuronNumber = nData.CellInfo(4);
        neuronNumber = nData.CellInfo(4);
        
       
        %sustained response
        anovaMatrix = zeros(4*repeats,4);
        for i = 1:uniqueEvents
            row = (ceil(i/4)-1)*repeats+1;
            column = mod(i-1,4)+1;
            anovaMatrix(row:row+(repeats-1),column) = sustTrialResponses(i,:);
        end
        [p, table, stats] = anova2(anovaMatrix,repeats,'off');
        if p(2)<0.05
            comparison = multcompare(stats,'estimate','row','display','off');
            if comparison(3,6)<0.05 || comparison(5,6)<0.05 || comparison(6,6)<0.05
                sustSoundResponsiveUnit = [sustSoundResponsiveUnit neuronNumber];
            end
            if comparison(1,6)<0.05 || comparison(2,6)<0.05 || comparison(4,6)<0.05
                sustSoundTunedUnit = [sustSoundTunedUnit neuronNumber];
            end
        end
        if p(1)<0.05
            comparison = multcompare(stats,'estimate','column','display','off');
            if comparison(3,6)<0.05 || comparison(5,6)<0.05 || comparison(6,6)<0.05
                sustLightResponsiveUnit = [sustLightResponsiveUnit neuronNumber];
            end
            if comparison(1,6)<0.05 || comparison(2,6)<0.05 || comparison(4,6)<0.05
                sustLightTunedUnit = [sustLightTunedUnit neuronNumber];
            end
        end
        if p(3)<0.05
            sustInteractUnit = [sustInteractUnit neuronNumber];
        end
        
        %peak response
        anovaMatrix = zeros(4*repeats,4);
        for i = 1:uniqueEvents
            row = (ceil(i/4)-1)*repeats+1;
            column = mod(i-1,4)+1;
            anovaMatrix(row:row+(repeats-1),column) = peakTrialResponses(i,:);
        end
        [p, table, stats] = anova2(anovaMatrix,repeats,'off');
        if p(2)<0.05
            comparison = multcompare(stats,'estimate','row','display','off');
            if comparison(3,6)<0.05 || comparison(5,6)<0.05 || comparison(6,6)<0.05
                peakSoundResponsiveUnit = [peakSoundResponsiveUnit neuronNumber];
            end
            if comparison(1,6)<0.05 || comparison(2,6)<0.05 || comparison(4,6)<0.05
                peakSoundTunedUnit = [peakSoundTunedUnit neuronNumber];
            end
        end
        if p(1)<0.05
            comparison = multcompare(stats,'estimate','column','display','off');
            if comparison(3,6)<0.05 || comparison(5,6)<0.05 || comparison(6,6)<0.05
                peakLightResponsiveUnit = [peakLightResponsiveUnit neuronNumber];
            end
            if comparison(1,6)<0.05 || comparison(2,6)<0.05 || comparison(4,6)<0.05
                peakLightTunedUnit = [peakLightTunedUnit neuronNumber];
            end
        end
        if p(3)<0.05
            peakInteractUnit = [peakInteractUnit neuronNumber];
        end
        
        
        f1 = figure;
        set(f1,'Position',[150 80 1250 700]);
        for sp = 1:uniqueEvents
            subplot(4,4,sp);
            scatter(spikeRaster{sp,1},spikeRaster{sp,2},'.');
            xlim([analysisWindow(1) analysisWindow(end)]);
        end
        suptitle(['Unit ' num2str(nData.CellInfo(4)) ' (n = ' num2str(n) '), Raster Plot, unitType = ' num2str(nData.CellInfo(6))]);
        
        f2 = figure;
        set(f2,'Position',[150 80 1250 700]);
        modType = {'Sound'; 'Light'};
        direcType = {'app';'rec';'stat';'none'};
        for modality1=1:4
            for modality2 = 1:4
                subplot(2,4,modality1);hold on;
                audInd = find(stimInfo.index(:,2)==modality1);
                visInd = find(stimInfo.index(:,3)==modality2);
                ind = intersect(audInd,visInd);
                
                plot((analysisWindow(1)+frBinWidth):analysisWindow(2),smooth(frTrain(ind,:),50));
                xlim([analysisWindow(1) analysisWindow(end)]);
                title([modType{1} ' ' direcType{modality1}]);
                
                
                subplot(2,4,modality1+4);hold on;
                audInd = find(stimInfo.index(:,2)==modality2);
                visInd = find(stimInfo.index(:,3)==modality1);
                ind = intersect(audInd,visInd);
                
                plot((analysisWindow(1)+frBinWidth):analysisWindow(2),smooth(frTrain(ind,:),50));
                xlim([analysisWindow(1) analysisWindow(end)]);
                title([modType{2} ' ' direcType{modality1}]);
            end
            
        end
        legend('Approach','Recede','Static','None');
        suptitle(['Unit ' num2str(nData.CellInfo(4)) ' (n = ' num2str(n) '), Firing Rate, unitType = ' num2str(nData.CellInfo(6))]);
        
        f3 = figure;
        set(f3,'Position',[150 150 1250 500]);
        matrixPeak = reshape(meanPeakResponse(:,2),[4 4]);
        subplot(1,2,1);
        imagesc(matrixPeak');
        xticks(1:4);xticklabels({'App','Rec','Static','None'});xlabel('Light');
        yticks(1:4);yticklabels({'App','Rec','Static','None'});ylabel('Sound');
        colorbar
        title('Peak response');
        matrixSust = reshape(meanSustResponse(:,2),[4 4]);
        subplot(1,2,2);
        imagesc(matrixSust');
        xticks(1:4);xticklabels({'App','Rec','Static','None'});xlabel('Light');
        yticks(1:4);yticklabels({'App','Rec','Static','None'});ylabel('Sound');
        colorbar
        title('Sustained response');
        suptitle(['Unit ' num2str(nData.CellInfo(4)) ' (n = ' num2str(n) '), Response Heat Map, unitType = ' num2str(nData.CellInfo(6))]);
        
        saveas(f1,fullfile(figureDir,['Unit ' num2str(nData.CellInfo(4)) ' Raster plots.fig']));
        saveas(f1,fullfile(figureDir,['Unit ' num2str(nData.CellInfo(4)) ' Raster plots.jpg']));
        saveas(f2,fullfile(figureDir,['Unit ' num2str(nData.CellInfo(4)) ' PSTH.fig']));
        saveas(f2,fullfile(figureDir,['Unit ' num2str(nData.CellInfo(4)) ' PSTH.jpg']));
        saveas(f3,fullfile(figureDir,['Unit ' num2str(nData.CellInfo(4)) ' Response heat map.fig']));
        saveas(f3,fullfile(figureDir,['Unit ' num2str(nData.CellInfo(4)) ' Response heat map.jpg']));
        close(f1);
        close(f2);
        close(f3);
%     end
end


analysisParams.mouseID = mouseID;
analysisParams.session = session;
analysisParams.date = date;
analysisParams.stimPath = stimPath;
analysisParams.analysisWindow = analysisWindow;
analysisParams.peakQuantWindow = peakQuantWindow;
analysisParams.sustQuantWindow = sustQuantWindow;
analysisParams.frBinWidth = frBinWidth;
analysisParams.baselineWindow = baselineWindow;

responsiveUnits.peakSoundResponsiveUnit = peakSoundResponsiveUnit;
responsiveUnits.peakLightResponsiveUnit = peakLightResponsiveUnit;
responsiveUnits.peakSoundTunedUnit = peakSoundTunedUnit;
responsiveUnits.peakLightTunedUnit = peakLightTunedUnit;
responsiveUnits.peakInteractUnit = peakInteractUnit;
responsiveUnits.sustSoundResponsiveUnit = sustSoundResponsiveUnit;
responsiveUnits.sustLightResponsiveUnit = sustLightResponsiveUnit;
responsiveUnits.sustSoundTunedUnit = sustSoundTunedUnit;
responsiveUnits.sustLightTunedUnit = sustLightTunedUnit;
responsiveUnits.sustInteractUnit = sustInteractUnit;

save(fullfile(newDir,'AVLoomingWhiteNoiseData.mat'),'unitData','analysisParams','responsiveUnits');
responsiveUnits


