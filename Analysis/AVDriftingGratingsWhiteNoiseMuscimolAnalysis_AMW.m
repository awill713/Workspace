
clear;

experiment = 'EP009';
mouseID = 'AW149';
session = 'Session1';
date = '20200923';
preStimPath = fullfile('D:\Electrophysiology\',experiment,mouseID,date,'StimInfo',[date '_' mouseID '_AVsingleContrastDriftingGratingsWhiteNoise_stimInfo_pre']);
postStimPath = fullfile('D:\Electrophysiology\',experiment,mouseID,date,'StimInfo',[date '_' mouseID '_AVsingleContrastDriftingGratingsWhiteNoise_stimInfo_post']);

analysisWindow = [-100 1200]; %ms relative to stimulus onset
quantWindow = [0 300]; %ms relative to stimulus onset
baselineWindow = [-90 0]; %ms relative to stimulus onset

frBinWidth = 10; %ms

dataFolder = fullfile('D:\Electrophysiology\',experiment,mouseID,date,'SpikeMat');
dataFilesPre = dir(fullfile(dataFolder,'*avDriftingGratingsSingleContrast_whiteNoise_pre*'));
dataFilesPost = dir(fullfile(dataFolder,'*avDriftingGratingsSingleContrast_whiteNoise_post*'));

% newDir = fullfile('D:\KiloSort\',mouseID,session,folder,'OptoNoiseResponses');
newDir = fullfile('D:\Electrophysiology\',experiment,mouseID,date,'AVMultiContrastDriftingGratingsWhiteNoise_Recat');
figureDir = fullfile(newDir,'Figures');
if ~exist(newDir)
    mkdir(newDir);
end
if ~exist(figureDir)
    mkdir(figureDir);
end

load(preStimPath);
uniqueEvents = size(stimInfo.index,1);
indices = stimInfo.index(:,1);
repeats = stimInfo.repeats;
orientations = stimInfo.orientations;
contrasts = stimInfo.contrasts;

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

lightResponsiveUnits = [];
orientationSelectiveUnits = [];
directionSelectiveUnits = [];
soundResponsiveUnits = [];
singleUnits = [];
multiUnits = [];

rasterColorMap = 'hsv';
map = colormap(rasterColorMap);close;
rasterColors = map(round(linspace(1,length(map),length(orientations))),:);

for n = 1:totalUnits
    n
    nData = load(fullfile(dataFilesPre(n).folder,dataFilesPre(n).name));
    
    eventTimes = nData.Events;
    spikeTimes = nData.SpikeData(1,:);
    
    spikeRaster = cell(uniqueEvents,3); %spike times (1), trial number (y-axis) (2), color (3)
    meanResponse = zeros(uniqueEvents,5); %index (1),  mean baseline (2), baseline std (3), mean response (4), response std (5)
    trialResponse = zeros(uniqueEvents,repeats);
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
        trialResponse(u,:) = prePost(:,2)*quantScalar;
        
        frTrainTrials(u,:,:) = spikeFR; %JUST NUMBER OF SPIKES IN BIN, NOT FIRING RATE
        frTrain(u,:) = mean(spikeFR,1)*frScalar;
        
    end
    
    sparsity = [];
    for c = 1:length(contrasts)
        evIdx = 1 + (c-1)*length(orientations);
        idxOffset = length(orientations)*length(contrasts);
        
        v = meanResponse(evIdx:(evIdx+length(orientations)-1),4);
        av = meanResponse(evIdx+idxOffset:(evIdx+length(orientations)-1+idxOffset),4);
        
        vSparse = 1 - (sum(v./length(orientations))^2 / sum(v.^2 / length(orientations)));
        avSparse = 1 - (sum(av./length(orientations))^2 / sum(av.^2 / length(orientations)));
        
        sparsity = [sparsity [vSparse; avSparse]];
    end
    
    unitData(n).PREraster = spikeRaster;
    unitData(n).PREmeanResponse = meanResponse;
    unitData(n).PREfrTrain = frTrain;
    unitData(n).PREfrTrainSTD = squeeze(std(frTrainTrials,[],2));
    unitData(n).PREfrTrainTrials = frTrainTrials; %JUST NUMBER OF SPIKES IN BIN, NOT FIRING RATE
    unitData(n).PREtrialResponse = trialResponse;
    unitData(n).PREsparsity = sparsity;
    unitData(n).type = nData.CellInfo(6);
    unitData(n).neuronNumber = nData.CellInfo(4);
    neuronNumber = nData.CellInfo(4);
    
    if nData.CellInfo(6)==1
        singleUnits = [singleUnits neuronNumber];
    elseif nData.CellInfo(6)==2
        multiUnits = [multiUnits neuronNumber];
    end
    
    anovaMat = [];
    for c = 1:length(contrasts)
        baseInd = (c-1)*length(orientations);
        avBaseInd = (c-1)*length(orientations)+length(orientations)*length(contrasts);
        vColumn = reshape(trialResponse(baseInd+1:baseInd+length(orientations),:), [repeats*length(orientations) 1]);
        avColumn = reshape(trialResponse(avBaseInd+1:avBaseInd+length(orientations),:), [repeats*length(orientations), 1]);
        
        anovaMat = [anovaMat [vColumn; avColumn]];
    end
    pAV = anova2(anovaMat,repeats*length(orientations),'off');
    if pAV(1)<0.05
        lightResponsiveUnits = [lightResponsiveUnits neuronNumber];
    end
    if pAV(2)<0.05 || pAV(3)<0.05
        soundResponsiveUnits = [soundResponsiveUnits neuronNumber];
    end
    
    if pAV(1)<0.05
        c=2;
        baseInd = (c-1)*length(orientations);
        visMeanResp = meanResponse(baseInd+1:baseInd+length(orientations),4);
        [~, maxInd] = max(visMeanResp);
        oppInd = mod(maxInd+6,length(orientations));
        orthInd = [mod(maxInd+3,length(orientations)) mod(maxInd-3,length(orientations))];
        
        
        [~, pDir] = ttest2(trialResponse(baseInd+maxInd,:),trialResponse(baseInd+oppInd,:));
        if pDir<0.05
            directionSelectiveUnits = [directionSelectiveUnits neuronNumber];
        end
        
        [~, pOrient] = ttest2(trialResponse(baseInd+maxInd,:),[trialResponse(baseInd+orthInd(1),:) trialResponse(baseInd+orthInd(2),:)]);
        if pOrient<0.05
            orientationSelectiveUnits = [orientationSelectiveUnits neuronNumber];
        end
    else
        pDir = 1;
        pOrient=1;
    end
    
    %% post inj analysis
    postStim = load(postStimPath);
    nDataPost = load(fullfile(dataFilesPre(n).folder,dataFilesPost(n).name));
    
    eventTimes = nDataPost.Events;
    spikeTimes = nDataPost.SpikeData(1,:);
    
    spikeRasterPost = cell(uniqueEvents,3); %spike times (1), trial number (y-axis) (2), color (3)
    meanResponsePost = zeros(uniqueEvents,5); %index (1),  mean baseline (2), baseline std (3), mean response (4), response std (5)
    trialResponsePost = zeros(uniqueEvents,repeats);
    frTrainPost = zeros(uniqueEvents,binCount);
    frTrainTrialsPost = zeros(uniqueEvents,repeats,binCount);
    quantTrialsPost = zeros(uniqueEvents,repeats);
    
    for u = 1:uniqueEvents
        eventID = indices(u);
        eventsOfInterest = find(postStim.stimInfo.order==eventID);
        
        
        color = rasterColors(mod(u-1,size(rasterColors,1))+1,:);
        
        
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
        
        spikeRasterPost{u,1} = rasterX;
        spikeRasterPost{u,2} = rasterY;
        spikeRasterPost{u,3} = color;
        meanResponsePost(u,1) = eventID;
        meanResponsePost(u,2) = mean(prePost(:,1))*baselineScalar;
        meanResponsePost(u,3) = std(prePost(:,1))*baselineScalar/sqrt(repeats);
        meanResponsePost(u,4) = mean(prePost(:,2))*quantScalar;
        meanResponsePost(u,5) = std(prePost(:,2))*quantScalar/sqrt(repeats);
        trialResponsePost(u,:) = prePost(:,2)*quantScalar;
        
        frTrainTrialsPost(u,:,:) = spikeFR; %JUST NUMBER OF SPIKES IN BIN, NOT FIRING RATE
        frTrainPost(u,:) = mean(spikeFR,1)*frScalar;
        
    end
    
    sparsityPost = [];
    for c = 1:length(contrasts)
        evIdx = 1 + (c-1)*length(orientations);
        idxOffset = length(orientations)*length(contrasts);
        
        v = meanResponse(evIdx:(evIdx+length(orientations)-1),4);
        av = meanResponse(evIdx+idxOffset:(evIdx+length(orientations)-1+idxOffset),4);
        
        vSparse = 1 - (sum(v./length(orientations))^2 / sum(v.^2 / length(orientations)));
        avSparse = 1 - (sum(av./length(orientations))^2 / sum(av.^2 / length(orientations)));
        
        sparsityPost = [sparsityPost [vSparse; avSparse]];
    end
    
    unitData(n).POSTraster = spikeRasterPost;
    unitData(n).POSTmeanResponse = meanResponsePost;
    unitData(n).POSTfrTrain = frTrainPost;
    unitData(n).POSTfrTrainSTD = squeeze(std(frTrainTrialsPost,[],2));
    unitData(n).POSTfrTrainTrials = frTrainTrialsPost; %JUST NUMBER OF SPIKES IN BIN, NOT FIRING RATE
    unitData(n).POSTtrialResponse = trialResponsePost;
    unitData(n).POSTsparsity = sparsityPost;
    
    %% plot
        
    
    f1 = figure;
    set(f1,'Position',[150 80 1250 700]);
    for c = 1:length(contrasts)
        subplot(4,length(contrasts),c);
        for r = 1:length(orientations)
            evIdx = r + (c-1)*length(orientations);
            scatter(spikeRaster{evIdx,1}/1000,spikeRaster{evIdx,2},[],spikeRaster{evIdx,3},'.');
            hold on;
        end
        title(['Cont = ' num2str(contrasts(c)) ', light']);
        
        subplot(4,length(contrasts),c+length(contrasts));
        for r = 1:length(orientations)
            evIdx = r + (c-1)*length(orientations) + length(orientations)*length(contrasts);
            scatter(spikeRaster{evIdx,1},spikeRaster{evIdx,2},[],spikeRaster{evIdx,3},'.');
            hold on;
        end
        title(['Cont = ' num2str(contrasts(c)) ', light+sound']);
        
        subplot(4,length(contrasts),c+2*length(contrasts));
        for r = 1:length(orientations)
            evIdx = r + (c-1)*length(orientations);
            scatter(spikeRasterPost{evIdx,1}/1000,spikeRasterPost{evIdx,2},[],spikeRasterPost{evIdx,3},'.');
            hold on;
        end
        title(['Cont = ' num2str(contrasts(c)) ', light+musc']);
        
        subplot(4,length(contrasts),c+3*length(contrasts));
        for r = 1:length(orientations)
            evIdx = r + (c-1)*length(orientations) + length(orientations)*length(contrasts);
            scatter(spikeRasterPost{evIdx,1},spikeRasterPost{evIdx,2},[],spikeRasterPost{evIdx,3},'.');
            hold on;
        end
        title(['Cont = ' num2str(contrasts(c)) ', light+sound+musc']);
    end
    suptitle({['Unit ' num2str(nData.CellInfo(4)) ' (n = ' num2str(n)...
        '), unitType = ' num2str(nData.CellInfo(6))] ['Light responsive = ' num2str(pAV(1))...
        ', Sound responsive = ' num2str(pAV(2)) ', Orientation selective = ' num2str(pOrient)...
        ', Direction selective = ' num2str(pDir)]});
    
    saveas(f1,fullfile(figureDir,['Unit ' num2str(nData.CellInfo(4)) ' response rasters.fig']));
    saveas(f1,fullfile(figureDir,['Unit ' num2str(nData.CellInfo(4)) ' response rasters.jpg']));
    close(f1);
    
    
    f2 = figure;
    set(f2,'Position',[25 200 1500 350]);
    limMax = 0;
    for c = 1:length(contrasts)
        subplot(1,length(contrasts),c);
        
        evIdx = 1 + (c-1)*length(orientations);
        idxOffset = length(orientations)*length(contrasts);
        
        polarplot(2*pi*orientations./360,meanResponse(evIdx:(evIdx+length(orientations)-1),4),'Color',[0 0 0],'LineWidth',2);
        hold on;
        polarplot(2*pi*orientations./360,meanResponse((evIdx+idxOffset):(evIdx+idxOffset+length(orientations)-1),4),'Color',[0 0 1],'LineWidth',2);
        polarplot(2*pi*orientations./360,meanResponsePost((evIdx+idxOffset):(evIdx+idxOffset+length(orientations)-1),4),'Color',[0 0 1],'LineWidth',2,'LineStyle',':');
        title(['Cont = ' num2str(contrasts(c))]);
        tempLim = rlim;
        if max(tempLim)>limMax
            limMax = max(tempLim);
        end
    end
    for c = 1:length(contrasts)
        subplot(1,length(contrasts),c);
        rlim([0 limMax]);
    end
    suptitle({['Unit ' num2str(nData.CellInfo(4)) ' (n = ' num2str(n)...
        '), unitType = ' num2str(nData.CellInfo(6))] ['Light responsive = ' num2str(pAV(1))...
        ', Sound responsive = ' num2str(pAV(2)) ', Orientation selective = ' num2str(pOrient)...
        ', Direction selective = ' num2str(pDir)]});
    saveas(f2,fullfile(figureDir,['Unit ' num2str(nData.CellInfo(4)) ' response polar plots.fig']));
    saveas(f2,fullfile(figureDir,['Unit ' num2str(nData.CellInfo(4)) ' response polar plots.jpg']));
    close(f2);
    
   
end


analysisParams.mouseID = mouseID;
analysisParams.session = session;
analysisParams.date = date;
analysisParams.preStimPath = preStimPath;
analysisParams.postStimPath = postStimPath;
analysisParams.analysisWindow = analysisWindow;
analysisParams.frBinWidth = frBinWidth;
analysisParams.quantWindow = quantWindow;
analysisParams.baselineWindow = baselineWindow;

responsiveUnits.soundResponsiveUnits = intersect(soundResponsiveUnits,[singleUnits multiUnits]);
responsiveUnits.lightResponsiveUnits = intersect(lightResponsiveUnits,[singleUnits multiUnits]);
responsiveUnits.directionSelectiveUnits = intersect(directionSelectiveUnits,[singleUnits multiUnits]);
responsiveUnits.orientationSelectiveUnits = intersect(orientationSelectiveUnits,[singleUnits multiUnits]);
responsiveUnits.multiUnits = multiUnits;
responsiveUnits.singleUnits = singleUnits;

responsiveUnits.soundAndLight = intersect(responsiveUnits.soundResponsiveUnits,responsiveUnits.lightResponsiveUnits);
responsiveUnits.soundAndOrientation = intersect(responsiveUnits.soundResponsiveUnits,responsiveUnits.orientationSelectiveUnits);
responsiveUnits.lightAndOrientation = intersect(responsiveUnits.lightResponsiveUnits,responsiveUnits.orientationSelectiveUnits);



preL = [];
preB = [];
postB = [];
postL = [];
lblb = [];
for n = 1:totalUnits
    if unitData(n).type~=0 && ismember(unitData(n).neuronNumber,responsiveUnits.soundAndLight)
        
        preLight = unitData(n).PREmeanResponse(25:36,4);
        preBoth = unitData(n).PREmeanResponse(37:48,4);
        postBoth = unitData(n).POSTmeanResponse(37:48,4);
        postLight = unitData(n).POSTmeanResponse(25:36,4);
        
        preL = [preL; preLight'];
        preB = [preB; preBoth'];
        postB = [postB; postBoth'];
        postL = [postL; postLight'];
        
        lblb = [lblb; mean(preLight) mean(preBoth) mean(postLight) mean(postBoth)];
        
    end
end
mLPre = mean(preL);
mBPre = mean(preB);
mLPost = mean(postL);
mBPost = mean(postB);
     
orientationsL = [orientations orientations(1)];
f4 = figure;
    polarplot(2*pi*orientationsL./360,[mLPre mLPre(1)],'Color',[0 0 0],'LineWidth',2);
    hold on;
    polarplot(2*pi*orientationsL./360,[mBPre mBPre(1)],'Color',[0 0 1],'LineWidth',2);
    polarplot(2*pi*orientationsL./360,[mLPost mLPost(1)],'Color',[0 0 0],'LineWidth',2,'LineStyle',':');
    
    polarplot(2*pi*orientationsL./360,[mBPost mBPost(1)],'Color',[0 0 1],'LineWidth',2,'LineStyle',':');
    title(['Cont = ' num2str(contrasts(c))]);
suptitle({'Population orientation preference, +/- white noise, muscimol'});
saveas(f4,fullfile(newDir,'Population orientation preference'));

f6 = figure;
scatter(lblb(:,2)./lblb(:,1),lblb(:,4)./lblb(:,3));
hold on;
line([0 4],[0 4]);
xlabel('Audiovisual response ratio, pre-muscimol');
ylabel('Audiovisual response ratio, post-muscimol');
saveas(f6,fullfile(newDir,'Pre and post AV response ratio'));



save(fullfile(newDir,'AVMultiContrastDriftingGratingsWhiteNoiseData.mat'),'unitData','analysisParams','responsiveUnits');
responsiveUnits

