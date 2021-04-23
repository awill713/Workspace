
clear;

experiment = 'EP010';
mouseID = 'AW165';
session = 'Session2';
date = '20210106-2';
stimPath = fullfile('D:\Electrophysiology\',experiment,mouseID,date,'StimInfo',[date '_' mouseID '_AVmultiContrastDriftingGratingsWhiteNoise_stimInfo']);

analysisWindow = [-100 1200]; %ms relative to stimulus onset
quantWindow = [0 300]; %ms relative to stimulus onset
baselineWindow = [-90 0]; %ms relative to stimulus onset

frBinWidth = 10; %ms

dataFolder = fullfile('D:\Electrophysiology\',experiment,mouseID,date,'SpikeMat');
dataFiles = dir(fullfile(dataFolder,'*AV_driftingGratingsMultiContrast_whiteNoise*'));


newDir = fullfile('D:\Electrophysiology\',experiment,mouseID,date,'AVMultiContrastDriftingGratingsWhiteNoise_final');
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
contrasts = stimInfo.contrasts;
bootstrapRepeats = 200;

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
orientationSelectiveUnitsLight = [];
directionSelectiveUnitsLight = [];
orientationSelectiveUnitsOR = [];
directionSelectiveUnitsOR = [];
orientationSelectiveUnitsAND = [];
directionSelectiveUnitsAND = [];
soundResponsiveUnits = [];
singleUnits = [];
multiUnits = [];

lightResponsiveUnitsGLM = [];
soundResponsiveUnitsGLM = [];
lightSoundInteractGLM = [];

rasterColorMap = 'hsv';
map = colormap(rasterColorMap);close;
rasterColors = map(round(linspace(1,length(map),length(orientations))),:);

for n = 1:totalUnits
    n
    nData = load(fullfile(dataFiles(n).folder,dataFiles(n).name));
    
    eventTimes = nData.Events;
    spikeTimes = nData.SpikeData(1,:);
    
    spikeRaster = cell(uniqueEvents,3); %spike times (1), trial number (y-axis) (2), color (3)
    meanResponse = zeros(uniqueEvents,5); %index (1),  mean baseline (2), baseline std (3), mean response (4), response std (5)
    trialResponse = zeros(uniqueEvents,repeats);
    frTrain = zeros(uniqueEvents,binCount);
    frTrainTrials = zeros(uniqueEvents,repeats,binCount);
    quantTrials = zeros(uniqueEvents,repeats);
    
    glmSound = zeros(uniqueEvents,repeats);
    glmLight = zeros(uniqueEvents,repeats);
    
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
            
            glmSound(u,e) = stimInfo.index(eventID,4);
            glmLight(u,e) = stimInfo.index(eventID,3);
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
    
    unitData(n).raster = spikeRaster;
    unitData(n).meanResponse = meanResponse;
    unitData(n).frTrain = frTrain;
    unitData(n).frTrainSTD = squeeze(std(frTrainTrials,[],2));
    unitData(n).frTrainTrials = frTrainTrials; %JUST NUMBER OF SPIKES IN BIN, NOT FIRING RATE
    unitData(n).trialResponse = trialResponse;
    unitData(n).sparsity = sparsity;
    unitData(n).type = nData.CellInfo(6);
    unitData(n).neuronNumber = nData.CellInfo(4);
    neuronNumber = nData.CellInfo(4);
    
    if nData.CellInfo(6)==1
        singleUnits = [singleUnits neuronNumber];
    elseif nData.CellInfo(6)==2
        multiUnits = [multiUnits neuronNumber];
    end
    
    totalTrials = length(stimInfo.order);
    lightCol = reshape(glmLight,[totalTrials 1]);
    soundCol = reshape(glmSound,[totalTrials 1]);
    responseCol = reshape(trialResponse/quantScalar,[totalTrials 1]);
    glmOutput = fitglm([lightCol soundCol],responseCol,'interactions','distr','poisson','link','log');
    glmPval = glmOutput.Coefficients.pValue(2:end);
    
    if glmPval(1)<0.05
        lightResponsiveUnitsGLM = [lightResponsiveUnitsGLM neuronNumber];
    end
    if glmPval(2)<0.05
        soundResponsiveUnitsGLM = [soundResponsiveUnitsGLM neuronNumber];
    end
    if glmPval(3)<0.05
        lightSoundInteractGLM = [lightSoundInteractGLM neuronNumber];
    end
    unitData(n).glmPval = glmPval;
    unitData(n).glmOutput = glmOutput;
    
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
    
    if (glmPval(1)<0.05) || (glmPval(3)<0.05)
        c=5;
        baseInd = (c-1)*length(orientations);
        visMeanResp = meanResponse(baseInd+1:baseInd+length(orientations),4);
        [~, maxIndV] = max(visMeanResp);
        oppIndV = mod(maxIndV+6-1,length(orientations))+1;
        orthIndV = [mod(maxIndV+3-1,length(orientations)) mod(maxIndV-3-1,length(orientations))]+1;
        
        avBaseInd = (c-1)*length(orientations)+length(orientations)*length(contrasts);
        audvisMeanResp = meanResponse(avBaseInd+1:avBaseInd+length(orientations),4);
        [~, maxIndAV] = max(audvisMeanResp);
        oppIndAV = mod(maxIndAV+6-1,length(orientations))+1;
        orthIndAV = [mod(maxIndAV+3-1,length(orientations)) mod(maxIndAV-3-1,length(orientations))]+1;
        
        maxResponsesV = trialResponse(baseInd+maxIndV,:);
        orthResponsesV = reshape(trialResponse(baseInd+orthIndV,:),[1 2*repeats]);
        oppResponsesV = trialResponse(baseInd+oppIndV,:);
        maxResponsesAV = trialResponse(avBaseInd+maxIndAV,:);
        orthResponsesAV = reshape(trialResponse(avBaseInd+orthIndAV,:),[1 2*repeats]);
        oppResponsesAV = trialResponse(avBaseInd+oppIndAV,:);
        
        
        tempDSIv = zeros(1,bootstrapRepeats);
        realDSIv = (meanResponse(baseInd+maxIndV,4)-meanResponse(baseInd+oppIndV,4))...
            / (meanResponse(baseInd+maxIndV,4)+meanResponse(baseInd+oppIndV,4));
        poolv = [trialResponse(baseInd+maxIndV,:) trialResponse(baseInd+oppIndV,:)];
        tempDSIav = zeros(1,bootstrapRepeats);
        realDSIav = (meanResponse(avBaseInd+maxIndAV,4)-meanResponse(avBaseInd+oppIndAV,4))...
            / (meanResponse(avBaseInd+maxIndAV,4)+meanResponse(avBaseInd+oppIndAV,4));
        poolav = [trialResponse(avBaseInd+maxIndAV,:) trialResponse(avBaseInd+oppIndAV,:)];
        for bsr = 1:bootstrapRepeats
            bootMaxIndices = randsample(length(poolv),repeats,'false');
            bootMax = poolv(bootMaxIndices);
            bootOpp = poolv(~ismember(1:length(poolv),bootMaxIndices));
            tDSI = abs((mean(bootMax)-mean(bootOpp))/(mean(bootMax)+mean(bootOpp)));
            tempDSIv(bsr) = tDSI;
            
            bootMaxIndices = randsample(length(poolav),repeats,'false');
            bootMax = poolav(bootMaxIndices);
            bootOpp = poolav(~ismember(1:length(poolav),bootMaxIndices));
            tDSI = abs((mean(bootMax)-mean(bootOpp))/(mean(bootMax)+mean(bootOpp)));
            tempDSIav(bsr) = tDSI;
        end
        if realDSIv > prctile(tempDSIv,95)
            directionSelectiveUnitsLight = [directionSelectiveUnitsLight neuronNumber];
            unitData(n).DSI = [realDSIv realDSIav];
        end
        if realDSIv > prctile(tempDSIv,95) || realDSIav > prctile(tempDSIav,95)
            directionSelectiveUnitsOR = [directionSelectiveUnitsOR neuronNumber];
            unitData(n).DSI = [realDSIv realDSIav];
        end
        if realDSIv > prctile(tempDSIv,95) && realDSIav > prctile(tempDSIav,95)
            directionSelectiveUnitsAND = [directionSelectiveUnitsAND neuronNumber];
            unitData(n).DSI = [realDSIv realDSIav];
        end
        
        tempOSIv = zeros(1,bootstrapRepeats);
        realOSIv = (mean(maxResponsesV)-mean(orthResponsesV))...
            / (mean(maxResponsesV)+mean(orthResponsesV));
        poolv = [maxResponsesV orthResponsesV];
        tempOSIav = zeros(1,bootstrapRepeats);
        realOSIav = (mean(maxResponsesAV)-mean(orthResponsesAV))...
            / (mean(maxResponsesAV)+mean(orthResponsesAV));
        poolav = [maxResponsesAV orthResponsesAV];
        for bsr = 1:bootstrapRepeats
            bootMaxIndices = randsample(length(poolv),repeats,'false');
            bootMax = poolv(bootMaxIndices);
            bootOrth = poolv(~ismember(1:length(poolv),bootMaxIndices));
            tOSI = abs((mean(bootMax)-mean(bootOrth))/(mean(bootMax)+mean(bootOrth)));
            tempOSIv(bsr) = tOSI;
            
            bootMaxIndices = randsample(length(poolav),repeats,'false');
            bootMax = poolav(bootMaxIndices);
            bootOrth = poolav(~ismember(1:length(poolav),bootMaxIndices));
            tOSI = abs((mean(bootMax)-mean(bootOrth))/(mean(bootMax)+mean(bootOrth)));
            tempOSIav(bsr) = tOSI;
        end
        if realOSIv > prctile(tempOSIv,95)
            orientationSelectiveUnitsLight = [orientationSelectiveUnitsLight neuronNumber];
            unitData(n).OSI = [realOSIv realOSIav];
        end
        if realOSIv > prctile(tempOSIv,95) || realOSIav > prctile(tempOSIav,95)
            orientationSelectiveUnitsOR = [orientationSelectiveUnitsOR neuronNumber];
            unitData(n).OSI = [realOSIv realOSIav];
        end
        if realOSIv > prctile(tempOSIv,95) && realOSIav > prctile(tempOSIav,95)
            orientationSelectiveUnitsAND = [orientationSelectiveUnitsAND neuronNumber];
            unitData(n).OSI = [realOSIv realOSIav];
        end
    else
        pDir = 1;
        pOrient=1;
        realOSIv = 0;
        realDSIv = 0;
    end
        
    
    f1 = figure;
    set(f1,'Position',[150 80 1250 700]);
    for c = 1:length(contrasts)
        subplot(3,length(contrasts),c);
        for r = 1:length(orientations)
            evIdx = r + (c-1)*length(orientations);
            scatter(spikeRaster{evIdx,1}/1000,spikeRaster{evIdx,2},[],spikeRaster{evIdx,3},'.');
            hold on;
        end
        title(['Cont = ' num2str(contrasts(c))]);
        
        subplot(3,length(contrasts),c+length(contrasts));
        for r = 1:length(orientations)
            evIdx = r + (c-1)*length(orientations) + length(orientations)*length(contrasts);
            scatter(spikeRaster{evIdx,1},spikeRaster{evIdx,2},[],spikeRaster{evIdx,3},'.');
            hold on;
        end
        
        subplot(3,length(contrasts),c+2*length(contrasts));hold on;
        psthScalar = 20/max(max(frTrain));
        yMin = 0;
        for r = 1:length(orientations)
            evIdx = r + (c-1)*length(orientations);
            idxOffset = length(orientations)*length(contrasts);
            
            trace = (frTrain(evIdx,:)-mean(frTrain(evIdx,baselineFirstBin:baselineLastBin)))*psthScalar + 20*(r-1);
            plot(binEdges,trace,'Color',spikeRaster{evIdx,3});
            yMin = min([yMin min(trace)]);
            
            trace = (frTrain(evIdx+idxOffset,:)-mean(frTrain(evIdx+idxOffset,baselineFirstBin:baselineLastBin)))*psthScalar + 20*(r-1);
            plot(binEdges,trace,'Color',spikeRaster{evIdx+idxOffset,3},'LineStyle','--');
            yMin = min([yMin min(trace)]);
        end
        
        %         yLimits = ylim();
        %         ylim([yMin-5  yLimits(2)]);
        yticks([]);
        xlim([analysisWindow]);
        xlabel('Time (ms)');
    end
    suptitle({['Unit ' num2str(nData.CellInfo(4)) ' (n = ' num2str(n)...
        '), unitType = ' num2str(nData.CellInfo(6))] ['Light responsive = ' num2str(pAV(1))...
        ', Sound responsive = ' num2str(pAV(2)) ', Orientation selective = ' num2str(realOSIv)...
        ', Direction selective = ' num2str(realDSIv)] ['Light GLM = ' num2str(glmPval(1))...
        ', Sound GLM = ' num2str(glmPval(2))]});
    
%     saveas(f1,fullfile(figureDir,['Unit ' num2str(nData.CellInfo(4)) ' response rasters and PSTH.fig']));
%     saveas(f1,fullfile(figureDir,['Unit ' num2str(nData.CellInfo(4)) ' response rasters and PSTH.jpg']));
    close(f1);
    
    
    f2 = figure;
    set(f2,'Position',[25 200 1500 350]);
    limMax = 0;
    for c = 1:length(contrasts)
        subplot(1,length(contrasts),c);
        
        evIdx = 1 + (c-1)*length(orientations);
        idxOffset = length(orientations)*length(contrasts);
        
        orientCoords =[2*pi*orientations./360 2*pi*orientations(1)/360];
        vPolar = [meanResponse(evIdx:(evIdx+length(orientations)-1),4); meanResponse(evIdx,4)];
        avPolar = [meanResponse((evIdx+idxOffset):(evIdx+idxOffset+length(orientations)-1),4); meanResponse(evIdx+idxOffset,4)];
        polarplot(orientCoords,vPolar,'Color',[0 0 0],'LineWidth',2);
        hold on;
        polarplot(orientCoords,avPolar,'Color',[0 0 1],'LineWidth',2);
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
        ', Sound responsive = ' num2str(pAV(2)) ', Orientation selective = ' num2str(realOSIv)...
        ', Direction selective = ' num2str(realDSIv)] ['Light GLM = ' num2str(glmPval(1))...
        ', Sound GLM = ' num2str(glmPval(2))]});
%     saveas(f2,fullfile(figureDir,['Unit ' num2str(nData.CellInfo(4)) ' response polar plots.fig']));
%     saveas(f2,fullfile(figureDir,['Unit ' num2str(nData.CellInfo(4)) ' response polar plots.jpg']));
    close(f2);
    
    
    f5 = figure; hold on;
    set(f5,'Position',[400 150 800 600]);
    fanoFactor = [];
    for c = 1:length(contrasts)
        %         subplot(1,length(contrasts),c); hold on;
        
        evIdx = 1 + (c-1)*length(orientations);
        idxOffset = length(orientations)*length(contrasts);
        
        vis = meanResponse(evIdx:evIdx+length(orientations)-1,4)';
        visAud = meanResponse(evIdx+idxOffset:evIdx+idxOffset+length(orientations)-1,4)';
        stdVis = ((meanResponse(evIdx:evIdx+length(orientations)-1,5)*sqrt(repeats))');
        stdAudVis = ((meanResponse(evIdx+idxOffset:evIdx+idxOffset+length(orientations)-1,5)*sqrt(repeats))');
        
        sAlpha = scatter(vis,stdVis,[],[0 0 0],'filled');
        sAlpha.MarkerFaceAlpha = c/length(contrasts);
        sAlpha = scatter(visAud,stdAudVis,[],[0 0 1],'filled');
        sAlpha.MarkerFaceAlpha = c/length(contrasts);
        xlabel('Firing rate (Hz)');
        ylabel('Standard deviation');
        
        fanoVis = mean(stdVis)/mean(vis);
        fanoAudVis = mean(stdAudVis)/mean(visAud);
        fanoFactor = [fanoFactor [fanoVis; fanoAudVis]];
        
    end
    vis = meanResponse(1:length(orientations)*length(contrasts),4)';
    visAud = meanResponse(length(orientations)*length(contrasts)+1:length(orientations)*length(contrasts)*2,4)';
    stdVis = (meanResponse(1:length(orientations)*length(contrasts),5)'*sqrt(repeats));
    stdAudVis = (meanResponse(length(orientations)*length(contrasts)+1:length(orientations)*length(contrasts)*2,5)'*sqrt(repeats));
    pVis = polyfit(vis,stdVis,1);
    pVisAud = polyfit(visAud,stdAudVis,1);
    lineMeanV = [mean(vis) mean(stdVis)];
    lineMeanAV = [mean(visAud) mean(stdAudVis)];
    xax = xlim;
    yax = ylim;
    line([0 xax(2)],[0 xax(2)*lineMeanV(2)/lineMeanV(1)],'Color',[0 0 0]);
    line([0 xax(2)],[0 xax(2)*lineMeanAV(2)/lineMeanAV(1)],'Color',[0 0 1]);
    
    plot([min(vis) max(vis)],polyval(pVis,[min(vis) max(vis)]),'Color',[0 0 0]);
    plot([min(visAud) max(visAud)],polyval(pVisAud,[min(visAud) max(visAud)]),'Color',[0 0 1]);
    unitData(n).fanoFactor = fanoFactor;
    title(['m_V, b_V = ' num2str(pVis) ',    m_A_V, b_A_V = ' num2str(pVisAud)]);
    suptitle({['Unit ' num2str(nData.CellInfo(4)) ' (n = ' num2str(n)...
        '), unitType = ' num2str(nData.CellInfo(6))] ['Light responsive = ' num2str(pAV(1))...
        ', Sound responsive = ' num2str(pAV(2)) ', Orientation selective = ' num2str(realOSIv)...
        ', Direction selective = ' num2str(realDSIv)] ['Light GLM = ' num2str(glmPval(1))...
        ', Sound GLM = ' num2str(glmPval(2))]});
%     saveas(f5,fullfile(figureDir,['Unit ' num2str(nData.CellInfo(4)) ' response variance.fig']));
%     saveas(f5,fullfile(figureDir,['Unit ' num2str(nData.CellInfo(4)) ' response variance.jpg']));
    close(f5);
    
end


analysisParams.mouseID = mouseID;
analysisParams.session = session;
analysisParams.date = date;
analysisParams.stimPath = stimPath;
analysisParams.analysisWindow = analysisWindow;
analysisParams.frBinWidth = frBinWidth;
analysisParams.quantWindow = quantWindow;
analysisParams.baselineWindow = baselineWindow;
analysisParams.quantScalar = quantScalar;
analysisParams.binEdges = binEdges;

responsiveUnits.soundResponsiveUnits = intersect(soundResponsiveUnits,[singleUnits multiUnits]);
responsiveUnits.lightResponsiveUnits = intersect(lightResponsiveUnits,[singleUnits multiUnits]);
responsiveUnits.directionSelectiveUnits = intersect(directionSelectiveUnitsLight,[singleUnits multiUnits]);
responsiveUnits.directionSelectiveUnitsOR = intersect(directionSelectiveUnitsOR,[singleUnits multiUnits]);
responsiveUnits.directionSelectiveUnitsAND = intersect(directionSelectiveUnitsAND,[singleUnits multiUnits]);
responsiveUnits.orientationSelectiveUnits = intersect(orientationSelectiveUnitsLight,[singleUnits multiUnits]);
responsiveUnits.orientationSelectiveUnitsOR = intersect(orientationSelectiveUnitsOR,[singleUnits multiUnits]);
responsiveUnits.orientationSelectiveUnitsAND = intersect(orientationSelectiveUnitsAND,[singleUnits multiUnits]);
responsiveUnits.multiUnits = multiUnits;
responsiveUnits.singleUnits = singleUnits;

responsiveUnits.soundAndLight = intersect(responsiveUnits.soundResponsiveUnits,responsiveUnits.lightResponsiveUnits);
responsiveUnits.soundAndOrientation = intersect(responsiveUnits.soundResponsiveUnits,responsiveUnits.orientationSelectiveUnits);
responsiveUnits.lightAndOrientation = intersect(responsiveUnits.lightResponsiveUnits,responsiveUnits.orientationSelectiveUnits);

responsiveUnitsGLM.soundResponsiveUnits = intersect(soundResponsiveUnitsGLM,[singleUnits multiUnits]);
responsiveUnitsGLM.lightResponsiveUnits = intersect(lightResponsiveUnitsGLM,[singleUnits multiUnits]);
responsiveUnitsGLM.lightSoundInteractUnits = intersect(lightSoundInteractGLM,[singleUnits multiUnits]);

normalizedResponseCurves = cell(length(contrasts),2);
orientationCurves = cell(length(contrasts),2);
for n = 1:totalUnits
    if unitData(n).type~=0 && ismember(unitData(n).neuronNumber,responsiveUnits.lightResponsiveUnits)
        
        allResponses = unitData(n).meanResponse(:,4);
        highContVis = allResponses(((length(contrasts)-1)*length(orientations))+1:((length(contrasts)-1)*length(orientations)+length(orientations)));
        [peak ind] = max(highContVis);
        maxVal = max(allResponses(:));
        minVal = min(allResponses(:));
        
        for c = 1:length(contrasts)
            %             responseVis = unitData(n).meanResponse(((c-1)*length(orientations))+1:((c-1)*length(orientations)+length(orientations)),4);
            %             responseVisAud = unitData(n).meanResponse(((c-1)*length(orientations)+1+length(contrasts)*length(orientations)):((c-1)*length(orientations)+length(contrasts)*length(orientations)+length(orientations)),4);
            
            responseVis = allResponses(((c-1)*length(orientations))+1:((c-1)*length(orientations)+length(orientations)));
            responseVisAud = allResponses(((c-1)*length(orientations)+1+length(contrasts)*length(orientations)):((c-1)*length(orientations)+length(contrasts)*length(orientations)+length(orientations)));
            
            %             [peak ind] = max(responseVis);
            targetIndex = floor(length(orientations)/2) + 1;
            
            realignedVis = circshift(responseVis,targetIndex - ind);
            realignedVis = [realignedVis; realignedVis(1)];
            realignedVisAud = circshift(responseVisAud,targetIndex - ind);
            realignedVisAud = [realignedVisAud; realignedVisAud(1)];
            
            %             minVal = min(realignedVis); maxVal = max(realignedVis);
            %             if maxVal==minVal
            %                 maxVal = maxVal + 1;
            %             end
            normVis = (realignedVis - minVal) / (maxVal - minVal);
            normVisAud = (realignedVisAud - minVal) / (maxVal - minVal);
            
            normalizedResponseCurves{c,1} = [normalizedResponseCurves{c,1}; normVis'];
            normalizedResponseCurves{c,2} = [normalizedResponseCurves{c,2}; normVisAud'];
            
            vNorm = (responseVis - minVal) / (maxVal - minVal);
            vaNorm = (responseVisAud - minVal) / (maxVal - minVal);
            
            orientationCurves{c,1} = [orientationCurves{c,1}; vNorm'];
            orientationCurves{c,2} = [orientationCurves{c,2}; vaNorm'];
        end
    end
end

f3 = figure;
set(f3,'Position',[25 200 1500 400]);
for c = 1:length(contrasts)
    subplot(1,length(contrasts),c);hold on;
    plot(linspace(-180,180,length(orientations)+1),mean(normalizedResponseCurves{c,1}),'Color',[0 0 0],'LineWidth',2);
    plot(linspace(-180,180,length(orientations)+1),mean(normalizedResponseCurves{c,2}),'Color',[0 0 1],'LineWidth',2);
    ylabel('Normalized firing rate');
    xlabel('\Delta orientation (degrees)');
    xticks(-180:90:180);
    ylim([0 1])
    title(['Cont = ' num2str(contrasts(c))]);
end
legend({'Light','Light/Sound'});
suptitle({'Neuron-wise orientation preference, +/- white noise'});
% saveas(f3,fullfile(newDir,'Neuron-wise orientation preference'));

f4 = figure;
set(f4,'Position',[25 200 1500 400]);
for c = 1:length(contrasts)
    subplot(1,length(contrasts),c);
    polarplot(2*pi*orientations./360,mean(orientationCurves{c,1}),'Color',[0 0 0],'LineWidth',2);
    hold on;
    polarplot(2*pi*orientations./360,mean(orientationCurves{c,2}),'Color',[0 0 1],'LineWidth',2);
    rlim([0 1]);
    title(['Cont = ' num2str(contrasts(c))]);
end
suptitle({'Population orientation preference, +/- white noise'});
% saveas(f4,fullfile(newDir,'Population orientation preference'));
%
% f6 = figure;hold on;
% slopes = zeros(0,2);
% yint = zeros(0,2);
% for n = 1:totalUnits
%     if ismember(unitData(n).neuronNumber,soundResponsiveUnits) && unitData(n).type~=0
% %     if unitData(n).type~=0
%         slopes = [slopes; unitData(n).fanoFit(1) unitData(n).fanoFit(3)];
%         yint = [yint; unitData(n).fanoFit(2) unitData(n).fanoFit(4)];
%     end
% end
% scatter(slopes(:,1),yint(:,1),[],[0 0 0],'filled');
% scatter(slopes(:,2),yint(:,2),[],[0 0 1],'filled');
% [hslopes pslopes] = ttest(slopes(:,1),slopes(:,2));
% [hint pint] = ttest(yint(:,1),yint(:,2));
% xlabel('Slope');
% ylabel('Y-intercept');
% title(['Slopes p = ' num2str(pslopes) ', y-int p = ' num2str(pint)]);
% suptitle({'Mean response to variance relationship'});
% saveas(f6,fullfile(newDir,'Population variance analysis'));
% subplot(1,2,1);hold on;
% histogram(slopes(:,1));
% histogram(slopes(:,2));
% [~, pSlopes] = ttest2(slopes(:,1),slopes(:,2));
% title(['p = ' num2str(pSlopes)]);
%
% subplot(1,2,2);hold on;
% histogram(yint(:,1));
% histogram(yint(:,2));
% [~, pYint] = ttest2(yint(:,1),yint(:,2));
% title(['p = ' num2str(pYint)]);



save(fullfile(newDir,'AVMultiContrastDriftingGratingsWhiteNoiseData_selectivity.mat'),'unitData','analysisParams','responsiveUnits','responsiveUnitsGLM');
responsiveUnitsGLM



%%
% theOSIs = [];
% theDSIs = [];
% for n = 1:length(unitData)
%     neuronNumber = unitData(n).neuronNumber;
%     if ismember(neuronNumber,responsiveUnits.orientationSelectiveUnits)
%         v = unitData(n).OSI(1);
%         av = unitData(n).OSI(2);
%         theOSIs = [theOSIs; v av];
%     end
%     if ismember(neuronNumber,responsiveUnits.directionSelectiveUnits)
%         v = unitData(n).DSI(1);
%         av = unitData(n).DSI(2);
%         theDSIs = [theDSIs; v av];
%     end
% end
% figure;
% scatter(theOSIs(:,1),theOSIs(:,2));
% title('OSI');
% figure;
% scatter(theDSIs(:,1),theDSIs(:,2));
% title('DSI');
