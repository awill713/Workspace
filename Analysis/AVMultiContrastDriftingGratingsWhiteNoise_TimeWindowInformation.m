
clear all;
% close all;

saveDir = fullfile('E:\Electrophysiology','EP004','MetaData - AVMultiContrastDriftingGratingsWhiteNoise');
if ~exist(saveDir)
    mkdir(saveDir);
end

dataPaths{1} = fullfile('EP004','AW117','20200221-1');
dataPaths{2} = fullfile('EP004','AW117','20200221-2');
dataPaths{3} = fullfile('EP004','AW118','20200221-1');
dataPaths{4} = fullfile('EP004','AW118','20200221-2');
dataPaths{5} = fullfile('EP004','AW121','20200226-1');
dataPaths{6} = fullfile('EP004','AW121','20200226-2');
dataPaths{7} = fullfile('EP004','AW124','20200303-1');
dataPaths{8} = fullfile('EP004','AW124','20200303-2');

%which neurons to include
onlySingleUnits = 0;
allNeurons = 0;
lightResponsive = 1;
orientationSelective = 1;
soundResponsive = 1;


for dp = 1:length(dataPaths)
    dp
    %Load data
    dataFile = fullfile('E:\Electrophysiology\',dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise','AVMultiContrastDriftingGratingsWhiteNoiseData.mat');
    load(dataFile);
    stimPath = dir(fullfile('E:\Electrophysiology\',dataPaths{dp},'StimInfo','*AVmultiContrastDriftingGratingsWhiteNoise_stimInfo*'));
    load(fullfile(stimPath.folder,stimPath.name));
    
    orientations = stimInfo.orientations;
    contrasts = stimInfo.contrasts;
    
    binCount = analysisParams.analysisWindow(2)-analysisParams.analysisWindow(1)-analysisParams.frBinWidth+1; %sliding window
    binEdges = analysisParams.analysisWindow(1)+analysisParams.frBinWidth:1:analysisParams.analysisWindow(2);
    
    if ~exist('orientationCurves','var')
        orientationInformation = cell(length(contrasts),2);
        lightInformation = cell(1,2);
    end
    
    for u  = 1:size(unitData,2)
        neuronNumber = unitData(u).neuronNumber;
        %% Light information
        if (soundResponsive == ismember(neuronNumber,responsiveUnits.soundResponsiveUnits))...
                && (lightResponsive == ismember(neuronNumber,responsiveUnits.lightResponsiveUnits))...
                && unitData(u).type~=0
            
            responses = unitData(u).frTrainTrials;
            
            lightInfo = zeros(2,length(binEdges));
            for b = 1:binCount
                [dp 1 u b]
                
                vMax = max(max(responses(1:length(contrasts)*length(orientations))));
                avMax = max(max(responses(1+length(contrasts)*length(orientations):2*length(contrasts)*length(orientations))));
                
                vHist = zeros(length(contrasts),vMax+1);
                avHist = zeros(length(orientations),avMax+1);
                negEntropy = zeros(2,length(contrasts));
                for c = 1:length(contrasts)
                    inds = (c-1)*length(orientations)+1 : c*length(orientations);
                    indsOffset = (inds) + length(orientations)*length(contrasts);
                    
                    [val vMaxOrient] = max(unitData(u).meanResponse(inds,4));
                    [val avMaxOrient] = max(unitData(u).meanResponse(indsOffset,4));
                    
                    vResponses = responses(inds(vMaxOrient),:,b);
                    avResponses = responses(indsOffset(avMaxOrient),:,b);
                    
                    vHist(c,:) = histcounts(vResponses,0:vMax+1,'Normalization','probability');
                    avHist(c,:) = histcounts(avResponses,0:avMax+1,'Normalization','probability');
                    
                    vEnt = vHist(c,:).*log2(vHist(c,:)); vEnt(isnan(vEnt))=0;
                    avEnt = avHist(c,:).*log2(avHist(c,:)); avEnt(isnan(avEnt))=0;
                    
                    negEntropy(1,c) = sum(vEnt);
                    negEntropy(2,c) = sum(avEnt);
                end
                vTotalHist = sum(vHist)/length(contrasts); %divided by length(orientaitons) to make sum = 1 as a probability
                avTotalHist = sum(avHist)/length(contrasts); %divided by length(orientaitons) to make sum = 1 as a probability
                
                vEnt = vTotalHist.*log2(vTotalHist); vEnt(isnan(vEnt))=0;
                avEnt = avTotalHist.*log2(avTotalHist); avEnt(isnan(avEnt))=0;
                
                vTotalEntropy = -sum(vEnt);
                avTotalEntropy = -sum(avEnt);
                
                vInfo = vTotalEntropy - (-sum(1/length(contrasts) * negEntropy(1,:)));
                avInfo = avTotalEntropy - (-sum(1/length(contrasts) * negEntropy(2,:)));
                
                lightInfo(:,b) = [vInfo; avInfo];
            end
            lightInformation{1} = [lightInformation{1}; smooth(lightInfo(1,:))'];
            lightInformation{2} = [lightInformation{2}; smooth(lightInfo(2,:))'];
        end
        
        
        
        %% Orientation information
        if (soundResponsive == ismember(neuronNumber,responsiveUnits.soundResponsiveUnits))...
                && (orientationSelective == ismember(neuronNumber,responsiveUnits.orientationSelectiveUnits))...
                && unitData(u).type~=0
            
            responses = unitData(u).frTrainTrials;
            
            for c = 1:length(contrasts)
                inds = (c-1)*length(orientations)+1 : c*length(orientations);
                indsOffset = (inds) + length(orientations)*length(contrasts);
                
                
                vResponses = responses(inds,:,:);
                avResponses = responses(indsOffset,:,:);
                
                info = zeros(2,binCount);
                for b = 1:binCount
                    [dp 2 u c b]
                    vMax = max(max(vResponses(:,:,b)));
                    avMax = max(max(avResponses(:,:,b)));
                    
                    vHist = zeros(length(orientations),vMax+1);
                    avHist = zeros(length(orientations),avMax+1);
                    negEntropy = zeros(2,length(orientations));
                    for o = 1:length(orientations)
                        vHist(o,:) = histcounts(vResponses(o,:,b),0:vMax+1,'Normalization','probability');
                        avHist(o,:) = histcounts(avResponses(o,:,b),0:avMax+1,'Normalization','probability');
                        
                        vEnt = vHist(o,:).*log2(vHist(o,:)); vEnt(isnan(vEnt))=0;
                        avEnt = avHist(o,:).*log2(avHist(o,:)); avEnt(isnan(avEnt))=0;
                        
                        negEntropy(1,o) = sum(vEnt);
                        negEntropy(2,o) = sum(avEnt);
                    end
                    vTotalHist = sum(vHist)/length(orientations); %divided by length(orientaitons) to make sum = 1 as a probability
                    avTotalHist = sum(avHist)/length(orientations); %divided by length(orientaitons) to make sum = 1 as a probability
                    
                    vEnt = vTotalHist.*log2(vTotalHist); vEnt(isnan(vEnt))=0;
                    avEnt = avTotalHist.*log2(avTotalHist); avEnt(isnan(avEnt))=0;
                    
                    vTotalEntropy = -sum(vEnt);
                    avTotalEntropy = -sum(avEnt);
                    
                    vInfo = vTotalEntropy - (-sum(1/length(orientations) * negEntropy(1,:)));
                    avInfo = avTotalEntropy - (-sum(1/length(orientations) * negEntropy(2,:)));
                    
                    info(:,b) = [vInfo; avInfo];
                end
                orientationInformation{c,1} = [orientationInformation{c,1}; smooth(info(1,:))'];
                orientationInformation{c,2} = [orientationInformation{c,2}; smooth(info(2,:))'];
            end
        end
    end
end

%% Light information figure
f1 = figure;hold on;
vv = mean(lightInformation{1});av = mean(lightInformation{2});
vvSTD = std(lightInformation{1},[],1)./sqrt(size(lightInformation{1},1));
avSTD = std(lightInformation{2},[],1)./sqrt(size(lightInformation{2},1));

hA = area(binEdges,[vv-vvSTD; 2*vvSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
hA = area(binEdges,[av-avSTD; 2*avSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 1];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(binEdges,smooth(mean(lightInformation{1})),'Color',[0 0 0]);
plot(binEdges,smooth(mean(lightInformation{2})),'Color',[0 0 1]);
%     ylim([0 0.2]);
xlabel('Time (ms)');
ylabel('Mutual information_{light} (bits)');
title(['Cont = ' num2str(contrasts(c))]);


%% Orientation information figure
f2 = figure;
for c = 1:length(contrasts)
    subplot(1,length(contrasts),c);hold on;
    
    vv = mean(orientationInformation{c,1});av = mean(orientationInformation{c,2});
    vvSTD = std(orientationInformation{c,1},[],1)./sqrt(size(orientationInformation{c,1},1));
    avSTD = std(orientationInformation{c,2},[],1)./sqrt(size(orientationInformation{c,2},1));
    
    hA = area(binEdges,[vv-vvSTD; 2*vvSTD]');
    hA(1).FaceAlpha = 0;
    hA(1).EdgeColor = [1 1 1];
    hA(2).FaceColor = [0 0 0];
    hA(2).FaceAlpha = 0.5;
    hA(2).EdgeColor = [1 1 1];
    hA = area(binEdges,[av-avSTD; 2*avSTD]');
    hA(1).FaceAlpha = 0;
    hA(1).EdgeColor = [1 1 1];
    hA(2).FaceColor = [0 0 1];
    hA(2).FaceAlpha = 0.5;
    hA(2).EdgeColor = [1 1 1];
    plot(binEdges,smooth(mean(orientationInformation{c,1})),'Color',[0 0 0]);
    plot(binEdges,smooth(mean(orientationInformation{c,2})),'Color',[0 0 1]);
    ylim([0 0.2]);
    xlabel('Time (ms)');
    ylabel('Mutual information_{orientation} (bits)');
    title(['Cont = ' num2str(contrasts(c))]);
end

