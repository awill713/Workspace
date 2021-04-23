
clear all;
% close all;

saveDir = fullfile('D:\Electrophysiology','EP010','MetaData - AVMultiContrastDriftingGratingsWhiteNoise');
if ~exist(saveDir)
    mkdir(saveDir);
end

% dataPaths{1} = fullfile('EP010','AW159','20201213-1');
% dataPaths{2} = fullfile('EP010','AW159','20201213-2');
% dataPaths{3} = fullfile('EP010','AW162','20210102-1');
% dataPaths{4} = fullfile('EP010','AW162','20210102-2');
% dataPaths{5} = fullfile('EP010','AW163','20210102-1');
% dataPaths{6} = fullfile('EP010','AW163','20210102-2');
% dataPaths{7} = fullfile('EP010','AW164','20210105');
dataPaths{1} = fullfile('EP010','AW165','20210106-1');
% dataPaths{9} = fullfile('EP010','AW165','20210106-2');

%which neurons to include
onlySingleUnits = 0;
allNeurons = 0;
lightResponsive = 1;
orientationSelective = 1;
soundResponsive = 1;

infoBins = 10;
lambda = cell(1,2);

for dp = 1:length(dataPaths)
    dp
    %Load data
    dataFile = fullfile('D:\Electrophysiology\',dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise_Video_NEW_BSR_25bin','AVMultiContrastDriftingGratingsWhiteNoiseData.mat');
    load(dataFile);
    stimPath = dir(fullfile('D:\Electrophysiology\',dataPaths{dp},'StimInfo','*AVmultiContrastDriftingGratingsWhiteNoise_stimInfo*'));
    load(fullfile(stimPath(end).folder,stimPath(end).name));
    
    orientations = stimInfo.orientations;
    contrasts = stimInfo.contrasts;
    
    binCount = analysisParams.analysisWindow(2)-analysisParams.analysisWindow(1)-analysisParams.frBinWidth+1; %sliding window
    binEdges = analysisParams.analysisWindow(1)+analysisParams.frBinWidth:1:analysisParams.analysisWindow(2);
    
    if ~exist('orientationInformation','var')
        orientationInformation = cell(length(contrasts),3);
        lightInformation = cell(1,2);
    end
    
    for u  = 1:size(unitData,2)
        neuronNumber = unitData(u).neuronNumber;
        %% Light information
        %          if (((ismember(neuronNumber,responsiveUnitsGLM.lightResponsiveUnits) &&...
        %                 ismember(neuronNumber,responsiveUnitsGLM.soundResponsiveUnits))) ||...
        %                 ismember(neuronNumber,responsiveUnitsGLM.lightSoundInteractUnits))
        if ismember(neuronNumber,responsiveUnitsGLM.lightSoundInteractUnits)
            
            responses = unitData(u).frTrainTrials;
            
            lightInfo = zeros(2,length(binEdges));
            for b = 1:binCount
                [dp 1 u b]
                
                vMax = max(max(responses(1:length(contrasts)*length(orientations))));
                avMax = max(max(responses(1+length(contrasts)*length(orientations):2*length(contrasts)*length(orientations))));
                vMin = min(min(responses(1:length(contrasts)*length(orientations))));
                avMin = min(min(responses(1+length(contrasts)*length(orientations):2*length(contrasts)*length(orientations))));
                theMax = max([vMax avMax]);
                theMin = min([vMin avMin]);
                
                vHist = zeros(length(contrasts),infoBins);
                avHist = zeros(length(orientations),infoBins);
                negEntropy = zeros(2,length(contrasts));
                for c = 1:length(contrasts)
                    inds = (c-1)*length(orientations)+1 : c*length(orientations);
                    indsOffset = (inds) + length(orientations)*length(contrasts);
                    
                    [val vMaxOrient] = max(unitData(u).meanResponse(inds,4));
                    [val avMaxOrient] = max(unitData(u).meanResponse(indsOffset,4));
                    
                    vResponses = responses(inds(vMaxOrient),:,b);
                    avResponses = responses(indsOffset(avMaxOrient),:,b);
                    
                    vHist(c,:) = histcounts(vResponses,linspace(theMin,theMax,infoBins+1),'Normalization','probability');
                    avHist(c,:) = histcounts(avResponses,linspace(theMin,theMax,infoBins+1),'Normalization','probability');
                    
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
        %         if (((ismember(neuronNumber,responsiveUnitsGLM.lightResponsiveUnits) &&...
        %                 ismember(neuronNumber,responsiveUnitsGLM.soundResponsiveUnits))) ||...
        %                 ismember(neuronNumber,responsiveUnitsGLM.lightSoundInteractUnits)) &&...
        %                 ismember(neuronNumber,responsiveUnits.orientationSelectiveUnits)
        if ismember(neuronNumber,responsiveUnitsGLM.lightSoundInteractUnits) &&...
                ismember(neuronNumber,responsiveUnits.orientationSelectiveUnits)
            
            responses = unitData(u).frTrainTrials;
            
            for c = 1:length(contrasts)
                inds = (c-1)*length(orientations)+1 : c*length(orientations);
                indsOffset = (inds) + length(orientations)*length(contrasts);
                
                
                vResponses = responses(inds,:,:);
                avResponses = responses(indsOffset,:,:);
                
                info = zeros(2,binCount);
                tempLambda = zeros(2,binCount);
                for b = 1:binCount
                    [dp 2 u c b]
                    vMax = max(max(vResponses(:,:,b)));
                    avMax = max(max(avResponses(:,:,b)));
                    vMin = min(min(vResponses(:,:,b)));
                    avMin = min(min(avResponses(:,:,b)));
                    theMax = max([vMax avMax]);
                    theMin = min([vMin avMin]);
%                     vMax = max(max(mean(vResponses(:,:,b-9:b),3)));
%                     avMax = max(max(mean(avResponses(:,:,b-9:b),3)));

                    tempLambda(:,b) = [poissfit(reshape(vResponses(:,:,b),[1 120]));...
                        poissfit(reshape(avResponses(:,:,b),[1 120]))];
                    
                    vHist = zeros(length(orientations),infoBins);
                    avHist = zeros(length(orientations),infoBins);
                    negEntropy = zeros(2,length(orientations));
                    for o = 1:length(orientations)
                        vHist(o,:) = histcounts(vResponses(o,:,b),linspace(theMin,theMax,infoBins+1),'Normalization','probability');
                        avHist(o,:) = histcounts(avResponses(o,:,b),linspace(theMin,theMax,infoBins+1),'Normalization','probability');
%                         vHist(o,:) = histcounts(mean(vResponses(o,:,b-9:b),3),0:vMax+1,'Normalization','probability');
%                         avHist(o,:) = histcounts(mean(avResponses(o,:,b-9:b),3),0:avMax+1,'Normalization','probability');
                        
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
                orientationInformation{c,3} = [orientationInformation{c,3}; [dp u neuronNumber]];
                
                lambda{1,1} = [lambda{1,1}; tempLambda(1,:)];
                lambda{1,2} = [lambda{1,2}; tempLambda(2,:)];
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
%     ylim([0.05 0.22]);
    xlabel('Time (ms)');
    ylabel('Mutual information_{orientation} (bits)');
    title(['Cont = ' num2str(contrasts(c))]);
end

%% 
f3 = figure;
for c = 1:length(contrasts)
    subplot(1,length(contrasts),c);hold on;
    
    vv = mean(orientationInformation{c,1}-mean(orientationInformation{1,1},2));
    av = mean(orientationInformation{c,2}-mean(orientationInformation{1,1},2));
    vv = smooth(vv,70);
    av = smooth(av,70);
    plot(binEdges,vv,'Color',[0 0 0]);
    plot(binEdges,av,'Color',[0 0 1]);
    ylim([-.02 0.1]);
    xlabel('Time (ms)');
    ylabel('Mutual information_{orientation} (bits)');
    title(['Cont = ' num2str(contrasts(c))]);
end