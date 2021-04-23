
clear all;
% close all;

saveDir = fullfile('D:\Electrophysiology','EP010','MetaData - AVMultiContrastDriftingGratingsWhiteNoise');
if ~exist(saveDir)
    mkdir(saveDir);
end

dataPaths{1} = fullfile('EP010','AW159','20201213-1');
dataPaths{2} = fullfile('EP010','AW159','20201213-2');
dataPaths{3} = fullfile('EP010','AW162','20210102-1');
dataPaths{4} = fullfile('EP010','AW162','20210102-2');
dataPaths{5} = fullfile('EP010','AW163','20210102-1');
dataPaths{6} = fullfile('EP010','AW163','20210102-2');
dataPaths{7} = fullfile('EP010','AW164','20210105');
dataPaths{8} = fullfile('EP010','AW165','20210106-1');
dataPaths{9} = fullfile('EP010','AW165','20210106-2');

%which neurons to include
onlySingleUnits = 0;
allNeurons = 0;
lightResponsive = 1;
orientationSelective = 1;
soundResponsive = 1;


for dp = 1:length(dataPaths)
    dp
    %Load data
    dataFile = fullfile('D:\Electrophysiology\',dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise_Video_NEW','AVMultiContrastDriftingGratingsWhiteNoiseData.mat');
    load(dataFile);
    stimPath = dir(fullfile('D:\Electrophysiology\',dataPaths{dp},'StimInfo','*AVmultiContrastDriftingGratingsWhiteNoise_stimInfo*'));
    load(fullfile(stimPath.folder,stimPath.name));
    
    orientations = stimInfo.orientations;
    contrasts = stimInfo.contrasts;
    
    if ~exist('orientationInformation','var')
        orientationInformation = cell(1,3); %first cell is v info across contrasts,second is av info, third cell is unit info
    end
    
    for u  = 1:size(unitData,2)
        neuronNumber = unitData(u).neuronNumber;
        
        %% Orientation information
        if (soundResponsive == ismember(neuronNumber,responsiveUnits.soundResponsiveUnits))...
                && (orientationSelective == ismember(neuronNumber,responsiveUnits.orientationSelectiveUnits))...
                && (lightResponsive == ismember(neuronNumber,responsiveUnits.lightResponsiveUnits))...
                && unitData(u).type~=0
            
            responses = unitData(u).trialResponse;
            
            infoContrasts = zeros(2,length(contrasts));
            for c = 1:length(contrasts)
                inds = (c-1)*length(orientations)+1 : c*length(orientations);
                indsOffset = (inds) + length(orientations)*length(contrasts);
                
                
                vResponses = responses(inds,:);
                avResponses = responses(indsOffset,:);
                
                %                 info = zeros(2,binCount);
                
                vMax = ceil(max(max(vResponses)))+1;
                avMax = ceil(max(max(avResponses)))+1;
                
                vHist = zeros(length(orientations),vMax);
                avHist = zeros(length(orientations),avMax);
                negEntropy = zeros(2,length(orientations));
                for o = 1:length(orientations)
                    vHist(o,:) = histcounts(vResponses(o,:),0:vMax,'Normalization','probability');
                    avHist(o,:) = histcounts(avResponses(o,:),0:avMax,'Normalization','probability');
                    
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
                
                infoContrasts(:,c) = [vInfo; avInfo];
            end
            
            orientationInformation{1,1} = [orientationInformation{1,1}; infoContrasts(1,:)];
            orientationInformation{1,2} = [orientationInformation{1,2}; infoContrasts(2,:)];
            orientationInformation{1,3} = [orientationInformation{1,3}; dp u neuronNumber];
        end
    end
end


%% Orientation information figure

f1 = figure; hold on;

vv = mean(orientationInformation{1,1});av = mean(orientationInformation{1,2});
vvSTD = std(orientationInformation{1,1},[],1)./sqrt(size(orientationInformation{1,1},1));
avSTD = std(orientationInformation{1,2},[],1)./sqrt(size(orientationInformation{1,2},1));

hA = area(contrasts,[vv-vvSTD; 2*vvSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
hA = area(contrasts,[av-avSTD; 2*avSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 1];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(contrasts,vv,'Color',[0 0 0]);
plot(contrasts,av,'Color',[0 0 1]);
%     ylim([0 0.2]);
xlabel('Contrast');
ylabel('Mutual information_{orientation} (bits)');

f2 = figure; hold on;

diffInfo = orientationInformation{1,2} - orientationInformation{1,1};
dd = mean(diffInfo);
ddSTD = std(diffInfo,[],1)./sqrt(size(diffInfo,1));

hA = area(contrasts,[dd-ddSTD; 2*ddSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(contrasts,dd,'Color',[0 0 0]);
%     ylim([0 0.2]);
xlabel('Contrast');
ylabel('\Delta Mutual information_{orientation} (bits)');
