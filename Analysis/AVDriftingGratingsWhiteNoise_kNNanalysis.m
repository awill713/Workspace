
clear all;
% close all;

savePath = fullfile('E:\Electrophysiology','EP004');

dataPaths{1} = fullfile('EP004','AW118','20200201');

quantWindow = [50 1000]; %ms after stimulus onset

trainingTrialsPercent = 0.75; % 0 < x < 1, will test on the remaining percent (1-x)
randomizationRepeats = 1000; %repetitions of training/testing blocks, taking the average accuracy
k_value = 35; %nearest neighbors to find

%which neurons to include in the kNN algorithm.
onlySingleUnits = 0;
allNeurons = 0;
lightResponsive = 0;
orientationSelective = 1;

% kk = 30:40;
% pvalue = zeros(1,length(kk));
% maxThing = zeros(1,length(kk));
% for kkk = 1:length(kk)
%     k_value = kk(kkk);

percentVis = zeros(1,randomizationRepeats);
percentVisAud = zeros(1,randomizationRepeats);


for dp = 1:length(dataPaths)
    dp
    %Load data
    dataFile = fullfile('E:\Electrophysiology\',dataPaths{dp},'AVDriftingGratingsWhiteNoise','AVDriftingGratingsWhiteNoiseData.mat');
    load(dataFile);
    stimPath = dir(fullfile('E:\Electrophysiology\',dataPaths{dp},'StimInfo','*AVdriftingGratingsWhiteNoise_stimInfo*'));
    load(fullfile(stimPath.folder,stimPath.name));
    
    quantScalar = 1000/(quantWindow(2)-quantWindow(1));
    quantFirstBin = quantWindow(1)-analysisParams.analysisWindow(1)-analysisParams.frBinWidth+1;
    quantLastBin = quantWindow(2)-analysisParams.analysisWindow(1)-analysisParams.frBinWidth+1;
    
    for r = 1:randomizationRepeats
        r
        orientationCount = stimInfo.orientationCount;
        trainingTrialsCount = round(trainingTrialsPercent * stimInfo.repeats);
        testTrialsCount = stimInfo.repeats - trainingTrialsCount;
        trainingTrials = zeros(orientationCount*2,trainingTrialsCount);
        testTrials = zeros(orientationCount*2,testTrialsCount);
        for t = 1:stimInfo.orientationCount*2
            trainingTrials(t,:) = randperm(stimInfo.repeats,trainingTrialsCount);
            trials = 1:stimInfo.repeats;
            testTrials(t,:) = trials(~ismember(trials,trainingTrials(t,:)));
        end
        
        
        neuronColumnsTraining = [];
        neuronColumnsTest = [];
        for n = 1:length(unitData)
            neuronNumber = unitData(n).neuronNumber;
            if (allNeurons || (lightResponsive && ismember(neuronNumber,responsiveUnits.lightResponsiveUnits))...
                    || (orientationSelective && ismember(neuronNumber,responsiveUnits.orientationSelectiveUnits)))...
                    && (~onlySingleUnits || unitData(n).type==1)
                
                trainingNeuronData = [];
                testNeuronData = [];
                testData = [];
                for t = 1:size(trainingTrials,1)
                    data = mean(unitData(n).frTrainTrials(t,trainingTrials(t,:),quantFirstBin:quantLastBin),3);
                    trainingNeuronData = [trainingNeuronData; data'];
                    testData = [testData; repmat(t,[size(trainingTrials,2) 1]) data'];
                    
                    data = mean(unitData(n).frTrainTrials(t,testTrials(t,:),quantFirstBin:quantLastBin),3);
                    testNeuronData = [testNeuronData; data'];
                end
                
                neuronColumnsTraining = [neuronColumnsTraining trainingNeuronData];
                neuronColumnsTest = [neuronColumnsTest testNeuronData];
            end
        end
        direction = stimInfo.orientations(floor(((1:orientationCount*trainingTrialsCount)-1)./trainingTrialsCount)+1)';
        
        tableVis = array2table(neuronColumnsTraining(1:orientationCount*trainingTrialsCount,:));
        tableVis = [table(direction) tableVis];
        mdlVis = fitcknn(tableVis,'direction','NumNeighbors',k_value,'Standardize',1);
        
        tableVisAud = array2table(neuronColumnsTraining((orientationCount*trainingTrialsCount+1):end,:));
        tableVisAud = [table(direction) tableVisAud];
        mdlVisAud = fitcknn(tableVisAud,'direction','NumNeighbors',k_value,'Standardize',1);
        
        testTableVis = array2table(neuronColumnsTest(1:orientationCount*testTrialsCount,:));
        testTableVisAud = array2table(neuronColumnsTest((orientationCount*testTrialsCount+1):end,:));
        
        correctLabels = stimInfo.orientations(floor(((1:orientationCount*testTrialsCount)-1)./testTrialsCount)+1)';
        visLabels = predict(mdlVis,testTableVis);
        visAudLabels = predict(mdlVisAud,testTableVisAud);
        
        percentVis(r) = sum(visLabels==correctLabels)/length(correctLabels);
        percentVisAud(r) = sum(visAudLabels==correctLabels)/length(correctLabels);
    end
end

f = figure;hold on;histogram(percentVis);histogram(percentVisAud);
title(['Train on ' num2str(trainingTrialsPercent) ', k = ' num2str(k_value) ', vis = ' num2str(mean(percentVis)) ', visAud = ' num2str(mean(percentVisAud))]);
legend({'Light','Light/Sound'});
% saveas(f,fullfile(savePath,'kNN analysis drifting gratings white noise'));
% [h p] = ttest2(percentVis,percentVisAud);
% pvalue(kkk) = p;
% maxThing(kkk) = mean([percentVis percentVisAud]);
% end
% figure;plot(kk,pvalue);
% figure;plot(kk,maxThing);