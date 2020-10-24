
clear all;
% close all;

saveDir = fullfile('E:\Electrophysiology','EP004','MetaData - AVMultiContrastDriftingGratingsWhiteNoise_Orientations','SVM - population');
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
lightResponsive = 1;
orientationSelective = 1;
soundResponsive = 1;


for dp = 1:length(dataPaths)
    dp
    %     [randomize dp]
    %Load data
    dataFile = fullfile('E:\Electrophysiology\',dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise_Orientations','AVMultiContrastDriftingGratingsWhiteNoiseData_Orientations.mat');
    load(dataFile);
    stimPath = dir(fullfile('E:\Electrophysiology\',dataPaths{dp},'StimInfo','*AVmultiContrastDriftingGratingsWhiteNoise_stimInfo*'));
    load(fullfile(stimPath.folder,stimPath.name));
    
    orientations = stimInfo.orientations(1:length(stimInfo.orientations)/2);
    contrasts = stimInfo.contrasts;
    repeats = stimInfo.repeats * 2;
    
    
    if ~exist('performance','var')
        performance = cell(length(contrasts),2);
        %         highPerformance = cell(length(contrasts),2);
        %         lowPerformance = cell(length(contrasts),2);
    end
    
    unitsOfInterest = [];
    for u = 1:length(unitData)
        neuronNumber = unitData(u).neuronNumber;
        if (soundResponsive == ismember(neuronNumber,responsiveUnits.soundResponsiveUnits))...
                && (orientationSelective == ismember(neuronNumber,responsiveUnits.orientationSelectiveUnits))...
                && (lightResponsive == ismember(neuronNumber,responsiveUnits.lightResponsiveUnits))...
                && unitData(u).type~=0
            unitsOfInterest = [unitsOfInterest u];
        end
    end
    
    if length(unitsOfInterest)~=0
        
        for c = 1:length(contrasts)
            ind = (c-1)*length(orientations)+1;
            indOffset = length(orientations)*length(contrasts);
            
            accuracyMat = zeros(length(orientations),length(orientations),2);
            for iDir = 1:length(orientations)
                for jDir = 1:length(orientations)
                    [dp c iDir jDir]
                    
                    accuracy = zeros(2,repeats);
                    for r = 1:repeats
                        iTestTrials = r;
                        iTrainingTrials = setdiff(1:repeats,iTestTrials);
                        jTrainingTrials = 1:repeats;
                        
                        trainingAnswers = [ones(length(iTrainingTrials),1); ones(length(jTrainingTrials),1)*-1];
                        trainingVec = zeros(length(iTrainingTrials)+length(jTrainingTrials),length(unitsOfInterest),1);
                        testingVec = zeros(2,length(unitsOfInterest));
                        for u = 1:length(unitsOfInterest)
                            uu = unitsOfInterest(u);
                            iDataV = unitData(uu).trialResponse(ind+iDir-1,iTrainingTrials);
                            jDataV = unitData(uu).trialResponse(ind+jDir-1,jTrainingTrials);
                            iDataAV = unitData(uu).trialResponse(ind+iDir-1+indOffset,iTrainingTrials);
                            jDataAV = unitData(uu).trialResponse(ind+jDir-1+indOffset,jTrainingTrials);
                            
                                trainingVec(:,u,1) = [iDataV'; jDataV'];
                                testingVec(1,u) = unitData(uu).trialResponse(ind+iDir-1,iTestTrials);
                            
                                trainingVec(:,u,2) = [iDataAV'; jDataAV'];
                                testingVec(2,u) = unitData(uu).trialResponse(ind+iDir-1+indOffset,iTestTrials);
                        end
                        
                        svmVis = fitcsvm(squeeze(trainingVec(:,:,1)),trainingAnswers,'KernelFunction','Linear','Standardize',true,'KernelScale','auto');
                        testPrediction = svmVis.predict(testingVec(1,:));
                        accuracy(1,r) = mean(testPrediction==1); %1 is iDir answer
                        
                        svmAudVis = fitcsvm(squeeze(trainingVec(:,:,2)),trainingAnswers,'KernelFunction','Linear','Standardize',true,'KernelScale','auto');
                        testPrediction = svmAudVis.predict(testingVec(2,:));
                        accuracy(2,r) = mean(testPrediction==1); %1 is iDir answer
                    end
                    
                    accuracyMat(iDir,jDir,1) = nanmean(accuracy(1,:));
                    accuracyMat(iDir,jDir,2) = nanmean(accuracy(2,:));
                end
            end
            performance{c,1} = cat(3,performance{c,1},squeeze(accuracyMat(:,:,1)));
            performance{c,2} = cat(3,performance{c,2},squeeze(accuracyMat(:,:,2)));
        end
    end
end



blueMap = [linspace(0,1,128)' linspace(0,1,128)' ones(128,1)];
redMap = [ones(128,1) linspace(1,0,128)' linspace(1,0,128)'];
map = [blueMap; redMap];
maxVal = 0;
for c = 1:length(contrasts)
    dif = abs(nanmean(performance{c,2},3) - nanmean(performance{c,1},3));
    maxVal = max([maxVal max(dif(:))]);
end
f1 = figure;
set(f1,'Position',[300 150 1000 600]);
for c = 1:length(contrasts)
    vv = nanmean(performance{c,1},3);
    av = nanmean(performance{c,2},3);
    
    subplot(4,length(contrasts),c);
    imagesc(vv,[0 1]);
    colormap(gca,'copper');
    xticks(1:3:6);
    yticks(1:3:6);
    xticklabels(orientations(xticks));
    yticklabels(orientations(yticks));
    title(['C = ' num2str(contrasts(c)) ', vis']);
    
    subplot(4,length(contrasts),c+length(contrasts));
    imagesc(av,[0 1]);
    colormap(gca,'copper');
    xticks(1:3:6);
    yticks(1:3:6);
    xticklabels(orientations(xticks));
    yticklabels(orientations(yticks));
    title(['C = ' num2str(contrasts(c)) ', audiovis']);
    
    subplot(4,length(contrasts),c+2*length(contrasts));
        imagesc(av-vv,[-maxVal maxVal]);
        colormap(gca,map);
    xticks(1:3:6);
    yticks(1:3:6);
    xticklabels(orientations(xticks));
    yticklabels(orientations(yticks));
    title(['C = ' num2str(contrasts(c)) ', difference']);
    
    subplot(4,length(contrasts),c+3*length(contrasts));hold on;
    iMat = logical(eye(length(orientations)));
    scatter(vv(~iMat),av(~iMat));
    line([0 1],[0 1],'Color',[0 0 0]);
    xlim([0 1]);ylim([0 1]);
    delt = mean(av(~iMat)-vv(~iMat));
    [h p] = ttest(vv(~iMat),av(~iMat));
    xlabel('SVM accuracy (visual)');
    ylabel('SVM accuracy (audiovisual)');
    title(['\Delta = ' num2str(delt) ', p = ' num2str(p)]);
end



accuracyMean = zeros(2,length(contrasts));
accuracySTD = zeros(2,length(contrasts));
difference = zeros(1,length(contrasts));
differenSTD = zeros(1,length(contrasts));
for c = 1:length(contrasts)
    acc = zeros(2,size(performance{c,1},3));
    for d = 1:size(performance{c,1},3)
        vvv = performance{c,1}(:,:,d);
        aaa = performance{c,2}(:,:,d);
        iMat = logical(eye(length(orientations)));
        
        acc(1,d) = nanmean(vvv(~iMat));
        acc(2,d) = nanmean(aaa(~iMat));
    end
    accuracyMean(1,c) = mean(acc(1,:));
    accuracyMean(2,c) = mean(acc(2,:));
    accuracySTD(1,c) = std(acc(1,:))/sqrt(size(acc,2));
    accuracySTD(2,c) = std(acc(2,:))/sqrt(size(acc,2));
    
    difference(1,c) = mean(acc(2,:) - acc(1,:));
    differenceSTD(1,c) = std(acc(2,:) - acc(1,:))/sqrt(length(acc(2,:)-acc(1,:)));
end
f2 = figure;hold on;
hA = area(contrasts,[accuracyMean(1,:)-accuracySTD(1,:); 2*accuracySTD(1,:)]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
hA = area(contrasts,[accuracyMean(2,:)-accuracySTD(2,:); 2*accuracySTD(2,:)]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 1];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(contrasts,accuracyMean(1,:),'Color',[0 0 0]);
plot(contrasts,accuracyMean(2,:),'Color',[0 0 1]);
ylim([0.3 0.7]);
xticks(contrasts)
ylabel('SVM accuracy');
xlabel('Contrast');


f3 = figure;hold on;
hA = area(contrasts,[difference-differenceSTD; 2*differenceSTD]')
    hA(1).FaceAlpha = 0;
    hA(1).EdgeColor = [1 1 1];
    hA(2).FaceColor = [0 0 0];
    hA(2).FaceAlpha = 0.5;
    hA(2).EdgeColor = [1 1 1];
plot(contrasts,difference,'Color',[0 0 0]);
xticks(contrasts);
xlabel('Contrast');
ylabel('\Delta SVM accuracy');


