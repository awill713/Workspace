
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
lightResponsive = 0;
orientationSelective = 0;
soundResponsive = 1;

trainingTrialsPercent = 0.8;
randomizations = 5;

for dp = 1:length(dataPaths)
    %     [randomize dp]
    %Load data
    dataFile = fullfile('E:\Electrophysiology\',dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise','AVMultiContrastDriftingGratingsWhiteNoiseData.mat');
    load(dataFile);
    stimPath = dir(fullfile('E:\Electrophysiology\',dataPaths{dp},'StimInfo','*AVmultiContrastDriftingGratingsWhiteNoise_stimInfo*'));
    load(fullfile(stimPath.folder,stimPath.name));
    
    orientations = stimInfo.orientations;
    contrasts = stimInfo.contrasts;
    repeats = stimInfo.repeats;
    
    if ~exist('bayesianDecoding','var')
        bayesianDecoding = zeros(0,randomizations,length(contrasts),2);
        randomDecoding = zeros(0,randomizations,length(contrasts),2);
%         dpCount = 0;
    end
    
    analysisWindow = analysisParams.analysisWindow;
    quantWindow = analysisParams.quantWindow;
    frBinWidth = analysisParams.frBinWidth;
    quantFirstBin = quantWindow(1)-analysisWindow(1)-frBinWidth+1;
    quantLastBin = quantWindow(2)-analysisWindow(1)-frBinWidth+1;
    
    
    for randomize = 1:randomizations
        
        if exist('neuronStats','var')
           clear neuronStats randStats
        end
        
        neuronCount = 0;
        neuronNumberLog = [];
        
        [dp randomize]
        
        trainingTrials = randperm(repeats,repeats*trainingTrialsPercent);
%                 trainingTrials = 1:10;
        allTrials = 1:repeats;
        testTrials = allTrials(~ismember(allTrials,trainingTrials));
%                 testTrials = 1:10;
%                 testTrials = randperm(10,2);
%                 trainingTrials = randperm(10,8);
                allTrials = 1:repeats;
%                 testTrials = allTrials(~ismember(allTrials,trainingTrials));
%                 testTrials = randperm(10,9);
        for u  = 1:size(unitData,2)
            neuronNumber = unitData(u).neuronNumber;
            if (allNeurons || (lightResponsive && ismember(neuronNumber,responsiveUnits.lightResponsiveUnits))...
                    || (orientationSelective && ismember(neuronNumber,responsiveUnits.orientationSelectiveUnits))...
                    || (soundResponsive && ismember(neuronNumber,responsiveUnits.soundResponsiveUnits)))...
                    && (~onlySingleUnits || unitData(u).type==1) && unitData(u).type~=0
                
                neuronCount = neuronCount+1;
                neuronNumberLog = [neuronNumberLog; neuronCount u];
                %             neuronNumber
                
                trialResponses = mean(unitData(u).frTrainTrials(:,:,quantFirstBin:quantLastBin),3);
                
                for c = 1:length(contrasts)
                    
                    randOrientTrials = zeros(length(orientations),length(trainingTrials));
                    potentialTrials = 1:length(orientations)*length(trainingTrials);
                    for oh = 1:length(orientations)
                        randOrientTrials(oh,:) = potentialTrials(randperm(length(potentialTrials),length(trainingTrials)));
                        potentialTrials = potentialTrials(~ismember(potentialTrials,randOrientTrials(oh,:)));
                    end
                    
                    vTrialResponses = trialResponses((c-1)*length(orientations)+1:(c-1)*length(orientations)+length(orientations),trainingTrials);
                    vTrialResponses = vTrialResponses(:);
                    avTrialResponses = trialResponses((c-1)*length(orientations)+1+length(contrasts)*length(orientations):(c-1)*length(orientations)+length(orientations)+length(contrasts)*length(orientations),trainingTrials);
                    avTrialResponses = avTrialResponses(:);
                    for ori = 1:length(orientations)
                        ind = (c-1)*length(orientations)+ori;
                        
                        vTrials = trialResponses(ind,trainingTrials);
                        pdV = fitdist(vTrials','Normal');
                        neuronStats(neuronCount,c,ori,1) = pdV;
                        
                        avTrials = trialResponses(ind+length(orientations)*length(contrasts),trainingTrials);
                        pdAV = fitdist(avTrials','Normal');
                        neuronStats(neuronCount,c,ori,2) = pdAV;
                        
                        randVTrials = vTrialResponses(randOrientTrials(ori,:));
                        randAVTrials = avTrialResponses(randOrientTrials(ori,:));
                        
                        pdRandV = fitdist(randVTrials,'Normal');
                        randStats(neuronCount,c,ori,1) = pdRandV;
                        pdRandAV = fitdist(randAVTrials,'Normal');
                        randStats(neuronCount,c,ori,2) = pdRandAV;
                    end
                end
            end
        end
        
%         if neuronCount~=0
%             dpCount = dpCount+1;
            
            for c = 1:length(contrasts)
                [dp randomize c]
                actual = zeros(1,length(orientations)*length(testTrials));
                realPredict = zeros(2,size(actual,2));
                randPredict = zeros(2,size(actual,2));
                
                trialCount = 0;
                for ori = 1:length(orientations)
                    for t = 1:length(testTrials)
                        probeTrial = testTrials(t);
                        trialCount = trialCount+1;
                        
                        realProbs = ones(length(orientations),2);
                        randProbs = ones(length(orientations),2);
                        
                        for n = 1:size(neuronStats,1)
                            u = neuronNumberLog(n,2);
                            
                            respV = mean(unitData(u).frTrainTrials((c-1)*length(orientations)+ori,probeTrial,quantFirstBin:quantLastBin));
                            respAV = mean(unitData(u).frTrainTrials((c-1)*length(orientations)+ori+length(orientations)*length(contrasts),probeTrial,quantFirstBin:quantLastBin));
                            for o = 1:length(orientations)
                                if neuronStats(n,c,o,1).sigma~=0
                                    prob = pdf(neuronStats(n,c,o,1),respV);
                                    realProbs(o,1) = realProbs(o,1)*prob;
                                end
                                
                                if neuronStats(n,c,o,2).sigma~=0
                                    prob = pdf(neuronStats(n,c,o,2),respAV);
                                    realProbs(o,2) = realProbs(o,2)*prob;
                                end
                                
                                if randStats(n,c,o,1).sigma~=0
                                    prob = pdf(randStats(n,c,o,1),respV);
                                    randProbs(o,1) = randProbs(o,1)*prob;
                                end
                                
                                if randStats(n,c,o,2).sigma~=0
                                    prob = pdf(randStats(n,c,o,2),respAV);
                                    randProbs(o,2) = randProbs(o,2)*prob;
                                end
                            end
                        end
                        
                        [val ind] = max(realProbs(:,1));
                        actual(1,trialCount) = ori;
                        realPredict(1,trialCount) = ind;
                        
                        [val ind] = max(realProbs(:,2));
                        realPredict(2,trialCount) = ind;
                        
                        [val ind] = max(randProbs(:,1));
                        randPredict(1,trialCount) = ind;
                        
                        [val ind] = max(randProbs(:,2));
                        randPredict(2,trialCount) = ind;
                    end
                end
                
                bayesianDecoding(dp,randomize,c,1) = sum(actual==realPredict(1,:))/length(actual);
                bayesianDecoding(dp,randomize,c,2) = sum(actual==realPredict(2,:))/length(actual);
                randomDecoding(dp,randomize,c,1) = sum(actual==randPredict(1,:))/length(actual);
                randomDecoding(dp,randomize,c,2) = sum(actual==randPredict(2,:))/length(actual);
            end
%         end
        
%         clear neuronStats randStats;
        
    end
end

meanBayes = squeeze(mean(bayesianDecoding,2));
meanRand = squeeze(mean(randomDecoding,2));
% meanBayes = bayesianDecoding;
% meanRand = randomDecoding;

anova2([squeeze(meanBayes(:,:,1)); squeeze(meanBayes(:,:,2))],8)


f1 = figure;hold on;

semV = std(meanBayes(:,:,1))/sqrt(size(meanBayes,1));
semVA = std(meanBayes(:,:,2))/sqrt(size(meanBayes,1));
hA = area(contrasts,[mean(meanBayes(:,:,1))-semV;2*semV]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.4;
hA(2).EdgeColor = [1 1 1];
hAV = area(contrasts,[mean(meanBayes(:,:,2))-semVA;2*semVA]');
hAV(1).FaceAlpha = 0;
hAV(1).EdgeColor = [1 1 1];
hAV(2).FaceColor = [0 0 1];
hAV(2).FaceAlpha = 0.4;
hAV(2).EdgeColor = [1 1 1];
semVRand = std(meanRand(:,:,1))/sqrt(size(meanRand,1));
semVARand = std(meanRand(:,:,2))/sqrt(size(meanRand,1));
hA = area(contrasts,[mean(meanRand(:,:,1))-semVRand;2*semVRand]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.4;
hA(2).EdgeColor = [1 1 1];
hAV = area(contrasts,[mean(meanRand(:,:,2))-semVARand;2*semVARand]');
hAV(1).FaceAlpha = 0;
hAV(1).EdgeColor = [1 1 1];
hAV(2).FaceColor = [0 0 1];
hAV(2).FaceAlpha = 0.4;
hAV(2).EdgeColor = [1 1 1];


plot(contrasts,mean(meanBayes(:,:,1)),'Color',[0 0 0],'LineWidth',2);
plot(contrasts,mean(meanBayes(:,:,2)),'Color',[0 0 1],'LineWidth',2);
plot(contrasts,mean(meanRand(:,:,1)),'Color',[0 0 0],'LineWidth',2,'LineStyle',':');
plot(contrasts,mean(meanRand(:,:,2)),'Color',[0 0 1],'LineWidth',2,'LineStyle',':');

ylabel('Accuracy');
xlabel('Contrast');
xticks(contrasts);
suptitle({'Bayesian maximum likelihood decoding accuracy'});
% saveas(f1,fullfile(saveDir,'Maximum likelihood decoding accuracy - sound responsive'));
