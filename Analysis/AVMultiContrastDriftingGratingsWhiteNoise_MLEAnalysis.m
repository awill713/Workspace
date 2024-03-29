
clear all;
% close all;

saveDir = fullfile('D:\Electrophysiology','EP010','MetaData - AVMultiContrastDriftingGratingsWhiteNoise');
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
dataPaths{9} = fullfile('EP010','AW157','20201212-1');
dataPaths{10} = fullfile('EP010','AW157','20201212-2');
dataPaths{11} = fullfile('EP010','AW158','20201212-1');
dataPaths{12} = fullfile('EP010','AW158','20201212-2');

% dataPaths{13} = fullfile('EP010','AW159','20201213-1');
% dataPaths{14} = fullfile('EP010','AW159','20201213-2');
% dataPaths{15} = fullfile('EP010','AW162','20210102-1');
% dataPaths{16} = fullfile('EP010','AW162','20210102-2');
% dataPaths{17} = fullfile('EP010','AW163','20210102-1');
% dataPaths{18} = fullfile('EP010','AW163','20210102-2');
% dataPaths{19} = fullfile('EP010','AW164','20210105');
% dataPaths{20} = fullfile('EP010','AW165','20210106-1');
% dataPaths{21} = fullfile('EP010','AW165','20210106-2');

%which neurons to include
onlySingleUnits = 0;
allNeurons = 0;
lightResponsive = 1;
orientationSelective = 1;
soundResponsive = 1;


for dp = 1:length(dataPaths)
    %     [randomize dp]
    %Load data
    dataFile = fullfile('D:\Electrophysiology\',dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise_final','AVMultiContrastDriftingGratingsWhiteNoiseData_selectivity.mat');
    load(dataFile);
    stimPath = dir(fullfile('D:\Electrophysiology\',dataPaths{dp},'StimInfo','*AVmultiContrastDriftingGratingsWhiteNoise_stimInfo*'));
    load(fullfile(stimPath(end).folder,stimPath(end).name));
    
    orientations = stimInfo.orientations;
    contrasts = stimInfo.contrasts;
    repeats = stimInfo.repeats;
    
    if ~exist('choiceMap','var')
        choiceMap = cell(2,length(contrasts)); %v and av, by contrasts
        randomMap = cell(2,length(contrasts)); %v and av, by contrasts, trials shuffled
        dpUnitCount = [];
    end
    quantScalar = 1000/(analysisParams.quantWindow(2)-analysisParams.quantWindow(1));
    
    %     analysisWindow = analysisParams.analysisWindow;
    %     quantWindow = analysisParams.quantWindow;
    %     quantWindow = [0 250];
    %     frBinWidth = analysisParams.frBinWidth;
    %     quantFirstBin = quantWindow(1)-analysisWindow(1)-frBinWidth+1;
    %     quantLastBin = quantWindow(2)-analysisWindow(1)-frBinWidth+1;
    
    neuronCount = 0;
    neuronNumberLog = [];
    for u  = 1:size(unitData,2)
        neuronNumber = unitData(u).neuronNumber;
        if (ismember(neuronNumber,responsiveUnitsGLM.lightSoundInteractUnits) ||...
                (ismember(neuronNumber,responsiveUnitsGLM.soundResponsiveUnits) &&...
                ismember(neuronNumber,responsiveUnitsGLM.lightResponsiveUnits))) &&...
                (ismember(neuronNumber,responsiveUnits.orientationSelectiveUnits) ||...
                ismember(neuronNumber,responsiveUnits.directionSelectiveUnits)) &&...
                ismember(neuronNumber,[responsiveUnits.singleUnits responsiveUnits.multiUnits])
            
            neuronCount = neuronCount+1;
            neuronNumberLog = [neuronNumberLog; neuronCount neuronNumber u];
        end
    end
    dpUnitCount = [dpUnitCount neuronCount];
    
    if neuronCount>0
        
        for c = 1:length(contrasts)
            tempChoiceMap = zeros(2,length(orientations),length(orientations));
            
            for o = 1:length(orientations)
                [dp c o]
                mleClassification = zeros(2,length(orientations));
                randomClassification = zeros(2,length(orientations));
                
                
                
                for trial = 1:repeats
                
                    if exist('neuronStats','var')
                        clear neuronStats randStats
                    end
                    probeData = [];
                    
                    for neuron = 1:neuronCount
                        neuronU = neuronNumberLog(neuron,3);
                        
                        for orient = 1:length(orientations)
                            vInd = (c-1)*length(orientations) + orient;
                            avInd = (c-1)*length(orientations) + orient +length(orientations)*length(contrasts);
                            
                            trialsIncluded = 1:repeats;
                            if orient==o
                                trialsIncluded  = setdiff(trialsIncluded,trial);
                                
                                probeTrialV = unitData(neuronU).trialResponse(vInd,trial)/quantScalar;
                                probeTrialAV = unitData(neuronU).trialResponse(avInd,trial)/quantScalar;
                                probeData = [probeData [probeTrialV; probeTrialAV]];
                            end
                            
                            vData = unitData(neuronU).trialResponse(vInd,trialsIncluded)/quantScalar;
                            avData = unitData(neuronU).trialResponse(avInd,trialsIncluded)/quantScalar;
                            
                            neuronStats(neuron,orient,1) = fitdist(vData','Poisson');
                            neuronStats(neuron,orient,2) = fitdist(avData','Poisson');
                        end
                    end
                    
                    vEstimate = maximumLikelihoodFunction(neuronStats(:,:,1),probeData(1,:));
                    mleClassification(1,vEstimate) = mleClassification(1,vEstimate) +1;
                    
                    avEstimate = maximumLikelihoodFunction(neuronStats(:,:,2),probeData(2,:));
                    mleClassification(2,avEstimate) = mleClassification(2,avEstimate) +1;
                    
                end
                
                tempChoiceMap(1,o,:) = mleClassification(1,:) / sum(mleClassification(1,:));
                tempChoiceMap(2,o,:) = mleClassification(2,:) / sum(mleClassification(2,:));
            end
            choiceMap{1,c} = cat(3,choiceMap{1,c},squeeze(tempChoiceMap(1,:,:)));
            size(choiceMap{1,c});
            choiceMap{2,c} = cat(3,choiceMap{2,c},squeeze(tempChoiceMap(2,:,:)));
        end
    end
end

viableDPs = size(choiceMap{1},3);

meanMap = squeeze(mean(choiceMap{1,5},3));
figure;imagesc(meanMap); 
statsMat = [];
Y=[];Sub=[];F1=[];F2=[];FACTNAMES = {'Contrast','Sound'};
for c = 1:5
    clear tempV tempAV tempDiff
    for s=1:size(choiceMap{1,c},3)
        tempV(s) = mean(diag(choiceMap{1,c}(:,:,s)));
        tempAV(s) = mean(diag(choiceMap{2,c}(:,:,s)));
        tempDiff(s) = tempAV(s) - tempV(s);
        
        Y = [Y tempV(s) tempAV(s)];
        Sub = [Sub s s];
        F1 = [F1 c c];
        F2 = [F2 0 1];
    end
    vAcc(c) = mean(tempV);
    vAccStd(c) = std(tempV)/sqrt(length(tempV));
    avAcc(c) = mean(tempAV);
    avAccStd(c) = std(tempAV)/sqrt(length(tempAV));
    diffAcc(c) = mean(tempDiff);
    diffAccStd(c) = std(diffAcc)/sqrt(length(diffAcc));
    diffAccStd(c) = std(tempDiff)/sqrt(length(tempDiff));
    
    statsMat = [statsMat; [tempV' tempAV']];
    
%     [c vAcc avAcc]
end
figure;hold on;
hA = area(contrasts,[vAcc-vAccStd; 2*vAccStd]')
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
hA = area(contrasts,[avAcc-avAccStd; 2*avAccStd]')
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 1];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];

plot(contrasts,vAcc,'Color',[0 0 0]);
plot(contrasts,avAcc,'Color',[0 0 1]);

figure;hold on;
hA = area(contrasts,[diffAcc-diffAccStd; 2*diffAccStd]')
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(contrasts,diffAcc,'Color',[0 0 0]);

[p,~,~] = anova2(statsMat,viableDPs,'off');
stats = rm_anova2(Y,Sub,F1,F2,FACTNAMES)

choiceTrace = zeros(2,length(contrasts),length(orientations)+1,size(choiceMap{1},3));
accTrace = zeros(2,length(contrasts),length(orientations)+1);
accSTD = zeros(2,length(contrasts),length(orientations)+1);
for c = 1:5
    for av = 1:2
        for dp = 1:size(choiceTrace,4)
            tempMap = choiceMap{av,c}(:,:,dp);
            collapsed = zeros(length(orientations));
            for o = 1:length(orientations)
                shift = 7-o;
                oTrace = tempMap(o,:);
                oTrace = circshift(oTrace,shift);
                collapsed(o,:) = oTrace;
            end
            collapsed = mean(collapsed);
            collapsed = [collapsed, collapsed(1)];
            choiceTrace(av,c,:,dp) = collapsed;
        end
    end
    
    accTrace(1,c,:) = squeeze(mean(choiceTrace(1,c,:,:),4));
    accSTD(1,c,:) = squeeze(std(choiceTrace(1,c,:,:),[],4))/sqrt(size(choiceTrace,4));
    accTrace(2,c,:) = squeeze(mean(choiceTrace(2,c,:,:),4));
    accSTD(2,c,:) = squeeze(std(choiceTrace(2,c,:,:),[],4))/sqrt(size(choiceTrace,4));
end

xs = -180:30:180;
figure;
for c = 1:5
    vv = squeeze((accTrace(1,c,:)))';
    vvSTD = squeeze((accSTD(1,c,:)))';
    aavv = squeeze((accTrace(2,c,:)))';
    aavvSTD = squeeze((accSTD(2,c,:)))';
    
    subplot(1,5,c);hold on;
    errorbar(xs,vv,vvSTD,'Color',[0 0 0],'LineWidth',1);
    errorbar(xs,aavv,aavvSTD,'Color',[0 0 1],'LineWidth',1);
    
    ylim([0 0.3]);
    title(['Contrast = ' num2str(contrasts(c))]);
    
    hA = area(xs,[vv-vvSTD; 2*vvSTD]')
    hA(1).FaceAlpha = 0;
    hA(1).EdgeColor = [1 1 1];
    hA(2).FaceColor = [0 0 0];
    hA(2).FaceAlpha = 0.5;
    hA(2).EdgeColor = [1 1 1];
    hA = area(xs,[aavv-aavvSTD; 2*aavvSTD]')
    hA(1).FaceAlpha = 0;
    hA(1).EdgeColor = [1 1 1];
    hA(2).FaceColor = [0 0 1];
    hA(2).FaceAlpha = 0.5;
    hA(2).EdgeColor = [1 1 1];
    
    plot(xs,vv,'Color',[0 0 0]);
    plot(xs,aavv,'Color',[0 0 1]);
    
end

