
clear all;
% close all;

saveDir = fullfile('D:\Electrophysiology','EP004','MetaData - AVMultiContrastDriftingGratingsWhiteNoise','MLE model','Individual neurons - recat');
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
    
    %Load data
    dataFile = fullfile('D:\Electrophysiology\',dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise_Recat','AVMultiContrastDriftingGratingsWhiteNoiseData.mat');
    load(dataFile);
    stimPath = dir(fullfile('D:\Electrophysiology\',dataPaths{dp},'StimInfo','*AVmultiContrastDriftingGratingsWhiteNoise_stimInfo*'));
    load(fullfile(stimPath.folder,stimPath.name));
    
    orientations = stimInfo.orientations;
    contrasts = stimInfo.contrasts;
    repeats = stimInfo.repeats;
    
    if ~exist('choiceMap','var')
        choiceMap = cell(2,length(contrasts)); %v and av, by contrasts
        accuracy = cell(2,length(contrasts));
        orientationMap = cell(1,length(contrasts));
        directionMap = cell(1,length(contrasts));
    end
    
    for u = 1:size(unitData,2)
        neuronNumber = unitData(u).neuronNumber;
        if ismember(neuronNumber, responsiveUnits.soundResponsiveUnits)...
                && (ismember(neuronNumber, responsiveUnits2.directionSelectiveUnits)...
                || ismember(neuronNumber, responsiveUnits2.orientationSelectiveUnits))
            
            
            
            for c = 1:length(contrasts)
                
                %%Overall accuracy
                tempChoiceMap = zeros(2,length(orientations),length(orientations));
                for testDir = 1:length(orientations)
                    [dp u c testDir]
                    
                    mleClassification = zeros(2,length(orientations));
                    
                    
                    for trial = 1:repeats
                        if exist('neuronStats','var')
                            clear neuronStats
                        end
                        
                        
                        for orient = 1:length(orientations)
                            
                            vInd = (c-1)*length(orientations) + orient;
                            avInd = (c-1)*length(orientations) + orient + length(orientations)*length(contrasts);
                            
                            trialsIncluded = 1:repeats;
                            if orient==testDir
                                trialsIncluded  = setdiff(trialsIncluded,trial);
                                
                                probeTrialV = unitData(u).trialResponse(vInd,trial);
                                probeTrialAV = unitData(u).trialResponse(avInd,trial);
                                probeData = [probeTrialV; probeTrialAV];
                            end
                            
                            vData = unitData(u).trialResponse(vInd,trialsIncluded);
                            avData = unitData(u).trialResponse(avInd,trialsIncluded);
                            
                            neuronStats(orient,1) = fitdist(vData','Normal');
                            neuronStats(orient,2) = fitdist(avData','Normal');
                        end
                        
                        vEstimate = maximumLikelihoodFunctionSingle(neuronStats(:,1),probeData(1));
                        mleClassification(1,vEstimate) = mleClassification(1,vEstimate)+1;
                        
                        avEstimate = maximumLikelihoodFunctionSingle(neuronStats(:,2),probeData(2));
                        mleClassification(2,avEstimate) = mleClassification(2,avEstimate) +1;
                    end
                    
                    tempChoiceMap(1,testDir,:) = mleClassification(1,:) / sum(mleClassification(1,:));
                    tempChoiceMap(2,testDir,:) = mleClassification(2,:) / sum(mleClassification(2,:));
                end
                choiceMap{1,c} = cat(3,choiceMap{1,c},squeeze(tempChoiceMap(1,:,:)));
                choiceMap{2,c} = cat(3,choiceMap{2,c},squeeze(tempChoiceMap(2,:,:)));
                
                accuracy{1,c} = [accuracy{1,c}; mean(diag(squeeze(tempChoiceMap(1,:,:)))) dp u neuronNumber size(accuracy{1,c},1)+1];
                accuracy{2,c} = [accuracy{2,c}; mean(diag(squeeze(tempChoiceMap(2,:,:)))) dp u neuronNumber size(accuracy{2,c},1)+1];
                
                %%Orientation accuracy
                vIndBase = (c-1)*length(orientations) + 1;
                avIndBase = (c-1)*length(orientations) + 1 + length(orientations)*length(contrasts);
                
                vResp = unitData(u).meanResponse(vIndBase:vIndBase+length(orientations)-1,4);
                avResp = unitData(u).meanResponse(avIndBase:avIndBase+length(orientations)-1,4);
                [~, prefOrientV] = max(vResp);
                [~, prefOrientAV] = max(avResp);
                
                orthIndV = [mod(prefOrientV+3-1,length(orientations)) mod(prefOrientV-3-1,length(orientations))]+1;
                orthIndAV = [mod(prefOrientAV+3-1,length(orientations)) mod(prefOrientAV-3-1,length(orientations))]+1;
                
                orientClassification = zeros(2,2);
                for trial = 1:repeats
                    if exist('neuronStats','var')
                        clear neuronStats
                    end
                    vIndBase = (c-1)*length(orientations);
                    avIndBase = (c-1)*length(orientations) + length(orientations)*length(contrasts);
                    
                    probeTrialV = unitData(u).trialResponse(vIndBase + prefOrientV,trial);
                    probeTrialAV = unitData(u).trialResponse(avIndBase +prefOrientAV,trial);
                    probeData = [probeTrialV; probeTrialAV];
                    
                    prefDataV = unitData(u).trialResponse(vIndBase+prefOrientV,setdiff(1:repeats,trial));
                    prefDataAV = unitData(u).trialResponse(avIndBase+prefOrientAV,setdiff(1:repeats,trial));
                    neuronStats(1,1) = fitdist(prefDataV','Normal');
                    neuronStats(1,2) = fitdist(prefDataAV','Normal');
                    
                    orthDataV = reshape(unitData(u).trialResponse(vIndBase+orthIndV,:),[1 2*repeats]);
                    orthDataAV = reshape(unitData(u).trialResponse(avIndBase+orthIndAV,:),[1 2*repeats]);
                    neuronStats(2,1) = fitdist(orthDataV','Normal');
                    neuronStats(2,2) = fitdist(orthDataAV','Normal');
                    
                    vEstimate = maximumLikelihoodFunctionSingle(neuronStats(:,1),probeData(1));
                    orientClassification(1,vEstimate) = orientClassification(1,vEstimate)+1;
                    
                    avEstimate = maximumLikelihoodFunctionSingle(neuronStats(:,2),probeData(2));
                    orientClassification(2,avEstimate) = orientClassification(2,avEstimate)+1;
                end
                if ismember(neuronNumber, responsiveUnits2.orientationSelectiveUnits)
                    orientationMap{1,c} = [orientationMap{1,c}; orientClassification(:,1)'./repeats dp u neuronNumber];
                end
                
                %%Direction accuracy
                vIndBase = (c-1)*length(orientations) + 1;
                avIndBase = (c-1)*length(orientations) + 1 + length(orientations)*length(contrasts);
                
                vResp = unitData(u).meanResponse(vIndBase:vIndBase+length(orientations)-1,4);
                avResp = unitData(u).meanResponse(avIndBase:avIndBase+length(orientations)-1,4);
                [~, prefOrientV] = max(vResp);
                [~, prefOrientAV] = max(avResp);
                
                oppIndV = mod(prefOrientV+6-1,length(orientations))+1;
                oppIndAV = mod(prefOrientAV+6-1,length(orientations))+1;
                
                directClassification = zeros(2,2);
                for trial = 1:repeats
                    if exist('neuronStats','var')
                        clear neuronStats
                    end
                    vIndBase = (c-1)*length(orientations);
                    avIndBase = (c-1)*length(orientations) + length(orientations)*length(contrasts);
                    
                    probeTrialV = unitData(u).trialResponse(vIndBase + prefOrientV,trial);
                    probeTrialAV = unitData(u).trialResponse(avIndBase +prefOrientAV,trial);
                    probeData = [probeTrialV; probeTrialAV];
                    
                    prefDataV = unitData(u).trialResponse(vIndBase+prefOrientV,setdiff(1:repeats,trial));
                    prefDataAV = unitData(u).trialResponse(avIndBase+prefOrientAV,setdiff(1:repeats,trial));
                    neuronStats(1,1) = fitdist(prefDataV','Normal');
                    neuronStats(1,2) = fitdist(prefDataAV','Normal');
                    
                    oppDataV = unitData(u).trialResponse(vIndBase+oppIndV,:);
                    oppDataAV = unitData(u).trialResponse(avIndBase+oppIndAV,:);
                    neuronStats(2,1) = fitdist(oppDataV','Normal');
                    neuronStats(2,2) = fitdist(oppDataAV','Normal');
                    
                    vEstimate = maximumLikelihoodFunctionSingle(neuronStats(:,1),probeData(1));
                    directClassification(1,vEstimate) = directClassification(1,vEstimate)+1;
                    
                    avEstimate = maximumLikelihoodFunctionSingle(neuronStats(:,2),probeData(2));
                    directClassification(2,avEstimate) = directClassification(2,avEstimate)+1;
                end
                if ismember(neuronNumber, responsiveUnits2.directionSelectiveUnits)
                    directionMap{1,c} = [directionMap{1,c}; directClassification(:,1)'./repeats dp u neuronNumber];
                end
                
            end
%             
%             f1 = figure;
%             set(f1,'Position',[150 100 1250 500]);
%             for cc = 1:length(contrasts)
%                 subplot(2,length(contrasts),cc);
%                 imagesc(squeeze(choiceMap{1,cc}(:,:,end)));
%                 title(['cont = ' num2str(contrasts(cc)) ' vis']);
%                 
%                 subplot(2,length(contrasts),cc+length(contrasts));
%                 imagesc(squeeze(choiceMap{2,cc}(:,:,end)));
%                 title(['cont = ' num2str(contrasts(cc)) ' audvis']);
%             end
%             suptitle([num2str(dp) ' unit ' num2str(neuronNumber) ' (u=' num2str(u) ')']);
%             saveas(f1,fullfile(saveDir,[num2str(dp) ' unit ' num2str(neuronNumber) ' (u=' num2str(u) ') heat map.fig']));
%             saveas(f1,fullfile(saveDir,[num2str(dp) ' unit ' num2str(neuronNumber) ' (u=' num2str(u) ') heat map.jpg']));
%             close(f1);
        end
    end
end

save(fullfile(saveDir,'MLE_IndividualNeuron_Data.mat'),'choiceMap','accuracy','orientationMap','directionMap');

%% Overall accuracy
% meanMap = squeeze(mean(choiceMap{1,5},3));
% figure;imagesc(meanMap);
statsMat = [];
for c = 1:5
    clear tempV tempAV tempDiff
    for s=1:size(choiceMap{1,c},3)
        tempV(s) = mean(diag(choiceMap{1,c}(:,:,s)));
        tempAV(s) = mean(diag(choiceMap{2,c}(:,:,s)));
        tempDiff(s) = tempAV(s) - tempV(s);
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

[p,~,~] = anova2(statsMat,118,'off')

figure;
scatter(accuracy{3}(:,1),accuracy{3}(:,2)-accuracy{3}(:,1));
xlabel('MLE accuracy (visual)');
ylabel('\Delta MLE accuracy (AV - V)');
title('Overall accuracy');

%% Orientation accuracy
for c = 1:5
    clear tempV tempAV tempDiff
    
    ooV = orientationMap{c}(:,1);
    ooAV = orientationMap{c}(:,2);
    ooDiff = ooAV - ooV;
    vAcc(c) = mean(ooV);
    vAccStd(c) = std(ooV)/sqrt(length(ooV));
    avAcc(c) = mean(ooAV);
    avAccStd(c) = std(ooAV)/sqrt(length(ooAV));
    diffAcc(c) = mean(ooDiff);
    diffAccStd(c) = std(ooDiff)/sqrt(length(ooDiff));
    
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
title('Orientation');

figure;hold on;
hA = area(contrasts,[diffAcc-diffAccStd; 2*diffAccStd]')
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(contrasts,diffAcc,'Color',[0 0 0]);
title('Orientation');

figure;
scatter(orientationMap{3}(:,1),orientationMap{3}(:,2)-orientationMap{3}(:,1));
xlabel('MLE accuracy (visual)');
ylabel('\Delta MLE accuracy (AV - V)');
title('Orientation');

%% Direction accuracy
for c = 1:5
    clear tempV tempAV tempDiff
    
    ddV = directionMap{c}(:,1);
    ddAV = directionMap{c}(:,2);
    ddDiff = ddAV - ddV;
    vAcc(c) = mean(ddV);
    vAccStd(c) = std(ddV)/sqrt(length(ddV));
    avAcc(c) = mean(ddAV);
    avAccStd(c) = std(ddAV)/sqrt(length(ddAV));
    diffAcc(c) = mean(ddDiff);
    diffAccStd(c) = std(ddDiff)/sqrt(length(ddDiff));
    
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
title('Direction');

figure;hold on;
hA = area(contrasts,[diffAcc-diffAccStd; 2*diffAccStd]')
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(contrasts,diffAcc,'Color',[0 0 0]);
title('Direction');



