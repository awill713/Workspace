
clear all;
% close all;

saveDir = fullfile('E:\Electrophysiology','EP004','MetaData - AVMultiContrastDriftingGratingsWhiteNoise','SVM - individual neuron');
if ~exist(saveDir)
    mkdir(saveDir);
end

% dataPaths{1} = fullfile('EP004','AW117','20200221-1');
% dataPaths{2} = fullfile('EP004','AW117','20200221-2');
% dataPaths{3} = fullfile('EP004','AW118','20200221-1');
% dataPaths{4} = fullfile('EP004','AW118','20200221-2');
% dataPaths{5} = fullfile('EP004','AW121','20200226-1');
% dataPaths{6} = fullfile('EP004','AW121','20200226-2');
dataPaths{1} = fullfile('EP004','AW124','20200303-1');
dataPaths{2} = fullfile('EP004','AW124','20200303-2');

%which neurons to include
onlySingleUnits = 0;
allNeurons = 0;
lightResponsive = 0;
orientationSelective = 1;
soundResponsive = 1;


trainingTrialsPercent = 0.8;
randomizations = 50;

colorMap = 'parula';
map = colormap(colorMap);close;
% rasterColorMap = 'cool';
% map = colormap(rasterColorMap);close;
% optoColors = map(round(linspace(1,length(map),uniqueEvents-2)),:);
% rasterColors = [0 0 0; 0 0 1; optoColors];


for dp = 1:length(dataPaths)
    dp
    %     [randomize dp]
    %Load data
    dataFile = fullfile('E:\Electrophysiology\',dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise','AVMultiContrastDriftingGratingsWhiteNoiseData.mat');
    load(dataFile);
    stimPath = dir(fullfile('E:\Electrophysiology\',dataPaths{dp},'StimInfo','*AVmultiContrastDriftingGratingsWhiteNoise_stimInfo*'));
    load(fullfile(stimPath.folder,stimPath.name));
    
    orientations = stimInfo.orientations;
    contrasts = stimInfo.contrasts;
    repeats = stimInfo.repeats;
    
    
    if ~exist('performance','var')
        performance = cell(length(contrasts),2);
        highPerformance = cell(length(contrasts),2);
        lowPerformance = cell(length(contrasts),2);
    end
    
    
    for u = 1:length(unitData)
        neuronNumber = unitData(u).neuronNumber;
        if (soundResponsive == ismember(neuronNumber,responsiveUnits.soundResponsiveUnits))...
                && (orientationSelective == ismember(neuronNumber,responsiveUnits.orientationSelectiveUnits))...
                && unitData(u).type~=0
               
%             data = mean(unitData(u).frTrainTrials(:,:,quantFirstBin:quantLastBin),3);
            data = unitData(u).trialResponse;
            
            accuracyMat = zeros(length(contrasts),length(orientations),length(orientations),2);
            alignedMat = cell(1,length(contrasts));
            for c = 1:length(contrasts)
                ind = (c-1)*length(orientations)+1;
                indOffset = length(orientations)*length(contrasts);
                
                [maxValV maxIndV] = max(unitData(u).meanResponse(ind:ind+length(orientations)-1,4));
                [minValV minIndV] = min(unitData(u).meanResponse(ind:ind+length(orientations)-1,4));
                [maxValAV maxIndAV] = max(unitData(u).meanResponse(ind+indOffset:ind+length(orientations)-1+indOffset,4));
                [minValAV minIndAV] = min(unitData(u).meanResponse(ind+indOffset:ind+length(orientations)-1+indOffset,4));
                targetIndex = floor(length(orientations)/2) + 1;
                
                dirOfInterest = unique([maxIndV minIndV maxIndAV minIndAV]);
%                 for iDir = 1:length(orientations)
                for ii = 1:length(dirOfInterest)
                    iDir = dirOfInterest(ii);
                    iDataV = data(ind+iDir-1,:);
                    iDataAV = data(ind+iDir-1+indOffset,:);
                    
                    for jDir = 1:length(orientations)
                        jDataV = data(ind+jDir-1,:);
                        jDataAV = data(ind+jDir-1+indOffset,:);
                        
                        accuracy = zeros(2,repeats);
%                         for r = 1:randomizations
                        for r = 1:repeats
                            [dp length(unitData) u c iDir jDir r]
                            
                            iTestTrials = r;
                            iTrainingTrials = setdiff(1:repeats,iTestTrials);
                            jTrainingTrials = 1:10;
                            jTestTrials = [];
                            
%                             trainingCount = floor(trainingTrialsPercent*repeats);
%                             testCount = repeats - trainingCount;
%                             trainingCount = 8;
%                             testCount = 5;
%                             
%                             iTrainingTrials = randperm(repeats,trainingCount);
%                             jTrainingTrials = randperm(repeats,trainingCount);
%                             iTestTrials = setdiff(1:repeats,iTrainingTrials);
%                             jTestTrials = setdiff(1:repeats,jTrainingTrials);
%                             iTestTrials = randperm(repeats,testCount);
%                             jTestTrials = randperm(repeats,testCount);

                            trainingVec = [iDataV(iTrainingTrials)' iDataAV(iTrainingTrials)';...
                                jDataV(jTrainingTrials)' jDataAV(jTrainingTrials)'];
                            trainingAnswers = [ones(length(iTrainingTrials),1); ones(length(jTrainingTrials),1)*-1];
                            testingVec = [iDataV(iTestTrials)' iDataAV(iTestTrials)';...
                                jDataV(jTestTrials)' jDataAV(jTestTrials)'];
                            testAnswers = [ones(length(iTestTrials),1); ones(length(jTestTrials),1)*-1];
                            
                            
                            svmVis = fitcsvm(trainingVec(:,1),trainingAnswers,'KernelFunction','Linear','Standardize',false,'KernelScale','auto');
                            testPrediction = svmVis.predict(testingVec(:,1));
                            accuracy(1,r) = mean(testPrediction==testAnswers);
                            
                            svmAudVis = fitcsvm(trainingVec(:,2),trainingAnswers,'KernelFunction','Linear','Standardize',false,'KernelScale','auto');
                            testPrediction = svmAudVis.predict(testingVec(:,2));
                            accuracy(2,r) = mean(testPrediction==testAnswers);
                        end
                        
                        accuracyMat(c,iDir,jDir,:) = mean(accuracy,2);
                    end
                end
                alignedV = squeeze(accuracyMat(c,:,:,1));
                    alignedV= circshift(alignedV,[targetIndex-maxIndV targetIndex-maxIndV]);
                    alignedV = [alignedV alignedV(:,1); alignedV(1,:) alignedV(1,1)];
                    
                    alignedAV = squeeze(accuracyMat(c,:,:,2));
                    alignedAV= circshift(alignedAV,[targetIndex-maxIndAV targetIndex-maxIndAV]);
                    alignedAV = [alignedAV alignedAV(:,1); alignedAV(1,:) alignedAV(1,1)];
                    
%                     lowVMat = squeeze(accuracyMat(c,:,:,1));
%                     lowVMat= circshift(lowVMat,[targetIndex-minIndV targetIndex-minIndV]);
%                     lowVMat = [lowVMat lowVMat(:,1); lowVMat(1,:) lowVMat(1,1)];
%                     
%                     lowAVMat = squeeze(accuracyMat(c,:,:,2));
%                     lowAVMat= circshift(lowAVMat,[targetIndex-minIndAV targetIndex-minIndAV]);
%                     lowAVMat = [lowAVMat lowAVMat(:,1); lowAVMat(1,:) lowAVMat(1,1)];
%                     
                    alignedMat{c} = cat(3,alignedV,alignedAV);
                    
                    performance{c,1} = cat(3,performance{c,1},alignedMat{c}(:,:,1));
                    performance{c,2} = cat(3,performance{c,2},alignedMat{c}(:,:,2));
                
                singleOrient = circshift(squeeze(accuracyMat(c,maxIndV,:,1))',targetIndex-maxIndV);
                singleOrient = [singleOrient singleOrient(1)];
                highPerformance{c,1} = [highPerformance{c,1}; singleOrient];
                
                singleOrient = circshift(squeeze(accuracyMat(c,maxIndAV,:,2))',targetIndex-maxIndAV);
                singleOrient = [singleOrient singleOrient(1)];
                highPerformance{c,2} = [highPerformance{c,2}; singleOrient];
                
                singleOrient = circshift(squeeze(accuracyMat(c,minIndV,:,1))',targetIndex-minIndV);
                singleOrient = [singleOrient singleOrient(1)];
                lowPerformance{c,1} = [lowPerformance{c,1}; singleOrient];
                
                singleOrient = circshift(squeeze(accuracyMat(c,minIndAV,:,2))',targetIndex-minIndAV);
                singleOrient = [singleOrient singleOrient(1)];
                lowPerformance{c,2} = [lowPerformance{c,2}; singleOrient];
            end
            
%             performance{1} = [performance{1}; squeeze(mean(mean(accuracyMat(:,:,:,1),2),3))'];
%             performance{2} = [performance{2}; squeeze(mean(mean(accuracyMat(:,:,:,2),2),3))'];
            
            f1 = figure;
            set(f1,'Position',[300 150 1000 600]);
            for c = 1:length(contrasts)
                subplot(3,length(contrasts),c);
                imagesc(alignedMat{c}(:,:,1),[0.5 1]);
%                 imagesc(squeeze(accuracyMat(c,:,:,1)),[0.5 1]);
                colormap(gca,'copper');
                xticks(1:3:13);
                yticks(1:3:13);
                xticklabels(-180:90:180);
                yticklabels(-180:90:180);
                title(['C = ' num2str(c) ', vis']);
%                 imshow(squeeze(heatMat(c,:,:,1,:)));
                
                subplot(3,length(contrasts),c+length(contrasts));
                imagesc(alignedMat{c}(:,:,2),[0.5 1]);
%                 imagesc(squeeze(accuracyMat(c,:,:,2)),[0.5 1]);
                colormap(gca,'copper');
                xticks(1:3:13);
                yticks(1:3:13);
                xticklabels(-180:90:180);
                yticklabels(-180:90:180);
                title(['C = ' num2str(c) ', audiovis']);
%                 imshow(squeeze(heatMat(c,:,:,2,:)));
                
                subplot(3,length(contrasts),c+2*length(contrasts));
                imagesc(alignedMat{c}(:,:,2) - alignedMat{c}(:,:,1),[-0.5 0.5]);
%                 imagesc(squeeze(accuracyMat(c,:,:,2) - accuracyMat(c,:,:,1)),[-0.5 0.5]);
                colormap(gca,'parula');
                xticks(1:3:13);
                yticks(1:3:13);
                xticklabels(-180:90:180);
                yticklabels(-180:90:180);
                title(['C = ' num2str(c) ', difference']);
            end
            suptitle([num2str(dp) ' unit ' num2str(neuronNumber) ' (u=' num2str(u) ')']);
%             saveas(f1,fullfile(saveDir,[num2str(dp) ' unit ' num2str(neuronNumber) ' (u=' num2str(u) ') heat map.fig']));
%             saveas(f1,fullfile(saveDir,[num2str(dp) ' unit ' num2str(neuronNumber) ' (u=' num2str(u) ') heat map.jpg']));
            close(f1);
            
            f2 = figure;
            set(f2,'Position',[200 250 1100 400]);
            for c = 1:length(contrasts)
                v = alignedMat{c}(:,:,1);
                av = alignedMat{c}(:,:,2);
                
                subplot(2,length(contrasts),c);
                imagesc(av-v,[-0.5 0.5]);
                colormap(gca,'parula');
                xticks(1:3:13);
                yticks(1:3:13);
                xticklabels(-180:90:180);
                yticklabels(-180:90:180);
                title(['C = ' num2str(c) ', difference']);
                
                subplot(2,length(contrasts),c+length(contrasts));
                hold on;
                scatter(v(:),av(:));
                line([0.3 1],[0.3 1],'Color',[0 0 0],'LineStyle',':');
                xlabel('Accuracy (visual)');
                ylabel('Accuracy (audiovisual)');
                [h p] = ttest2(v(:),av(:));
                deltAccuracy = mean(av(:)-v(:));
                title(['\Delta = ' num2str(deltAccuracy) 'p = ' num2str(p)]);
                
            end
            suptitle([num2str(dp) ' unit ' num2str(neuronNumber) ' (u=' num2str(u) ')']);
%             saveas(f2,fullfile(saveDir,[num2str(dp) ' unit ' num2str(neuronNumber) ' (u=' num2str(u) ') diff map.fig']));
%             saveas(f2,fullfile(saveDir,[num2str(dp) ' unit ' num2str(neuronNumber) ' (u=' num2str(u) ') diff map.jpg']));
            close(f2);
            
            f3 = figure;
            set(f3,'Position',[200 250 1100 400]);
            for c = 1:length(contrasts)
                subplot(2,length(contrasts),c);hold on;
                plot(-180:30:180,highPerformance{c,1}(end,:),'Color',[0 0 0]);
                plot(-180:30:180,highPerformance{c,2}(end,:),'Color',[0 0 1]);
                
                subplot(2,length(contrasts),c+length(contrasts));hold on;
                plot(-180:30:180,lowPerformance{c,1}(end,:),'Color',[0 0 0]);
                plot(-180:30:180,lowPerformance{c,2}(end,:),'Color',[0 0 1]);
            end
            suptitle([num2str(dp) ' unit ' num2str(neuronNumber) ' (u=' num2str(u) ')']);
%             saveas(f3,fullfile(saveDir,[num2str(dp) ' unit ' num2str(neuronNumber) ' (u=' num2str(u) ') preferred and anti dir.fig']));
%             saveas(f3,fullfile(saveDir,[num2str(dp) ' unit ' num2str(neuronNumber) ' (u=' num2str(u) ') preferred and anti dir.jpg']));
            close(f3);
        end
    end
end  

f4 = figure;
for c = 1:length(contrasts)
    subplot(2,length(contrasts),c);hold on;
    vv = mean(highPerformance{c,1});av = mean(highPerformance{c,2});
    vvSTD = std(highPerformance{c,1})./sqrt(size(vv,2));
    avSTD = std(highPerformance{c,2})./sqrt(size(av,2));
    
    hA = area(-180:30:180,[vv-vvSTD; 2*vvSTD]');
    hA(1).FaceAlpha = 0;
    hA(1).EdgeColor = [1 1 1];
    hA(2).FaceColor = [0 0 0];
    hA(2).FaceAlpha = 0.5;
    hA(2).EdgeColor = [1 1 1];
    hA = area(-180:30:180,[av-avSTD; 2*avSTD]');
    hA(1).FaceAlpha = 0;
    hA(1).EdgeColor = [1 1 1];
    hA(2).FaceColor = [0 0 1];
    hA(2).FaceAlpha = 0.5;
    hA(2).EdgeColor = [1 1 1];
    plot(-180:30:180,mean(highPerformance{c,1}),'Color',[0 0 0]);
    plot(-180:30:180,mean(highPerformance{c,2}),'Color',[0 0 1]);
%     ylim([0.5 0.8]);
    xticks(-180:90:180);
    title(['Cont = ' num2str(c)]);
    
    
    subplot(2,length(contrasts),c+length(contrasts));hold on;
    vv = mean(lowPerformance{c,1});av = mean(lowPerformance{c,2});
    vvSTD = std(lowPerformance{c,1})./sqrt(size(vv,2));
    avSTD = std(lowPerformance{c,2})./sqrt(size(av,2));
    
    hA = area(-180:30:180,[vv-vvSTD; 2*vvSTD]');
    hA(1).FaceAlpha = 0;
    hA(1).EdgeColor = [1 1 1];
    hA(2).FaceColor = [0 0 0];
    hA(2).FaceAlpha = 0.5;
    hA(2).EdgeColor = [1 1 1];
    hA = area(-180:30:180,[av-avSTD; 2*avSTD]');
    hA(1).FaceAlpha = 0;
    hA(1).EdgeColor = [1 1 1];
    hA(2).FaceColor = [0 0 1];
    hA(2).FaceAlpha = 0.5;
    hA(2).EdgeColor = [1 1 1];
    plot(-180:30:180,mean(lowPerformance{c,1}),'Color',[0 0 0]);
    plot(-180:30:180,mean(lowPerformance{c,2}),'Color',[0 0 1]);
%     ylim([0.5 0.8]);
    xticks(-180:90:180);
end
% saveas(f4,fullfile(saveDir,['Preferred and lowest orientations.fig']));
            

f5 = figure;
for c = 1:length(contrasts)
    vv = mean(performance{c,1},3);
    av = mean(performance{c,2},3);
    
    subplot(3,length(contrasts),c);
    imagesc(vv,[0.5 1]);
    colormap(gca,'copper');
    xticks(1:3:13);
    yticks(1:3:13);
    xticklabels(-180:90:180);
    yticklabels(-180:90:180);
    title(['C = ' num2str(c) ', vis']);
    
    subplot(3,length(contrasts),c+length(contrasts));
    imagesc(av,[0.5 1]);
    colormap(gca,'copper');
    xticks(1:3:13);
    yticks(1:3:13);
    xticklabels(-180:90:180);
    yticklabels(-180:90:180);
    title(['C = ' num2str(c) ', audiovis']);
    
    subplot(3,length(contrasts),c+2*length(contrasts));
    imagesc(av-vv,[-0.5 0.5]);
    colormap(gca,'parula');
    xticks(1:3:13);
    yticks(1:3:13);
    xticklabels(-180:90:180);
    yticklabels(-180:90:180);
    title(['C = ' num2str(c) ', difference']);
end
% saveas(f5,fullfile(saveDir,['Average heat map.fig']));
            
    
    
    
    
    
    % figure;hold on;
    % vv = mean(performance{1}); av = mean(performance{2});
    % vvSTD = std(performance{1})./sqrt(size(vv,2));
    % avSTD = std(performance{2})./sqrt(size(av,2));
    % hA = area(contrasts,[vv-vvSTD; 2*vvSTD]')
    %     hA(1).FaceAlpha = 0;
    %     hA(1).EdgeColor = [1 1 1];
    %     hA(2).FaceColor = [0 0 0];
    %     hA(2).FaceAlpha = 0.5;
    %     hA(2).EdgeColor = [1 1 1];
    %     hA = area(contrasts,[av-avSTD; 2*avSTD]')
    %     hA(1).FaceAlpha = 0;
    %     hA(1).EdgeColor = [1 1 1];
    %     hA(2).FaceColor = [0 0 1];
    %     hA(2).FaceAlpha = 0.5;
    %     hA(2).EdgeColor = [1 1 1];
    % plot(contrasts,mean(performance{1}),'Color',[0 0 0]);
    % plot(contrasts,mean(performance{2}),'Color',[0 0 1]);
