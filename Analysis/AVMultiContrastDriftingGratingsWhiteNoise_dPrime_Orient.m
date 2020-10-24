
clear all;
% close all;

saveDir = fullfile('E:\Electrophysiology','EP004','MetaData - AVMultiContrastDriftingGratingsWhiteNoise_Orientations','D prime - individual neuron');
if ~exist(saveDir)
    mkdir(saveDir);
end
% 
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

blueMap = [linspace(0,1,128)' linspace(0,1,128)' ones(128,1)];
redMap = [ones(128,1) linspace(1,0,128)' linspace(1,0,128)'];
map = [blueMap; redMap];

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
    
    
    if ~exist('dprime','var')
        dprime = cell(length(contrasts),2);
        highPerformance = cell(length(contrasts),2);
        lowPerformance = cell(length(contrasts),2);
        highVsLow = cell(1,2);
    end
    
    
    for u = 1:length(unitData)
        neuronNumber = unitData(u).neuronNumber;
        if (soundResponsive == ismember(neuronNumber,responsiveUnits.soundResponsiveUnits))...       
                && (orientationSelective == ismember(neuronNumber,responsiveUnits.orientationSelectiveUnits))...
                && (lightResponsive == ismember(neuronNumber,responsiveUnits.lightResponsiveUnits))...
                && unitData(u).type~=0

             [dp length(unitData) u]
%             data = mean(unitData(u).frTrainTrials(:,:,quantFirstBin:quantLastBin),3);
            data = unitData(u).trialResponse;
            
            dPrimeMat = zeros(length(contrasts),length(orientations),length(orientations),2);
            alignedMat = cell(1,length(contrasts));
            tempHighVsLow = zeros(2,length(contrasts));
            for c = 1:length(contrasts)
                
                ind = (c-1)*length(orientations)+1;
                indOffset = length(orientations)*length(contrasts);
                
                [maxValV maxIndV] = max(unitData(u).meanResponse(ind:ind+length(orientations)-1,4));
                [minValV minIndV] = min(unitData(u).meanResponse(ind:ind+length(orientations)-1,4));
                [maxValAV maxIndAV] = max(unitData(u).meanResponse(ind+indOffset:ind+length(orientations)-1+indOffset,4));
                [minValAV minIndAV] = min(unitData(u).meanResponse(ind+indOffset:ind+length(orientations)-1+indOffset,4));
                targetIndex = floor(length(orientations)/2) + 1;
                
%                 dirOfInterest = unique([maxIndV minIndV maxIndAV minIndAV]);
                for iDir = 1:length(orientations)
%                 for ii = 1:length(dirOfInterest)
%                     iDir = dirOfInterest(ii);
                    iDataV = data(ind+iDir-1,:);
                    iDataAV = data(ind+iDir-1+indOffset,:);
                    
                    for jDir = 1:length(orientations)
                        jDataV = data(ind+jDir-1,:);
                        jDataAV = data(ind+jDir-1+indOffset,:);
                        
                        iMeanV = mean(iDataV);
                        jMeanV = mean(jDataV);
                        iSTDV = std(iDataV);
                        jSTDV = std(jDataV);
                        dPrimeV = (iMeanV - jMeanV) / sqrt(0.5*(iSTDV^2 + jSTDV^2));
                        
                        iMeanAV = mean(iDataAV);
                        jMeanAV = mean(jDataAV);
                        iSTDAV = std(iDataAV);
                        jSTDAV = std(jDataAV);
                        dPrimeAV = (iMeanAV - jMeanAV) / sqrt(0.5*(iSTDAV^2 + jSTDAV^2));
                        
                        dPrimeMat(c,iDir,jDir,:) = [dPrimeV dPrimeAV];
                    end
                end
                tempHighVsLow(1,c) = dPrimeMat(c,maxIndV,minIndV,1);
                tempHighVsLow(2,c) = dPrimeMat(c,maxIndAV,minIndAV,2);
                
                
                    alignedV = squeeze(dPrimeMat(c,:,:,1));
                    alignedV= circshift(alignedV,[targetIndex-maxIndV targetIndex-maxIndV]);
                    alignedV = [alignedV alignedV(:,1); alignedV(1,:) alignedV(1,1)];
                    
                    alignedAV = squeeze(dPrimeMat(c,:,:,2));
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
                    
                    dprime{c,1} = cat(3,dprime{c,1},alignedMat{c}(:,:,1));
                    dprime{c,2} = cat(3,dprime{c,2},alignedMat{c}(:,:,2));
                
                singleOrient = circshift(squeeze(dPrimeMat(c,maxIndV,:,1))',targetIndex-maxIndV);
                singleOrient = [singleOrient singleOrient(1)];
                highPerformance{c,1} = [highPerformance{c,1}; singleOrient];
                
                singleOrient = circshift(squeeze(dPrimeMat(c,maxIndAV,:,2))',targetIndex-maxIndAV);
                singleOrient = [singleOrient singleOrient(1)];
                highPerformance{c,2} = [highPerformance{c,2}; singleOrient];
                
                singleOrient = circshift(squeeze(dPrimeMat(c,minIndV,:,1))',targetIndex-minIndV);
                singleOrient = [singleOrient singleOrient(1)];
                lowPerformance{c,1} = [lowPerformance{c,1}; singleOrient];
                
                singleOrient = circshift(squeeze(dPrimeMat(c,minIndAV,:,2))',targetIndex-minIndAV);
                singleOrient = [singleOrient singleOrient(1)];
                lowPerformance{c,2} = [lowPerformance{c,2}; singleOrient];
            end
            
            highVsLow{1} = [highVsLow{1}; tempHighVsLow(1,:)];
            highVsLow{2} = [highVsLow{2}; tempHighVsLow(2,:)];

            
%             performance{1} = [performance{1}; squeeze(mean(mean(accuracyMat(:,:,:,1),2),3))'];
%             performance{2} = [performance{2}; squeeze(mean(mean(accuracyMat(:,:,:,2),2),3))'];
            
            f1 = figure;
            set(f1,'Position',[300 150 1000 600]);
            for c = 1:length(contrasts)
                subplot(3,length(contrasts),c);
                imagesc(alignedMat{c}(:,:,1));
%                 imagesc(squeeze(accuracyMat(c,:,:,1)),[0.5 1]);
                colormap(gca,'copper');
                xticks(1:3:7);
                yticks(1:3:7);
                xticklabels(-90:90:90);
                yticklabels(-90:90:90);
                title(['C = ' num2str(contrasts(c)) ', vis']);
%                 imshow(squeeze(heatMat(c,:,:,1,:)));
                
                subplot(3,length(contrasts),c+length(contrasts));
                imagesc(alignedMat{c}(:,:,2));
%                 imagesc(squeeze(accuracyMat(c,:,:,2)),[0.5 1]);
                colormap(gca,'copper');
                xticks(1:3:7);
                yticks(1:3:7);
                xticklabels(-90:90:90);
                yticklabels(-90:90:90);
                title(['C = ' num2str(contrasts(c)) ', audiovis']);
%                 imshow(squeeze(heatMat(c,:,:,2,:)));
                
                aa = alignedMat{c}(:,:,1);bb = alignedMat{c}(:,:,2);
                maxVal = max(abs(bb(:)-aa(:)));
                subplot(3,length(contrasts),c+2*length(contrasts));
                imagesc(alignedMat{c}(:,:,2) - alignedMat{c}(:,:,1),[-maxVal maxVal]);
%                 imagesc(squeeze(accuracyMat(c,:,:,2) - accuracyMat(c,:,:,1)),[-0.5 0.5]);
                colormap(gca,map);
                xticks(1:3:7);
                yticks(1:3:7);
                xticklabels(-90:90:90);
                yticklabels(-90:90:90);
                title(['C = ' num2str(contrasts(c)) ', difference']);
            end
            suptitle([num2str(dp) ' unit ' num2str(neuronNumber) ' (u=' num2str(u) ')']);
            saveas(f1,fullfile(saveDir,[num2str(dp) ' unit ' num2str(neuronNumber) ' (u=' num2str(u) ') heat map.fig']));
            saveas(f1,fullfile(saveDir,[num2str(dp) ' unit ' num2str(neuronNumber) ' (u=' num2str(u) ') heat map.jpg']));
            close(f1);
            
            f2 = figure;
            set(f2,'Position',[200 250 1100 400]);
            for c = 1:length(contrasts)
                v = alignedMat{c}(:,:,1);
                av = alignedMat{c}(:,:,2);
                
                maxVal = max(abs(av(:)-v(:)));
                subplot(2,length(contrasts),c);
                imagesc(av-v,[-maxVal maxVal]);
                colormap(gca,map);
                xticks(1:3:7);
                yticks(1:3:7);
                xticklabels(-90:90:90);
                yticklabels(-90:90:90);
                title(['C = ' num2str(contrasts(c)) ', difference']);
                
                subplot(2,length(contrasts),c+length(contrasts));
                hold on;
                scatter(v(:),av(:));
                theMax = max([max(v(:)) max(av(:))]);
                line([-theMax theMax],[-theMax theMax],'Color',[0 0 0]);
                xlabel('D prime (visual)');
                ylabel('D prime (audiovisual)');
                [h p] = ttest2(v(:),av(:));
                deltAccuracy = mean(av(:)-v(:));
                title(['\Delta = ' num2str(deltAccuracy) 'p = ' num2str(p)]);
                
            end
            suptitle([num2str(dp) ' unit ' num2str(neuronNumber) ' (u=' num2str(u) ')']);
            saveas(f2,fullfile(saveDir,[num2str(dp) ' unit ' num2str(neuronNumber) ' (u=' num2str(u) ') diff map.fig']));
            saveas(f2,fullfile(saveDir,[num2str(dp) ' unit ' num2str(neuronNumber) ' (u=' num2str(u) ') diff map.jpg']));
            close(f2);
            
            f3 = figure;
            set(f3,'Position',[200 250 1100 400]);
            for c = 1:length(contrasts)
                subplot(2,length(contrasts),c);hold on;
                plot(-90:30:90,highPerformance{c,1}(end,:),'Color',[0 0 0]);
                plot(-90:30:90,highPerformance{c,2}(end,:),'Color',[0 0 1]);
                xticks(-90:90:90);
                xlabel('\Delta orientation');
                ylabel('D prime');
                
                subplot(2,length(contrasts),c+length(contrasts));hold on;
                plot(-90:30:90,lowPerformance{c,1}(end,:),'Color',[0 0 0]);
                plot(-90:30:90,lowPerformance{c,2}(end,:),'Color',[0 0 1]);
                xticks(-90:90:90);
                xlabel('\Delta orientation');
                ylabel('D prime');
            end
            suptitle([num2str(dp) ' unit ' num2str(neuronNumber) ' (u=' num2str(u) ')']);
            saveas(f3,fullfile(saveDir,[num2str(dp) ' unit ' num2str(neuronNumber) ' (u=' num2str(u) ') preferred and anti dir.fig']));
            saveas(f3,fullfile(saveDir,[num2str(dp) ' unit ' num2str(neuronNumber) ' (u=' num2str(u) ') preferred and anti dir.jpg']));
            close(f3);
        end
    end
end 

f6 = figure;hold on;
vv = nanmean(highVsLow{1,1});av = nanmean(highVsLow{1,2});
    vvSTD = nanstd(highVsLow{1,1})./sqrt(size(highVsLow{1,1},1));
    avSTD = nanstd(highVsLow{1,2})./sqrt(size(highVsLow{1,2},1));
    
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
plot(contrasts,mean(highVsLow{1}),'Color',[0 0 0]);
plot(contrasts,mean(highVsLow{2}),'Color',[0 0 1]);


diffData = zeros(length(contrasts),size(highPerformance{c,1},1));
f4 = figure;
set(f4,'Position',[300 150 1000 600]);
for c = 1:length(contrasts)
    subplot(2,length(contrasts),c);hold on;
    vv = nanmean(highPerformance{c,1});av = nanmean(highPerformance{c,2});
    vvSTD = nanstd(highPerformance{c,1})./sqrt(size(highPerformance{c,1},1));
    avSTD = nanstd(highPerformance{c,2})./sqrt(size(highPerformance{c,2},1));
    diffData(c,:) = mean(highPerformance{c,2}(:,[1:3, 4:7]) - highPerformance{c,1}(:,[1:3, 4:7]),2)';
    
    hA = area(-90:30:90,[vv-vvSTD; 2*vvSTD]');
    hA(1).FaceAlpha = 0;
    hA(1).EdgeColor = [1 1 1];
    hA(2).FaceColor = [0 0 0];
    hA(2).FaceAlpha = 0.5;
    hA(2).EdgeColor = [1 1 1];
    hA = area(-90:30:90,[av-avSTD; 2*avSTD]');
    hA(1).FaceAlpha = 0;
    hA(1).EdgeColor = [1 1 1];
    hA(2).FaceColor = [0 0 1];
    hA(2).FaceAlpha = 0.5;
    hA(2).EdgeColor = [1 1 1];
    plot(-90:30:90,nanmean(highPerformance{c,1}),'Color',[0 0 0]);
    plot(-90:30:90,nanmean(highPerformance{c,2}),'Color',[0 0 1]);
%     ylim([0 1.2]);
    xticks(-90:90:90);
    xlabel('\Delta orientation');
    ylabel('D prime');
    title(['Cont = ' num2str(contrasts(c))]);
    
    subplot(2,length(contrasts),c+length(contrasts));hold on;
    
    dd = nanmean(highPerformance{c,2} - highPerformance{c,1});
    ddSTD = nanstd(highPerformance{c,2} - highPerformance{c,1})./sqrt(size(highPerformance{c,2} - highPerformance{c,1},1));
    
    hA = area(-90:30:90,[dd-ddSTD; 2*ddSTD]');
    hA(1).FaceAlpha = 0;
    hA(1).EdgeColor = [1 1 1];
    hA(2).FaceColor = [0 0 0];
    hA(2).FaceAlpha = 0.5;
    hA(2).EdgeColor = [1 1 1];
    plot(-90:30:90,dd,'Color',[0 0 0]);
%     ylim([-0.3 0.3]);
    xticks(-90:90:90);
    xlabel('\Delta orientation');
    ylabel('\Delta D prime');
    
    
%     subplot(2,length(contrasts),c+length(contrasts));hold on;
%     vv = nanmean(lowPerformance{c,1});av = nanmean(lowPerformance{c,2});
%     vvSTD = nanstd(lowPerformance{c,1})./sqrt(size(lowPerformance{c,1},1));
%     avSTD = nanstd(lowPerformance{c,2})./sqrt(size(lowPerformance{c,2},1));
%     
%     hA = area(-180:30:180,[vv-vvSTD; 2*vvSTD]');
%     hA(1).FaceAlpha = 0;
%     hA(1).EdgeColor = [1 1 1];
%     hA(2).FaceColor = [0 0 0];
%     hA(2).FaceAlpha = 0.5;
%     hA(2).EdgeColor = [1 1 1];
%     hA = area(-180:30:180,[av-avSTD; 2*avSTD]');
%     hA(1).FaceAlpha = 0;
%     hA(1).EdgeColor = [1 1 1];
%     hA(2).FaceColor = [0 0 1];
%     hA(2).FaceAlpha = 0.5;
%     hA(2).EdgeColor = [1 1 1];
%     plot(-180:30:180,nanmean(lowPerformance{c,1}),'Color',[0 0 0]);
%     plot(-180:30:180,nanmean(lowPerformance{c,2}),'Color',[0 0 1]);
% %     ylim([0.5 0.8]);
%     xticks(-180:90:180);
end
saveas(f4,fullfile(saveDir,['Preferred orientations traces.fig']));

figure;hold on;
    ddData = mean(diffData,2)';
    ddDataSTD = std(diffData,[],2)'/sqrt(size(diffData,2));
    
    hA = area(contrasts,[ddData-ddDataSTD; 2*ddDataSTD]');
    hA(1).FaceAlpha = 0;
    hA(1).EdgeColor = [1 1 1];
    hA(2).FaceColor = [0 0 0];
    hA(2).FaceAlpha = 0.5;
    hA(2).EdgeColor = [1 1 1];
    plot(contrasts,ddData,'Color',[0 0 0]);
%     ylim([-0.3 0.3]);
    xticks(contrasts);
    xlabel('Contrast');
    ylabel('\Delta D prime');





maxVal = 0;
for c = 1:length(contrasts)
    dif = abs(nanmean(dprime{c,2},3) - nanmean(dprime{c,1},3));
    maxVal = max([maxVal max(dif(:))]);
end
f5 = figure;
set(f5,'Position',[300 150 1000 600]);
for c = 1:length(contrasts)
    vv = nanmean(dprime{c,1},3);
    av = nanmean(dprime{c,2},3);
    
    subplot(3,length(contrasts),c);
    imagesc(vv);
    colormap(gca,'parula');
    xticks(1:3:7);
    yticks(1:3:7);
    xticklabels(-90:90:90);
    yticklabels(-90:90:90);
    title(['C = ' num2str(contrasts(c)) ', vis']);
    
    subplot(3,length(contrasts),c+length(contrasts));
    imagesc(av);
    colormap(gca,'parula');
    xticks(1:3:7);
    yticks(1:3:7);
    xticklabels(-90:90:90);
    yticklabels(-90:90:90);
    title(['C = ' num2str(contrasts(c)) ', audiovis']);
    
    subplot(3,length(contrasts),c+2*length(contrasts));
    imagesc(av-vv,[-maxVal maxVal]);
    colormap(gca,map);
    xticks(1:3:7);
    yticks(1:3:7);
    xticklabels(-90:90:90);
    yticklabels(-90:90:90);
    title(['C = ' num2str(contrasts(c)) ', difference']);
end
% saveas(f5,fullfile(saveDir,['Average heat map.fig']));

f7 = figure;
set(f7,'Position',[200 150 1200 600]);
for c = 1:length(contrasts)
    vv = nanmean(dprime{c,1},3);
    av = nanmean(dprime{c,2},3);
    
    subplot(2,length(contrasts),c);
    imagesc(av-vv,[-maxVal maxVal]);
    colormap(gca,map);
    xticks(1:3:7);
    yticks(1:3:7);
    xticklabels(-90:90:90);
    yticklabels(-90:90:90);
    title(['C = ' num2str(contrasts(c)) ', difference']);
    
    subplot(2,length(contrasts),c+length(contrasts));hold on;
    scatter(vv(:),av(:));
%     line([0 max([max(vv(:)) max(av(:))])],[0 max([max(vv(:)) max(av(:))])],'Color',[0 0 0]);
    line([-1 1],[-1 1],'Color',[0 0 0]);
    xlabel('D prime (visual)');
    ylabel('D prime (audiovisual)');
%     xlim([0 1]);ylim([0 1]);
    [hh pp] = ttest(vv(:),av(:));
    deltDprime = mean(av(:)-vv(:));
    title(['\Delta = ' num2str(deltDprime) ', p = ' num2str(pp)]);
    
end
% saveas(f7,fullfile(saveDir,['Heat map scatterplot.fig']));

