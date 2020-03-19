
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


for dp = 1:length(dataPaths)
    dp
    %Load data
    dataFile = fullfile('E:\Electrophysiology\',dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise','AVMultiContrastDriftingGratingsWhiteNoiseData.mat');
    load(dataFile);
    stimPath = dir(fullfile('E:\Electrophysiology\',dataPaths{dp},'StimInfo','*AVmultiContrastDriftingGratingsWhiteNoise_stimInfo*'));
    load(fullfile(stimPath.folder,stimPath.name));
    
    orientations = stimInfo.orientations;
    contrasts = stimInfo.contrasts;
    
    if ~exist('orientationCurves','var')
        orientationCurves = cell(length(contrasts),2);
        alignedCurves = cell(length(contrasts),2);
        bestOrientShift = cell(1,length(contrasts));
    end
    
    aligned = cell(length(contrasts),2);
    unaligned = cell(length(contrasts),2);
    
    for u  = 1:size(unitData,2)
        neuronNumber = unitData(u).neuronNumber;
        if (allNeurons || (lightResponsive && ismember(neuronNumber,responsiveUnits.lightResponsiveUnits))...
                || (orientationSelective && ismember(neuronNumber,responsiveUnits.orientationSelectiveUnits))...
                || (soundResponsive && ismember(neuronNumber,responsiveUnits.soundResponsiveUnits)))...
                && (~onlySingleUnits || unitData(u).type==1) && unitData(u).type~=0
            responses = unitData(u).meanResponse(:,4);
            trialResponses = unitData(u).frTrainTrials;
            
            alignCont = 5;
            highContVis = responses(((alignCont-1)*length(orientations))+1:((alignCont-1)*length(orientations)+length(orientations)));
            highContVisAud = responses(((alignCont-1)*length(orientations))+1+length(contrasts)*length(orientations):((alignCont-1)*length(orientations)+length(orientations)+length(contrasts)*length(orientations)));
            [peak ind] = max(highContVis);
            maxVal = max(responses(:));
            minVal = min(responses(:));
            
            for c = 1:length(contrasts)

                responseVis = responses(((c-1)*length(orientations))+1:((c-1)*length(orientations)+length(orientations)));
                responseVisAud = responses(((c-1)*length(orientations)+1+length(contrasts)*length(orientations)):((c-1)*length(orientations)+length(contrasts)*length(orientations)+length(orientations)));
                
%                 [peak ind] = max(responseVis);
                targetIndex = floor(length(orientations)/2) + 1;
                
                alignedVis = circshift(responseVis,targetIndex - ind);
                alignedVis = [alignedVis; alignedVis(1)];
                alignedVisAud = circshift(responseVisAud,targetIndex - ind);
                alignedVisAud = [alignedVisAud; alignedVisAud(1)];
                
                
                normVis = (alignedVis - minVal) / (maxVal - minVal);
                normVisAud = (alignedVisAud - minVal) / (maxVal - minVal);
                
                aligned{c,1} = [aligned{c,1}; normVis'];
                aligned{c,2} = [aligned{c,2}; normVisAud'];
                
                vNorm = (responseVis - minVal) / (maxVal - minVal);
                vaNorm = (responseVisAud - minVal) / (maxVal - minVal);
                
                unaligned{c,1} = [unaligned{c,1}; vNorm'];
                unaligned{c,2} = [unaligned{c,2}; vaNorm'];
                
                [valV indV] = max(responseVis);
                [valAV indAV] = max(responseVisAud);
                orientDiff = wrapTo180(orientations(indAV) - orientations(indV));
                frDiff = valAV - valV;
                bestOrientShift{c} = [bestOrientShift{c}; indV valV indAV valAV orientDiff frDiff responseVisAud(indV)-valV];
%                 bestOrientShift{c} = [bestOrientShift{c}; indV valV indAV responseVisAud(indV) orientDiff frDiff];
                
            end
        end
    end
    
    
    for c = 1:length(contrasts)
        orientationCurves{c,1} = [orientationCurves{c,1}; mean(unaligned{c,1},1)];
        orientationCurves{c,2} = [orientationCurves{c,2}; mean(unaligned{c,2},1)];
        
        alignedCurves{c,1} = [alignedCurves{c,1}; mean(aligned{c,1},1)];
        alignedCurves{c,2} = [alignedCurves{c,2}; mean(aligned{c,2},1)];
        
    end
end

f1 = figure;
set(f1,'Position',[25 200 1500 400]);
for c = 1:length(contrasts)
    subplot(1,length(contrasts),c);hold on;
    
    semV = std(alignedCurves{c,1})/sqrt(size(alignedCurves{c,1},1));
    semVA = std(alignedCurves{c,2})/sqrt(size(alignedCurves{c,2},2));
    
    hA = area(linspace(-180,180,length(orientations)+1),[mean(alignedCurves{c,1})-semV;2*semV]');
    hA(1).FaceAlpha = 0;
    hA(1).EdgeColor = [1 1 1];
    hA(2).FaceColor = [0 0 0];
    hA(2).FaceAlpha = 0.5;
    hA(2).EdgeColor = [1 1 1];
    
    hAV = area(linspace(-180,180,length(orientations)+1),[mean(alignedCurves{c,2})-semVA;2*semVA]');
    hAV(1).FaceAlpha = 0;
    hAV(1).EdgeColor = [1 1 1];
    hAV(2).FaceColor = [0 0 1];
    hAV(2).FaceAlpha = 0.5;
    hAV(2).EdgeColor = [1 1 1];
    
    
    plot(linspace(-180,180,length(orientations)+1),mean(alignedCurves{c,1}),'Color',[0 0 0],'LineWidth',2);
    plot(linspace(-180,180,length(orientations)+1),mean(alignedCurves{c,2}),'Color',[0 0 1],'LineWidth',2);
    ylabel('Normalized firing rate');
    xlabel('\Delta orientation (degrees)');
    xticks(-180:90:180);
    ylim([0 1])
    title(['Cont = ' num2str(contrasts(c))]);
end
% legend({'Light','Light/Sound'});
suptitle({'Aligned orientation curves, +/- white noise'});
saveas(f1,fullfile(saveDir,'Neuron-wise orientation preference'));

f2 = figure;
set(f2,'Position',[25 200 1500 400]);
for c = 1:length(contrasts)
    subplot(1,length(contrasts),c);
    polarplot(2*pi*orientations./360,mean(orientationCurves{c,1}),'Color',[0 0 0],'LineWidth',2);
    hold on;
    polarplot(2*pi*orientations./360,mean(orientationCurves{c,2}),'Color',[0 0 1],'LineWidth',2);
    rlim([0 1]);
    title(['Cont = ' num2str(contrasts(c))]);
end
suptitle({'Population orientation preference, +/- white noise'});
saveas(f2,fullfile(saveDir,'Population orientation preference'));

% f3 = figure;
% set(f3,'Position',[25 200 1500 400]);
% for c = 1:length(contrasts)
%     subplot(1,length(contrasts),c);hold on;
%     
%     semV = std(orientationCurves{c,1})/sqrt(size(orientationCurves{c,1},1));
%     semVA = std(orientationCurves{c,2})/sqrt(size(orientationCurves{c,2},2));
%     
%     hA = area(orientations,[mean(orientationCurves{c,1})-semV;2*semV]');
%     hA(1).FaceAlpha = 0;
%     hA(1).EdgeColor = [1 1 1];
%     hA(2).FaceColor = [0 0 0];
%     hA(2).FaceAlpha = 0.5;
%     hA(2).EdgeColor = [1 1 1];
%     
%     hAV = area(orientations,[mean(orientationCurves{c,2})-semVA;2*semVA]');
%     hAV(1).FaceAlpha = 0;
%     hAV(1).EdgeColor = [1 1 1];
%     hAV(2).FaceColor = [0 0 1];
%     hAV(2).FaceAlpha = 0.5;
%     hAV(2).EdgeColor = [1 1 1];
%     
%     
%     plot(orientations,mean(orientationCurves{c,1}),'Color',[0 0 0],'LineWidth',2);
%     plot(orientations,mean(orientationCurves{c,2}),'Color',[0 0 1],'LineWidth',2);
%     ylabel('Normalized firing rate');
%     xlabel('Orientation (degrees)');
% %     xticks(-180:90:180);
%     ylim([0 1])
%     title(['Cont = ' num2str(contrasts(c))]);
% end
% % legend({'Light','Light/Sound'});
% suptitle({'Unaligned orientation curves, +/- white noise'});
% % saveas(f3,fullfile(saveDir,'Neuron-wise orientation preference'));

f4 = figure;
set(f4,'Position',[25 100 1500 600]);
for c = 1:length(contrasts)
    subplot(2,length(contrasts),c);
    histogram(bestOrientShift{c}(:,5),(-180-orientations(2)/2):orientations(2):(180+orientations(2)/2));
    xlabel('\Delta Orientation shift');
    ylabel('Number of neurons');
    title(['Contrast = ' num2str(contrasts(c))]);
    ylim([0 60]);
    
    subplot(2,length(contrasts),c+length(contrasts));hold on;
    scatter(bestOrientShift{c}(:,2),bestOrientShift{c}(:,4));
    xx = xlim; yy = ylim; limMax = max([xx(2) yy(2)]);
    line([0 limMax],[0 limMax],'Color',[0 0 0]);
    xlim([0 limMax]);
    ylim([0 limMax]);
    delta = mean(bestOrientShift{c}(:,4) - bestOrientShift{c}(:,2));
    ratio = mean(bestOrientShift{c}(:,4) ./ bestOrientShift{c}(:,2));
    xlabel('FR at best \theta_V');
    ylabel('FR at best \theta_A_V');
    title(['\Delta FR = ' num2str(delta) ' Hz (' num2str((ratio-1)*100) ' %)']);
end
suptitle({'Best orientation and FR differences'});
saveas(f4,fullfile(saveDir,'Orientation and FR shifts'));

f5 = figure;
set(f5,'Position',[25 100 1500 600]);
uniqueBestOrient = -180:orientations(2):180;
orientToFR = zeros(length(contrasts),2,length(uniqueBestOrient)); %mean and sem
orientToBestOrientFR = zeros(length(contrasts),2,length(uniqueBestOrient)); %mean and sem
for c = 1:length(contrasts)
    subplot(2,length(contrasts),c);hold on;
    scatter(bestOrientShift{c}(:,5),bestOrientShift{c}(:,6));
    
% 	orientToFR = zeros(2,length(uniqueBestOrient)); %mean and sem
    for d = 1:size(orientToFR,3)
        dir = uniqueBestOrient(d);
        
        orientToFR(c,1,d) = mean(bestOrientShift{c}(find(bestOrientShift{c}(:,5)==dir),6));
        orientToFR(c,2,d) = std(bestOrientShift{c}(find(bestOrientShift{c}(:,5)==dir),6))/sqrt(length(find(bestOrientShift{c}(:,5)==dir)));
        
        orientToBestOrientFR(c,1,d) = mean(bestOrientShift{c}(find(bestOrientShift{c}(:,5)==dir),7));
        orientToBestOrientFR(c,2,d) = std(bestOrientShift{c}(find(bestOrientShift{c}(:,5)==dir),7))/sqrt(length(find(bestOrientShift{c}(:,5)==dir)));
    
    end
    plot(uniqueBestOrient,squeeze(orientToFR(c,1,:)),'Color',[0 0 0]);
    ylabel('\Delta FR (Hz)');
    xlabel('\Delta best orientation (degrees)');
    title(['Cont = ' num2str(contrasts(c))]);
    
    
    subplot(2,length(contrasts),c+5);hold on;
    scatter(bestOrientShift{c}(:,5),bestOrientShift{c}(:,7));
    
%     orientToBestOrientFR = zeros(2,length(uniqueBestOrient)); %mean and sem
%     for d = 1:size(orientToBestOrientFR,2)
%         dir = uniqueBestOrient(d);
%         orientToBestOrientFR(1,d) = mean(bestOrientShift{c}(find(bestOrientShift{c}(:,5)==dir),7));
%         orientToBestOrientFR(2,d) = std(bestOrientShift{c}(find(bestOrientShift{c}(:,5)==dir),7))/sqrt(length(find(bestOrientShift{c}(:,5)==dir)));
%     end
    plot(uniqueBestOrient,squeeze(orientToBestOrientFR(c,1,:)),'Color',[0 0 0]);
    ylabel('\Delta FR (Hz)');
    xlabel('\Delta best orientation (degrees)');
    title(['Cont = ' num2str(contrasts(c))]);
end
suptitle({'\Delta Peak FR overall (top), and \Delta FR at Best \theta_v_i_s (bottom)'});
saveas(f5,fullfile(saveDir,'FR change at preferred orientations - scatter'));

f6 = figure;
set(f6,'Position',[25 100 1500 600]);
for c = 1:length(contrasts)
    subplot(2,length(contrasts),c);hold on;
    
    hA = area(uniqueBestOrient,squeeze([orientToFR(c,1,:)-orientToFR(c,2,:);2*orientToFR(c,2,:)])');
    hA(1).FaceAlpha = 0;
    hA(1).EdgeColor = [1 1 1];
    hA(2).FaceColor = [0 0 0];
    hA(2).FaceAlpha = 0.5;
    hA(2).EdgeColor = [1 1 1];

    plot(uniqueBestOrient,squeeze(orientToFR(c,1,:)),'Color',[0 0 0],'LineWidth',2);
    ylabel('\Delta FR (Hz)');
    xlabel('\Delta best orientation (degrees)');
    title(['Cont = ' num2str(contrasts(c))]);
    
    
    
    subplot(2,length(contrasts),c+length(contrasts));hold on;
    
    hA = area(uniqueBestOrient,squeeze([orientToBestOrientFR(c,1,:)-orientToBestOrientFR(c,2,:);2*orientToBestOrientFR(c,2,:)])');
    hA(1).FaceAlpha = 0;
    hA(1).EdgeColor = [1 1 1];
    hA(2).FaceColor = [0 0 0];
    hA(2).FaceAlpha = 0.5;
    hA(2).EdgeColor = [1 1 1];

    plot(uniqueBestOrient,squeeze(orientToBestOrientFR(c,1,:)),'Color',[0 0 0],'LineWidth',2);
    ylabel('\Delta FR (Hz)');
    xlabel('\Delta best orientation (degrees)');
    title(['Cont = ' num2str(contrasts(c))]);
end
suptitle({'\Delta Peak FR overall (top), and \Delta FR at Best \theta_v_i_s (bottom)'});
saveas(f6,fullfile(saveDir,'FR change at preferred orientations'));

% figure;
% for c = 1:5
%     subplot(1,5,c);
%     histogram(bestOrientShift{c}(:,1));
% end
