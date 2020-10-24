
clear all;
% close all;

saveDir = fullfile('D:\Electrophysiology','EP004','MetaData - AVMultiContrastDriftingGratingsWhiteNoise');
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

soundCount = 0;lightCount = 0;orientationCount=0;
soundLightCount=0;soundOrientationCount=0;lightOrientationCount=0;
allThree=0;
soundUp = 0;soundDown = 0;

orientationShuffles = 1000;

for dp = 1:length(dataPaths)
    dp
    %Load data
    dataFile = fullfile('D:\Electrophysiology\',dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise','AVMultiContrastDriftingGratingsWhiteNoiseData.mat');
    load(dataFile);
    stimPath = dir(fullfile('D:\Electrophysiology\',dataPaths{dp},'StimInfo','*AVmultiContrastDriftingGratingsWhiteNoise_stimInfo*'));
    load(fullfile(stimPath.folder,stimPath.name));
    
    orientations = stimInfo.orientations;
    contrasts = stimInfo.contrasts;
    repeats = stimInfo.repeats;
    
    soundCount = soundCount+length(responsiveUnits.soundResponsiveUnits);
    lightCount = lightCount+length(responsiveUnits.lightResponsiveUnits);
    orientationCount = orientationCount+length(responsiveUnits.orientationSelectiveUnits);
    soundLightCount = soundLightCount+length(responsiveUnits.soundAndLight);
    soundOrientationCount = soundOrientationCount+length(responsiveUnits.soundAndOrientation);
    lightOrientationCount = lightOrientationCount+length(responsiveUnits.lightAndOrientation);
    allThree = allThree+length(intersect(responsiveUnits.soundAndLight,responsiveUnits.orientationSelectiveUnits));
    
    if ~exist('orientationCurves','var')
        orientationCurves = cell(length(contrasts),2);
        alignedCurves = cell(length(contrasts),2);
        bestOrientShift = cell(1,length(contrasts));
        rawTC = cell(length(contrasts),2,2);
        sparsity = zeros(2,length(contrasts),0);
        shuffledOrientShift = zeros(0,orientationShuffles);
        normFRTrain = cell(length(contrasts),2,3);
        timingData = cell(length(contrasts),2,2); %neurons that increase vs neurons that decrease % [latency to peak, response slope, fwhm, latency to response])
        newAligned = cell(length(contrasts),2);
        responseLinearity = cell(length(contrasts),1);
    end
    
    aligned = cell(length(contrasts),2);
    unaligned = cell(length(contrasts),2);
    raw = cell(length(contrasts),2,2);
    
    for u  = 1:size(unitData,2)
        neuronNumber = unitData(u).neuronNumber;
        %         if (allNeurons || (lightResponsive && ismember(neuronNumber,responsiveUnits.lightResponsiveUnits))...
        %                 || (orientationSelective && ismember(neuronNumber,responsiveUnits.orientationSelectiveUnits))...
        %                 || (soundResponsive && ismember(neuronNumber,responsiveUnits.soundResponsiveUnits)))...
        %                 && (~onlySingleUnits || unitData(u).type==1) && unitData(u).type~=0
        if (soundResponsive == ismember(neuronNumber,responsiveUnits.soundResponsiveUnits))...
                && (lightResponsive && ismember(neuronNumber,responsiveUnits.lightResponsiveUnits))...
                && (orientationSelective == ismember(neuronNumber,responsiveUnits.orientationSelectiveUnits))...
                && unitData(u).type~=0
            
                            
            responses = unitData(u).meanResponse(:,4);
            %             trialResponses = unitData(u).frTrainTrials;
            
            alignCont = 5;
            highContVis = responses(((alignCont-1)*length(orientations))+1:((alignCont-1)*length(orientations)+length(orientations)));
            highContVisAud = responses(((alignCont-1)*length(orientations))+1+length(contrasts)*length(orientations):((alignCont-1)*length(orientations)+length(orientations)+length(contrasts)*length(orientations)));
            %             [peak ind] = max(highContVis);
            %             ind = mod((ind+6)-1,12)+1;
            maxValV = max(responses(:));
            minValV = min(responses(:));
            
            if mean(highContVisAud)>mean(highContVis)
                soundUp = soundUp+1;
            else
                soundDown = soundDown+1;
            end
            
            sparsity = cat(3,sparsity,unitData(u).sparsity);
            
            
            theAverages = unitData(u).meanResponse((length(contrasts)-1)*length(orientations)+1:length(contrasts)*length(orientations),4);
            theSTDs = unitData(u).meanResponse((length(contrasts)-1)*length(orientations)+1:length(contrasts)*length(orientations),5)*repeats;
            shuffledList = zeros(1,orientationShuffles);
            [~, thisUnitPeak] = max(theAverages);
            for shuff = 1:orientationShuffles
                thisShuffPeakVal = 0;
                for shuffDir = 1:length(theAverages)
                    dirVal = mean(normrnd(theAverages(shuffDir),theSTDs(shuffDir),1,repeats));
                    if dirVal>thisShuffPeakVal
                        thisShuffPeak = shuffDir;
                        thisShuffPeakVal = dirVal;
                    end
                end
                shuffDiff = wrapTo180(orientations(thisShuffPeak) - orientations(thisUnitPeak));
                if shuffDiff == 180
                    shuffDiff = -180;
                end
                shuffledList(shuff) = shuffDiff;
            end
            shuffledOrientShift = [shuffledOrientShift; shuffledList];
            
            for c = 1:length(contrasts)
                
                responseVis = responses(((c-1)*length(orientations))+1:((c-1)*length(orientations)+length(orientations)));
                responseVisAud = responses(((c-1)*length(orientations)+1+length(contrasts)*length(orientations)):((c-1)*length(orientations)+length(contrasts)*length(orientations)+length(orientations)));
                
                soundResponse = mean(responses(length(orientations)*length(contrasts)+1:length(orientations)*(length(contrasts)+1)));
                linearity = [mean(responseVis) soundResponse mean(responseVisAud) mean(responseVis)+soundResponse dp u neuronNumber];
                responseLinearity{c} = [responseLinearity{c}; linearity];
                
                maxValV = max(responseVis);
                minValV = min(responseVis);
                maxValAV = max(responseVisAud);
                minValAV = min(responseVisAud);
                [peak vInd] = max(responseVis);
                [peak avInd] = max(responseVisAud);
                
                xTime = analysisParams.analysisWindow(1)+analysisParams.frBinWidth : analysisParams.analysisWindow(2);
                baselineBins = [find(xTime==analysisParams.baselineWindow(1)) find(xTime==analysisParams.baselineWindow(2))];
                quantBins = [find(xTime==analysisParams.quantWindow(1)) find(xTime==analysisParams.quantWindow(2))];
                
                vFRTrain = smooth(unitData(u).frTrain(((c-1)*length(orientations))+vInd,:),30);
                avFRTrain = smooth(unitData(u).frTrain(((c-1)*length(orientations))+length(contrasts)*length(orientations)+avInd,:),30);
                aFRTrain = smooth(unitData(u).frTrain(vInd+length(contrasts)*length(orientations),:),30);
                
                vFRBaseline = mean(vFRTrain(baselineBins(1):baselineBins(2)));
                avFRBaseline = mean(vFRTrain(baselineBins(1):baselineBins(2)));
                aFRBaseline = mean(aFRTrain(baselineBins(1):baselineBins(2)));
                normVTrain = (vFRTrain - vFRBaseline)./(max(vFRTrain)-vFRBaseline);
                normAVTrain = (avFRTrain - avFRBaseline)./(max(vFRTrain)-vFRBaseline);
                normATrain = (aFRTrain - aFRBaseline)./(max(vFRTrain)-vFRBaseline);
                
                quantTrainV = vFRTrain(quantBins(1):quantBins(2));
                quantTrainAV = avFRTrain(quantBins(1):quantBins(2));
                
                [peakVal(1) peakLat(1)] = max(quantTrainV);
                [peakVal(2) peakLat(2)] = max(quantTrainAV);
                peakLatBin = peakLat - xTime(1);
                
                respLatV  = find(quantTrainV > (vFRBaseline+std(vFRTrain(baselineBins(1):baselineBins(2)))),1,'first');
                respLatAV = find(quantTrainAV > (avFRBaseline+std(avFRTrain(baselineBins(1):baselineBins(2)))),1,'first');
                
                frSlope(1) = (peakVal(1)-vFRBaseline) / peakLat(1); frSlope(2) = (peakVal(2)-avFRBaseline) / peakLat(2);
                
                %                 smoothedV = smooth(vFRTrain,30); smoothedAV = smooth(avFRTrain,30);
                smoothedV = vFRTrain; smoothedAV = avFRTrain;
                threshV = smoothedV > 0.5*smoothedV(peakLatBin(1));
                threshAV = smoothedAV > 0.5*smoothedAV(peakLatBin(2));
                
                firstV = peakLatBin(1) - find(threshV(1:peakLatBin(1))==0,1,'last');
                lastV = find(threshV(peakLatBin(1):end)==0,1,'first');
                firstAV = peakLatBin(2) - find(threshAV(1:peakLatBin(2))==0,1,'last');
                lastAV = find(threshAV(peakLatBin(2):end)==0,1,'first');
                
                
                
                if mean(vFRTrain(quantBins(1):quantBins(2)))>mean(vFRTrain(baselineBins(1):baselineBins(2)))
                    normFRTrain{c,1,1} = [normFRTrain{c,1,1}; normVTrain'];
                    normFRTrain{c,1,2} = [normFRTrain{c,1,2}; normAVTrain'];
                    normFRTrain{c,1,3} = [normFRTrain{c,1,3}; normATrain'];
                    
                    if (length(firstV) + length(lastV) + length(firstAV) + length(lastAV) + length(respLatV) + length(respLatAV) ==6)
                        fwhm(1) = firstV+lastV;
                        fwhm(2) = firstAV+lastAV;
                        
                        tempTimingV = [peakLat(1) frSlope(1) fwhm(1) respLatV];
                        tempTimingAV =[peakLat(2) frSlope(2) fwhm(2) respLatAV];
                        
                        timingData{c,1,1} = [timingData{c,1,1}; tempTimingV];
                        timingData{c,1,2} = [timingData{c,1,2}; tempTimingAV];
                    end
                    
                    
                else
                    normFRTrain{c,2,1} = [normFRTrain{c,2,1}; normVTrain'];
                    normFRTrain{c,2,2} = [normFRTrain{c,2,2}; normAVTrain'];
                    normFRTrain{c,2,3} = [normFRTrain{c,2,3}; normATrain'];
                    
                    %                     timingData{c,1,1} = [timingData{c,1,1}; tempTimingV];
                    %                     timingData{c,1,2} = [timingData{c,1,2}; tempTimingAV];
                end
                
                if mean(responseVisAud)>mean(responseVis)
                    raw{c,1,1} = [raw{c,1,1}; responseVis'];
                    raw{c,1,2} = [raw{c,1,2}; responseVisAud'];
                else
                    raw{c,2,1} = [raw{c,2,1}; responseVis'];
                    raw{c,2,2} = [raw{c,2,2}; responseVisAud'];
                end
                
                %                 [peak ind] = max(responseVis);
                targetIndex = floor(length(orientations)/2) + 1;
                
                alignedVis = circshift(responseVis,targetIndex - vInd);
                alignedVis = [alignedVis; alignedVis(1)];
                alignedVisAud = circshift(responseVisAud,targetIndex - avInd);
                alignedVisAud = [alignedVisAud; alignedVisAud(1)];
                
                newAligned{c,1} = [newAligned{c,1}; (alignedVis./alignedVis)'];
                newAligned{c,2} = [newAligned{c,2}; (alignedVisAud./alignedVis)'];
                
                
                normVis = (alignedVis - minValV) / (maxValV - minValV);
                normVisAud = (alignedVisAud - minValV) / (maxValV - minValV);
                
                aligned{c,1} = [aligned{c,1}; normVis'];
                aligned{c,2} = [aligned{c,2}; normVisAud'];
                
                
                vNorm = (responseVis - minValV) / (maxValV - minValV);
                vaNorm = (responseVisAud - minValV) / (maxValV - minValV);
                
                unaligned{c,1} = [unaligned{c,1}; vNorm'];
                unaligned{c,2} = [unaligned{c,2}; vaNorm'];
                
                [valV indV] = max(responseVis);
                [valAV indAV] = max(responseVisAud);
                orientDiff = wrapTo180(orientations(indAV) - orientations(indV));
                if orientDiff == 180
                    orientDiff = -180;
                end
                frDiff = valAV - valV;
                %                 if mean(responseVisAud)<mean(responseVis)
                bestOrientShift{c} = [bestOrientShift{c}; indV valV indAV valAV orientDiff frDiff responseVisAud(indV)-valV];
                %                 end
                %                                                 bestOrientShift{c} = [bestOrientShift{c}; indV valV indAV responseVisAud(indV) orientDiff frDiff responseVisAud(indV)-valV];
                %                 bestOrientShift{c} = [bestOrientShift{c}; indV valV indAV responseVisAud(indV) orientDiff frDiff];
                
            end
        end
    end
    
    
    for c = 1:length(contrasts)
        orientationCurves{c,1} = [orientationCurves{c,1}; mean(unaligned{c,1},1)];
        orientationCurves{c,2} = [orientationCurves{c,2}; mean(unaligned{c,2},1)];
        
        alignedCurves{c,1} = [alignedCurves{c,1}; mean(aligned{c,1},1)];
        alignedCurves{c,2} = [alignedCurves{c,2}; mean(aligned{c,2},1)];
        
        rawTC{c,1,1} = [rawTC{c,1,1}; mean(raw{c,1,1},1)];
        rawTC{c,1,2} = [rawTC{c,1,2}; mean(raw{c,1,2},1)];
        
        rawTC{c,2,1} = [rawTC{c,2,1}; mean(raw{c,2,1},1)];
        rawTC{c,2,2} = [rawTC{c,2,2}; mean(raw{c,2,2},1)];
        
    end
end

soundUpPercent = round(soundUp/(soundUp+soundDown)*100);
soundDownPercent = 100-soundUpPercent;
f7 = figure;
set(f7,'Position',[25 200 1500 400]);
for c = 1:length(contrasts)
    subplot(1,length(contrasts),c);
    polarplot(2*pi*orientations./360,mean(rawTC{c,1,1}),'Color',[0 0 0],'LineWidth',2);
    hold on;
    polarplot(2*pi*orientations./360,mean(rawTC{c,1,2}),'Color',[0 0 1],'LineWidth',2);
    rlim([0 25]);
    title(['Cont = ' num2str(contrasts(c))]);
    suptitle(['Sound increases light-evoked responses (' num2str(soundUpPercent) '%)']);
end
f8 = figure;
set(f8,'Position',[25 200 1500 400]);
for c = 1:length(contrasts)
    subplot(1,length(contrasts),c);
    polarplot(2*pi*orientations./360,mean(rawTC{c,2,1}),'Color',[0 0 0],'LineWidth',2);
    hold on;
    polarplot(2*pi*orientations./360,mean(rawTC{c,2,2}),'Color',[0 0 1],'LineWidth',2);
    rlim([0 15]);
    title(['Cont = ' num2str(contrasts(c))]);
    suptitle(['Sound decreases light-evoked responses (' num2str(soundDownPercent) '%)']);
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
%     ylim([0 1])
    title(['Cont = ' num2str(contrasts(c))]);
end
% legend({'Light','Light/Sound'});
suptitle({'Aligned orientation curves, +/- white noise'});
% saveas(f1,fullfile(saveDir,'Neuron-wise orientation preference'));

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
% saveas(f2,fullfile(saveDir,'Population orientation preference'));

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

shuffledOrientShift;
shuffPercent = zeros(orientationShuffles,length(orientations)+1);
shuffVector = zeros(1,orientationShuffles);
for sh = 1:orientationShuffles
    percent = histcounts(shuffledOrientShift(:,sh),(-180-orientations(2)/2):orientations(2):(180+orientations(2)/2),'Normalization','probability');
    percent(end) = percent(1);
    shuffPercent(sh,:) = percent;
    meanVector = mean(abs(shuffledOrientShift(:,sh)));
    shuffVector(sh) = meanVector;
end
percentMean = mean(shuffPercent);
percentStd = std(shuffPercent);
actualPercent = histcounts(bestOrientShift{length(contrasts)}(:,5),(-180-orientations(2)/2):orientations(2):(180+orientations(2)/2),'Normalization','probability');
actualPercent(end) = actualPercent(1);
actualVector = mean(abs(bestOrientShift{length(contrasts)}(:,5)));
figure;hold on;
hA = area(-180:30:180,[percentMean-percentStd; 2*percentStd]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(-180:30:180,percentMean,'Color',[0 0 0]);
plot(-180:30:180,actualPercent,'Color',[0 0 1]);
xticks(-180:90:180);
xlabel('\Delta direction (degrees)');
ylabel('Percentage of neurons');

figure;hold on;
hx = histogram(shuffVector);
line([actualVector actualVector],[0 120]);
xlabel('Mean \Delta orientation');
ylabel('Randomization instances');


f4 = figure;
set(f4,'Position',[25 100 1500 600]);
for c = 1:length(contrasts)
    subplot(2,length(contrasts),c);
    histogram(bestOrientShift{c}(:,5),(-180-orientations(2)/2):orientations(2):(180+orientations(2)/2),'Normalization','probability');
    xlabel('\Delta Orientation shift');
    ylabel('Percentage of neurons');
    title(['Contrast = ' num2str(contrasts(c))]);
    %     ylim([0 60]);
    
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
% saveas(f4,fullfile(saveDir,'Orientation and FR shifts'));

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
        
        if dir ==180
            orientToFR(c,1,d) = orientToFR(c,1,1);
            orientToFR(c,2,d) = orientToFR(c,2,1);
            
            orientToBestOrientFR(c,1,d) = orientToBestOrientFR(c,1,1);
            orientToBestOrientFR(c,2,d) = orientToBestOrientFR(c,2,1);
        end
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
% saveas(f5,fullfile(saveDir,'FR change at preferred orientations - scatter'));

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
% saveas(f6,fullfile(saveDir,'FR change at preferred orientations'));

% figure;
% for c = 1:5
%     subplot(1,5,c);
%     histogram(bestOrientShift{c}(:,1));
% end

for c = 1:length(contrasts)
    vv(c) = mean(bestOrientShift{c}(:,2));
    vvSTD(c) = std(bestOrientShift{c}(:,2))/sqrt(size(bestOrientShift{c}(:,2),1));
    av(c) = mean(bestOrientShift{c}(:,4));
    avSTD(c) = std(bestOrientShift{c}(:,4))/sqrt(size(bestOrientShift{c}(:,4),1));
    
    delt(c) = mean(bestOrientShift{c}(:,4) - bestOrientShift{c}(:,2));
    deltSTD(c) = std(bestOrientShift{c}(:,4) - bestOrientShift{c}(:,2))/sqrt(size(bestOrientShift{c}(:,4) - bestOrientShift{c}(:,2),1));
end
figure; hold on;

hA = area(contrasts,[vv-vvSTD; 2*vvSTD]')
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
hA = area(contrasts,[av-avSTD; 2*avSTD]')
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 1];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];

plot(contrasts,vv,'Color',[0 0 0]);
plot(contrasts,av,'Color',[0 0 1]);
xlabel('Contrast');
ylim([10 30]);
ylabel('Mean firing rate (Hz)');

figure;hold on;
hA = area(contrasts,[delt-deltSTD; 2*deltSTD]')
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.35;
hA(2).EdgeColor = [1 1 1];
plot(contrasts,delt,'Color',[0 0 1]);
xlabel('Contrast');
ylabel('\Delta firing rate (Hz)');


anovaMatSparse = [squeeze(sparsity(1,:,:))'; squeeze(sparsity(2,:,:))'];
p = anova2(anovaMatSparse,size(sparsity,3),'off');
figure;hold on;
vSpar = mean(sparsity(1,:,:),3);
avSpar = mean(sparsity(2,:,:),3);
vSparSTD = std(sparsity(1,:,:),[],3)/sqrt(size(sparsity,3));
avSparSTD = std(sparsity(2,:,:),[],3)/sqrt(size(sparsity,3));
hA = area(contrasts,[vSpar-vSparSTD; 2*vSparSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
hA = area(contrasts,[avSpar-avSparSTD; 2*avSparSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 1];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(contrasts,vSpar,'Color',[0 0 0]);
plot(contrasts,avSpar,'Color',[0 0 1]);
xlabel('Contrast');
ylabel('Sparseness');

figure;
xTime = analysisParams.analysisWindow(1)+analysisParams.frBinWidth : analysisParams.analysisWindow(2);
for c = 1:length(contrasts)
    subplot(1,5,c);hold on;
    
    meanTrainV = mean(normFRTrain{c,1,1});
    meanTrainAV = mean(normFRTrain{c,1,2});
    meanTrainA = mean(normFRTrain{c,1,3});
    stdTrainV = std(normFRTrain{c,1,1})./sqrt(size(normFRTrain{c,1,1},1));
    stdTrainAV = std(normFRTrain{c,1,2})./sqrt(size(normFRTrain{c,1,2},1));
%         hA = area(xTime,[meanTrainV-stdTrainV; 2*stdTrainV]');
%         hA(1).FaceAlpha = 0;
%         hA(1).EdgeColor = [1 1 1];
%         hA(2).FaceColor = [0 0 0];
%         hA(2).FaceAlpha = 0.5;
%         hA(2).EdgeColor = [1 1 1];
%         hA = area(xTime,[meanTrainAV-stdTrainAV; 2*stdTrainAV]');
%         hA(1).FaceAlpha = 0;
%         hA(1).EdgeColor = [1 1 1];
%         hA(2).FaceColor = [0 0 1];
%         hA(2).FaceAlpha = 0.5;
%         hA(2).EdgeColor = [1 1 1];
    plot(xTime/1000,meanTrainV,'Color',[0 0 0]);
    plot(xTime/1000,meanTrainAV,'Color',[0 0 1]);
    plot(xTime/1000,meanTrainA,'Color',[0 1 0]);
    xlabel('Time (ms)');
    ylabel('Normalized FR');
end

figure;hold on
bb = bar([mean(timingData{5,1,1}(:,1)) mean(timingData{5,1,2}(:,1))]);
errorbar([1 2],[mean(timingData{5,1,1}(:,1)) mean(timingData{5,1,2}(:,1))],[std(timingData{5,1,1}(:,1))/sqrt(length(timingData{5,1,1}(:,1))) std(timingData{5,1,2}(:,1))/sqrt(length(timingData{5,1,1}(:,1)))],'LineStyle','none','Color',[0 0 0]);
bb(1,1).FaceColor = 'flat';
bb(1,1).CData(1,:) = [0 0 0];
bb(1,1).CData(2,:) = [0 0 1];
ylabel('Latency to peak FR (ms)');

figure;hold on
bb = bar([mean(timingData{5,1,1}(:,2)) mean(timingData{5,1,2}(:,2))]);
errorbar([1 2],[mean(timingData{5,1,1}(:,2)) mean(timingData{5,1,2}(:,2))],[std(timingData{5,1,1}(:,2))/sqrt(length(timingData{5,1,1}(:,2))) std(timingData{5,1,2}(:,2))/sqrt(length(timingData{5,1,1}(:,2)))],'LineStyle','none','Color',[0 0 0]);
bb(1,1).FaceColor = 'flat';
bb(1,1).CData(1,:) = [0 0 0];
bb(1,1).CData(2,:) = [0 0 1];
ylabel('Rate of FR increase (Hz/ms)');

cof = 5;
figure;hold on
bb = bar([mean(timingData{cof,1,1}(:,3)) mean(timingData{cof,1,2}(:,3))]);
errorbar([1 2],[mean(timingData{cof,1,1}(:,3)) mean(timingData{cof,1,2}(:,3))],[std(timingData{cof,1,1}(:,3))/sqrt(length(timingData{cof,1,1}(:,3))) std(timingData{cof,1,2}(:,3))/sqrt(length(timingData{cof,1,1}(:,3)))],'LineStyle','none','Color',[0 0 0]);
bb(1,1).FaceColor = 'flat';
bb(1,1).CData(1,:) = [0 0 0];
bb(1,1).CData(2,:) = [0 0 1];
ylabel('Response FWHM (ms)');

figure;hold on
bb = bar([mean(timingData{cof,1,1}(:,4)) mean(timingData{cof,1,2}(:,4))]);
errorbar([1 2],[mean(timingData{cof,1,1}(:,4)) mean(timingData{cof,1,2}(:,4))],[std(timingData{cof,1,1}(:,4))/sqrt(length(timingData{cof,1,1}(:,4))) std(timingData{cof,1,2}(:,4))/sqrt(length(timingData{cof,1,1}(:,4)))],'LineStyle','none','Color',[0 0 0]);
bb(1,1).FaceColor = 'flat';
bb(1,1).CData(1,:) = [0 0 0];
bb(1,1).CData(2,:) = [0 0 1];
ylabel('Latency to response (ms)');

for c = 1:length(contrasts)
    latencyDiff(1,c) = mean(timingData{c,1,1}(:,4));
    latencyDiff(2,c) = mean(timingData{c,1,2}(:,4));
    latencySTD(1,c) = std(timingData{c,1,1}(:,4)) / sqrt(length(timingData{c,1,1}(:,4)));
    latencySTD(2,c) = std(timingData{c,1,2}(:,4)) / sqrt(length(timingData{c,1,2}(:,4)));
    latencyDiff(3,c) = mean(timingData{c,1,2}(:,4) - timingData{c,1,1}(:,4));
    latencySTD(3,c) = std(timingData{c,1,2}(:,4) - timingData{c,1,1}(:,4)) / sqrt(length(timingData{c,1,2}(:,4)));
    
    fwhmDiff(1,c) = mean(timingData{c,1,1}(:,3));
    fwhmDiff(2,c) = mean(timingData{c,1,2}(:,3));
    fwhmSTD(1,c) = std(timingData{c,1,1}(:,3)) / sqrt(length(timingData{c,1,1}(:,3)));
    fwhmSTD(2,c) = std(timingData{c,1,2}(:,3)) / sqrt(length(timingData{c,1,2}(:,3)));
    fwhmDiff(3,c) = mean(timingData{c,1,2}(:,3) - timingData{c,1,1}(:,3));
    fwhmSTD(3,c) = std(timingData{c,1,2}(:,3) - timingData{c,1,1}(:,3)) / sqrt(length(timingData{c,1,2}(:,3)));
end
figure; hold on;
hA = area(contrasts(2:end),[latencyDiff(1,2:end)-latencySTD(1,2:end); 2*latencySTD(1,2:end)]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
hA = area(contrasts(2:end),[latencyDiff(2,2:end)-latencySTD(2,2:end); 2*latencySTD(2,2:end)]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 1];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(contrasts(2:end),latencyDiff(1,2:end),'Color',[0 0 0]);
plot(contrasts(2:end),latencyDiff(2,2:end),'Color',[0 0 1]);
xlabel('Contrast');
ylabel('Latency (ms)');

figure; hold on;
hA = area(contrasts(2:end),[latencyDiff(3,2:end)-latencySTD(3,2:end); 2*latencySTD(3,2:end)]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(contrasts(2:end),latencyDiff(3,2:end),'Color',[0 0 0]);
xlabel('Contrast');
ylabel('\Delta latency to response (ms)');

figure; hold on;
hA = area(contrasts(2:end),[fwhmDiff(1,2:end)-fwhmSTD(1,2:end); 2*fwhmSTD(1,2:end)]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
hA = area(contrasts(2:end),[fwhmDiff(2,2:end)-fwhmSTD(2,2:end); 2*fwhmSTD(2,2:end)]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 1];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(contrasts(2:end),fwhmDiff(1,2:end),'Color',[0 0 0]);
plot(contrasts(2:end),fwhmDiff(2,2:end),'Color',[0 0 1]);
xlabel('Contrast');
ylabel('Response FWHM (ms)');

figure; hold on;
hA = area(contrasts(2:end),[fwhmDiff(3,2:end)-fwhmSTD(3,2:end); 2*fwhmSTD(3,2:end)]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(contrasts(2:end),fwhmDiff(3,2:end),'Color',[0 0 0]);
xlabel('Contrast');
ylabel('\Delta response FWHM (ms)');




figure;
for c = 1:length(contrasts)
    subplot(1,length(contrasts),c);
    hold on;
    
    newAligned{c,1}(isinf(newAligned{c,1})) = NaN;
    newAligned{c,2}(isinf(newAligned{c,2})) = NaN;
    plot(-180:30:180,nanmean(newAligned{c,1}),'Color',[0 0 0]);
    plot(-180:30:180,nanmean(newAligned{c,2}),'Color',[0 0 1]);
end


figure;
for c = 1:length(contrasts)
    subplot(1,length(contrasts),c);
    baseline = responseLinearity{1}(:,1);
    lightDiff = responseLinearity{c}(:,1) - baseline;
    soundDiff = responseLinearity{c}(:,2) - baseline;
    predictLinearity = lightDiff + soundDiff;
    actualLinearity = responseLinearity{c}(:,3) - baseline;
    predPercent = actualLinearity ./ predictLinearity;
    predPercent(predPercent<-1e10)=0/0;
    predictDiff(:,c) = [nanmean(predPercent) nanstd(predPercent)/sqrt(size(predPercent,1))]';
    scatter(predictLinearity,actualLinearity);
    title(['\Delta = ' num2str(mean(actualLinearity-predictLinearity)) 'Hz']);
    
    
    meanFR(1,c) = mean(responseLinearity{c}(:,1));
    meanFR(2,c) = mean(responseLinearity{c}(:,3));
    meanFRStd(1,c) = std(responseLinearity{c}(:,1))/sqrt(size(responseLinearity{c}(:,1),1));
    meanFRStd(2,c) = std(responseLinearity{c}(:,3))/sqrt(size(responseLinearity{c}(:,3),1));
    deltaFR(c) = mean(responseLinearity{c}(:,3) - responseLinearity{c}(:,1));
    deltaFRStd(c) = (std(responseLinearity{c}(:,3) - responseLinearity{c}(:,1)))/sqrt(size(responseLinearity{c}(:,3),1));
end
figure;
    hold on;
    hA = area(contrasts,[predictDiff(1,:)-predictDiff(2,:); 2*predictDiff(2,:)]')
    hA(1).FaceAlpha = 0;
    hA(1).EdgeColor = [1 1 1];
    hA(2).FaceColor = [0 0 0];
    hA(2).FaceAlpha = 0.35;
    hA(2).EdgeColor = [1 1 1];
    plot(contrasts,predictDiff(1,:),'Color',[0 0 1]);
    xlabel('Contrast');
    ylabel('Linearity ratio');


figure; hold on;

hA = area(contrasts,[meanFR(1,:)-meanFRStd(1,:); 2*meanFRStd(1,:)]')
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
hA = area(contrasts,[meanFR(2,:)-meanFRStd(2,:); 2*meanFRStd(2,:)]')
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 1];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];

plot(contrasts,meanFR(1,:),'Color',[0 0 0]);
plot(contrasts,meanFR(2,:),'Color',[0 0 1]);
xlabel('Contrast');
% ylim([10 30]);
ylabel('Mean firing rate (Hz)');

figure;hold on;
hA = area(contrasts,[deltaFR-deltaFRStd; 2*deltaFRStd]')
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.35;
hA(2).EdgeColor = [1 1 1];
plot(contrasts,deltaFR,'Color',[0 0 1]);
xlabel('Contrast');
ylabel('\Delta firing rate (Hz)');


