
clear all;
close all;

savePath = fullfile('E:\Electrophysiology','EP001');

dataPaths{1} = fullfile('EP001','AW086','20190815');
% dataPaths{1} = fullfile('EP001','AW087','20190709');
dataPaths{2} = fullfile('EP001','AW087','20190709');
dataPaths{3} = fullfile('EP001','AW088','20190711');
dataPaths{4} = fullfile('EP001','AW088','20190803');
dataPaths{5} = fullfile('EP001','AW089','20190717');
dataPaths{6} = fullfile('EP001','AW089','20190818');
dataPaths{7} = fullfile('EP001','AW091','20191021');
dataPaths{8} = fullfile('EP001','AW092','20191026');
dataPaths{9} = fullfile('EP001','AW093','20191029');
dataPaths{10} = fullfile('EP001','AW093','20191104-1');
dataPaths{11} = fullfile('EP001','AW093','20191104-2');

laserOffsetIndex = 6; % [-20 -10 -5 0 5 10 20] + 2

% baselineWindow = [-40 -20];
baselineWindow = [-20 0]; %ms prior to sound onset
onsetWindow = [0 35]; %ms after sound onset
sustainWindow = [35 0]; %ms after sound onset, ms after sound offset
offsetWindow = [10 50]; %ms after sound offset
thresholdSTD = 1; %std above baseline to trigger latency

%column 1 is dataPath, column 2 is u (unitN, not unitNumber), column 3 is sound, column 4 is sound/laser
onsetAmplitudeSingle = [];
onsetAmplitudeAll = [];
onsetLatencySingle = [];
onsetLatencyAll = [];
onsetDurationSingle = [];
onsetDurationAll = [];
sustainAmplitudeSingle = [];
sustainAmplitudeAll = [];
offsetAmplitudeSingle = [];
offsetAmplitudeAll = [];

%column 1 is dataPath, column 2 is u (unitN, not unitNumber), column 3 is
%sound, column 4 is sound/laser, column 5 is whether it's sound responsive,
%column 6 is whether it's light responsive, column 7 is its supra/sublinearity
lightModulatesOnsetAll = [];
lightModulatesOnsetSingle = [];
lightModulatesSustainAll = [];
lightModulatesSustainSingle = [];
lightModulatesOffsetAll = [];
lightModulatesOffsetSingle = [];



trace = [];
totalNeurons = [];
for dp = 1:length(dataPaths)
    %Load data
    dataFile = fullfile('E:\Electrophysiology\',dataPaths{dp},'OptoNoiseResponses','OptoNoiseResponseDataFR-based.mat');
    if exist(dataFile)
        load(dataFile);
    else
        dataFile = fullfile('E:\Electrophysiology\',dataPaths{dp},'OptoNoiseResponses_25ms','OptoNoiseResponseDataFR-based.mat');
        load(dataFile);
    end
    neurons = length(unitData);
    totalNeurons = [totalNeurons neurons];
end

for dp = 1:length(dataPaths)
    dp
    %Load data
    dataFile = fullfile('E:\Electrophysiology\',dataPaths{dp},'OptoNoiseResponses','OptoNoiseResponseDataFR-based.mat');
    if exist(dataFile)
        load(dataFile);
    else
        dataFile = fullfile('E:\Electrophysiology\',dataPaths{dp},'OptoNoiseResponses_25ms','OptoNoiseResponseDataFR-based.mat');
        load(dataFile);
    end
    
    sustainWindowAdjusted = [sustainWindow(1) sustainWindow(2)+stimInfo.burstDuration];
    offsetWindowAdjusted = offsetWindow+stimInfo.burstDuration;

%     baselineScalar = 1000/(baselineWindow(2)-baselineWindow(1));
%     onsetScalar = 1000/(onsetWindow(2)-onsetWindow(1));
%     sustainWindow = 1000/(sustainWindow(2)+stimInfo.burstDuration-sustainWindow(1));
%     offsetScalar = 1000/(offsetWindow(2)-offsetWindow(1));
%     
    baselineBinFirst = baselineWindow(1) - analysisParams.analysisWindow(1) - analysisParams.frBinWidth + 1;
    baselineBinLast = baselineWindow(2) - analysisParams.analysisWindow(1) - analysisParams.frBinWidth + 1;
    onsetBinFirst = onsetWindow(1) - analysisParams.analysisWindow(1) - analysisParams.frBinWidth + 1;
    onsetBinLast = onsetWindow(2) - analysisParams.analysisWindow(1) - analysisParams.frBinWidth + 1;
    sustainBinFirst = sustainWindowAdjusted(1) - analysisParams.analysisWindow(1) - analysisParams.frBinWidth + 1;
    sustainBinLast = sustainWindowAdjusted(2) - analysisParams.analysisWindow(1) - analysisParams.frBinWidth + stimInfo.burstDuration + 1;
    offsetBinFirst = offsetWindowAdjusted(1) - analysisParams.analysisWindow(1) - analysisParams.frBinWidth + stimInfo.burstDuration + 1;
    offsetBinLast = offsetWindowAdjusted(2) - analysisParams.analysisWindow(1) - analysisParams.frBinWidth + stimInfo.burstDuration + 1;
    
    repeats = stimInfo.repeats;
    %     baselineBinLast = - analysisParams.analysisWindow(1) - analysisParams.frBinWidth;
    
    %Cycle through each neuron
    for u = 1:unitClassification.totalUnits
        if unitData(u).type==1
            trace = [trace; unitData(u).frTrain(1,:)];
        end
        if unitData(u).type~=0
            
            trialData = unitData(u).frTrainTrials;
            meanData = unitData(u).frTrain;
            spikeData = unitData(u).raster;
            spikeData{laserOffsetIndex,2} = spikeData{laserOffsetIndex,2}-(laserOffsetIndex-1)*repeats;
            spikeData{2,2} = spikeData{2,2}-repeats;
            
            %ONSET AMPLITUDE
            baselineSound = extractRaster(spikeData{1,1},spikeData{1,2},repeats,baselineWindow);
            baselineBoth = extractRaster(spikeData{laserOffsetIndex,1},spikeData{laserOffsetIndex,2},repeats,baselineWindow);
            baselineLight = extractRaster(spikeData{2,1},spikeData{2,2},repeats,baselineWindow);
            onsetSound = extractRaster(spikeData{1,1},spikeData{1,2},repeats,onsetWindow);
            onsetBoth = extractRaster(spikeData{laserOffsetIndex,1},spikeData{laserOffsetIndex,2},repeats,onsetWindow);
            onsetLight = extractRaster(spikeData{2,1},spikeData{2,2},repeats,onsetWindow);
            
            [h pp] = ttest2(onsetSound,onsetBoth,'alpha',0.05/sum(totalNeurons)); h(isnan(h))=0;
            if h
                [hS p] = ttest2(baselineSound,onsetSound); hL(isnan(hS))=0;
                [hL p] = ttest2(baselineLight,onsetLight); hL(isnan(hL))=0;
                
                lightModulatesOnsetAll = [lightModulatesOnsetAll; dp u mean(onsetSound) mean(onsetBoth) hS hL pp];

                if ismember(unitData(u).neuronNumber,unitClassification.singleUnits)
                    lightModulatesOnsetSingle = [lightModulatesOnsetSingle; dp u mean(onsetSound) mean(onsetBoth) hS hL pp];
                end
            end
            
            %SUSTAIN AMPLITUDE
            sustainSound = extractRaster(spikeData{1,1},spikeData{1,2},repeats,sustainWindowAdjusted);
            sustainBoth = extractRaster(spikeData{laserOffsetIndex,1},spikeData{laserOffsetIndex,2},repeats,sustainWindowAdjusted);
            sustainLight = extractRaster(spikeData{2,1},spikeData{2,2},repeats,sustainWindowAdjusted);
            
            [h ps] = ttest2(sustainSound,sustainBoth,'alpha',0.05/sum(totalNeurons)); h(isnan(h))=0;
            if h
                [hS p] = ttest2(baselineSound,sustainSound); hL(isnan(hS))=0;
                [hL p] = ttest2(baselineLight,sustainLight); hL(isnan(hL))=0;
                
                lightModulatesSustainAll = [lightModulatesSustainAll; dp u mean(sustainSound) mean(sustainBoth) hS hL ps];

                if ismember(unitData(u).neuronNumber,unitClassification.singleUnits)
                    lightModulatesSustainSingle = [lightModulatesSustainSingle; dp u mean(sustainSound) mean(sustainBoth) hS hL ps];
                end
            end
            
            %OFFSET AMPLITUDE
            offsetSound = extractRaster(spikeData{1,1},spikeData{1,2},repeats,offsetWindowAdjusted);
            offsetBoth = extractRaster(spikeData{laserOffsetIndex,1},spikeData{laserOffsetIndex,2},repeats,offsetWindowAdjusted);
            offsetLight = extractRaster(spikeData{2,1},spikeData{2,2},repeats,offsetWindowAdjusted);
            
            [h poff] = ttest2(offsetSound,offsetBoth,'alpha',0.05/sum(totalNeurons)); h(isnan(h))=0;
            if h
                [hS p] = ttest2(baselineSound,offsetSound); hL(isnan(hS))=0;
                [hL p] = ttest2(baselineLight,offsetLight); hL(isnan(hL))=0;
                
                lightModulatesOffsetAll = [lightModulatesOffsetAll; dp u mean(offsetSound) mean(offsetBoth) hS hL poff];

                if ismember(unitData(u).neuronNumber,unitClassification.singleUnits)
                    lightModulatesOffsetSingle = [lightModulatesOffsetSingle; dp u mean(offsetSound) mean(offsetBoth) hS hL poff];
                end
            end
            
            %latency and duration of sound responsive neurons
            [hO p] = ttest2(baselineSound,onsetSound); hO(isnan(hO))=0;
            [hSu p] = ttest2(baselineSound,sustainSound); hSu(isnan(hSu))=0;
            if hO || hSu
                if hO
                    theSign = sign(mean(onsetSound)-mean(baselineSound));
                else
                    theSign = sign(mean(sustainSound) - mean(baselineSound));
                end
                meanData(1,:) = theSign*(meanData(1,:)-mean(baselineSound)) + mean(baselineSound);
                meanData(laserOffsetIndex,:) = theSign*(meanData(laserOffsetIndex,:)-mean(baselineBoth)) + mean(baselineBoth);
                
                %latency
                baselineMeanSound = mean(meanData(1,baselineBinFirst:baselineBinLast));
                baselineStdSound = std(meanData(1,baselineBinFirst:baselineBinLast));
                baselineMeanBoth = mean(meanData(laserOffsetIndex,baselineBinFirst:baselineBinLast));
                baselineStdBoth = std(meanData(laserOffsetIndex,baselineBinFirst:baselineBinLast));
                
                latencySound = find(meanData(1,onsetBinFirst:end)>(baselineMeanSound+thresholdSTD*baselineStdSound),1);
                latencyBoth = find(meanData(laserOffsetIndex,onsetBinFirst:end)>(baselineMeanBoth+thresholdSTD*baselineStdBoth),1);
                
                if isempty(latencyBoth)
                    latencyBoth=latencySound;
                end
                
                onsetLatencyAll = [onsetLatencyAll; dp u latencySound latencyBoth];

                %duration
                aboveSound = meanData(1,onsetBinFirst:end)>(0.5*(max(meanData(1,onsetBinFirst:onsetBinLast))-baselineMeanSound)+baselineMeanSound);
                aboveSoundDiff = diff(aboveSound);
                firstSound = find(aboveSoundDiff==1,1);
                lastSound = find(aboveSoundDiff==-1,1);
                durSound = lastSound-firstSound;
                
                aboveBoth = meanData(laserOffsetIndex,onsetBinFirst:end)>(0.5*(max(meanData(laserOffsetIndex,onsetBinFirst:onsetBinLast))-baselineMeanBoth)+baselineMeanBoth);
                aboveBothDiff = diff(aboveBoth);
                firstBoth = find(aboveBothDiff==1,1);
                lastBoth = find(aboveBothDiff==-1,1);
                durBoth = lastBoth-firstBoth;
                
                if isempty(firstBoth) && ~isempty(lastBoth)
                    durBoth = lastBoth;
                end
                
                onsetDurationAll = [onsetDurationAll; dp u durSound durBoth];
                
                if ismember(unitData(u).neuronNumber,unitClassification.singleUnits)
                    onsetLatencySingle = [onsetLatencySingle; dp u latencySound latencyBoth];
                    onsetDurationSingle = [onsetDurationSingle; dp u durSound durBoth];
                end
            end
%                 
%                 
% 
% %             baselineSound = mean(trialData(1,:,baselineBinFirst:baselineBinLast),3);
% %             onsetSound = mean(trialData(1,:,onsetBinFirst:onsetBinLast),3);
%             [h p] = ttest2(baselineSound,onsetSound); h(isnan(h))=0;
%             if h
%                 
%                 %onset amplitude
%                 amplitudeSound = mean(onsetSound) - mean(baselineSound);
%                 amplitudeBoth = mean(onsetBoth) - mean(baselineBoth);
%                 
%                 onsetAmplitudeAll = [onsetAmplitudeAll; dp u amplitudeSound amplitudeBoth];
% 
% %                 amplitudeSound = onsetMeanSound - baselineMeanSound;
% %                 amplitudeBoth = onsetMeanBoth - baselineMeanBoth;
% %                 
% %                 onsetAmplitudeAll = [onsetAmplitudeAll; dp u amplitudeSound amplitudeBoth];
%                 
%                 %onset latency
%                 %temporarily flip if this is a sound-suppressed neuron
%                 meanData(1,:) = sign(amplitudeSound)*(meanData(1,:)-mean(baselineSound)) + mean(baselineSound);
%                 meanData(laserOffsetIndex,:) = sign(amplitudeSound)*(meanData(laserOffsetIndex,:)-mean(baselineBoth)) + mean(baselineBoth);
%                 
%                 baselineMeanSound = mean(meanData(1,baselineBinFirst:baselineBinLast));
%                 baselineStdSound = std(meanData(1,baselineBinFirst:baselineBinLast));
%                 baselineMeanBoth = mean(meanData(laserOffsetIndex,baselineBinFirst:baselineBinLast));
%                 baselineStdBoth = std(meanData(laserOffsetIndex,baselineBinFirst:baselineBinLast));
%                 
% %                 onsetMeanSound = mean(meanData(1,onsetBinFirst:onsetBinLast));
% %                 onsetMeanBoth = mean(mean(laserOffsetIndex,onsetBinFirst:onsetBinLast));
%                 latencySound = find(meanData(1,onsetBinFirst:end)>(baselineMeanSound+thresholdSTD*baselineStdSound),1);
%                 latencyBoth = find(meanData(laserOffsetIndex,onsetBinFirst:end)>(baselineMeanBoth+thresholdSTD*baselineStdBoth),1);
%                 
%                 if isempty(latencyBoth)
%                     latencyBoth=latencySound;
%                 end
%                 
%                 onsetLatencyAll = [onsetLatencyAll; dp u latencySound latencyBoth];
%                 
%                 %onset duration
%                 aboveSound = meanData(1,onsetBinFirst:end)>(0.5*(max(meanData(1,onsetBinFirst:onsetBinLast))-baselineMeanSound)+baselineMeanSound);
%                 aboveSoundDiff = diff(aboveSound);
%                 firstSound = find(aboveSoundDiff==1,1);
%                 lastSound = find(aboveSoundDiff==-1,1);
%                 durSound = lastSound-firstSound;
%                 
%                 aboveBoth = meanData(laserOffsetIndex,onsetBinFirst:end)>(0.5*(max(meanData(laserOffsetIndex,onsetBinFirst:onsetBinLast))-baselineMeanBoth)+baselineMeanBoth);
%                 aboveBothDiff = diff(aboveBoth);
%                 firstBoth = find(aboveBothDiff==1,1);
%                 lastBoth = find(aboveBothDiff==-1,1);
%                 durBoth = lastBoth-firstBoth;
%                 
%                 if isempty(firstBoth) && ~isempty(lastBoth)
%                     durBoth = lastBoth;
%                 end
%                 onsetDurationAll = [onsetDurationAll; dp u durSound durBoth];
%                 
%                 if ismember(unitData(u).neuronNumber,unitClassification.singleUnits)
%                     onsetAmplitudeSingle = [onsetAmplitudeSingle; dp u amplitudeSound amplitudeBoth];
%                     onsetLatencySingle = [onsetLatencySingle; dp u latencySound latencyBoth];
%                     onsetDurationSingle = [onsetDurationSingle; dp u durSound durBoth];
%                 end
%                 
%                 %return to original orientation (see onset latency section)
%                 meanData(1,:) = sign(amplitudeSound)*(meanData(1,:)-mean(baselineSound)) + mean(baselineSound);
%                 meanData(laserOffsetIndex,:) = sign(amplitudeSound)*(meanData(laserOffsetIndex,:)-mean(baselineBoth)) + mean(baselineBoth);
%             end
%             
%             %SUSTAIN RESPONSE
%             sustainSound = extractRaster(spikeData{1,1},spikeData{1,2},repeats,sustainWindowAdjusted);
% %             sustainSound = mean(trialData(1,:,sustainBinFirst:sustainBinLast),3);
%             [h p] = ttest2(baselineSound,sustainSound); h(isnan(h))=0;
%             if h
%                 baselineBoth = extractRaster(spikeData{laserOffsetIndex,1},spikeData{laserOffsetIndex,2},repeats,baselineWindow);
%                 sustainBoth = extractRaster(spikeData{laserOffsetIndex,1},spikeData{laserOffsetIndex,2},repeats,sustainWindowAdjusted);
%                 
% %                 baselineMeanSound = mean(baselineSound);
% %                 baselineMeanBoth = mean(mean(trialData(laserOffsetIndex,:,baselineBinFirst:baselineBinLast),3));
% %                 
% %                 sustainMeanSound = mean(sustainSound);
% %                 sustainMeanBoth = mean(mean(trialData(laserOffsetIndex,:,sustainBinFirst:sustainBinLast),3));
% %                 
%                 %sustain amplitude
%                 amplitudeSound = mean(sustainSound) - mean(baselineSound);
%                 amplitudeBoth = mean(sustainBoth) - mean(baselineBoth);
%                 
% %                 amplitudeSound = sustainMeanSound - baselineMeanBoth;
% %                 amplitudeBoth = sustainMeanBoth - baselineMeanBoth;
%                 
%                 sustainAmplitudeAll = [sustainAmplitudeAll; dp u amplitudeSound amplitudeBoth];
%                 
%                 if ismember(unitData(u).neuronNumber,unitClassification.singleUnits)
%                     sustainAmplitudeSingle = [sustainAmplitudeSingle; dp u amplitudeSound amplitudeBoth];
%                 end
%             end
%             
%             %OFFSET RESPONSE
%             offsetSound = extractRaster(spikeData{1,1},spikeData{1,2},repeats,offsetWindowAdjusted);
% %             offsetSound = mean(trialData(1,:,offsetBinFirst:offsetBinLast),3);
%             [h p] = ttest2(baselineSound,offsetSound); h(isnan(h))=0;
%             if h
%                 baselineBoth = extractRaster(spikeData{laserOffsetIndex,1},spikeData{laserOffsetIndex,2},repeats,baselineWindow);
%                 offsetBoth = extractRaster(spikeData{laserOffsetIndex,1},spikeData{laserOffsetIndex,2},repeats,offsetWindowAdjusted);
%                 
% %                 baselineMeanSound = mean(baselineSound);
% %                 %             baselineStdSound = std(baselineSound);
% %                 baselineMeanBoth = mean(mean(trialData(laserOffsetIndex,:,baselineBinFirst:baselineBinLast),3));
% %                 %             baselineStdBoth = std(mean(trialData(laserOffsetIndex,:,baselineBinFirst:baselineBinLast),3));
% %                 
% %                 offsetMeanSound = mean(offsetSound);
% %                 offsetMeanBoth = mean(mean(trialData(laserOffsetIndex,:,offsetBinFirst:offsetBinLast),3));
% %                 
%                 %offset amplitude
%                 amplitudeSound = mean(offsetSound) - mean(baselineSound);
%                 amplitudeBoth = mean(offsetBoth) - mean(baselineBoth);
% %                 amplitudeSound = offsetMeanSound - baselineMeanSound;
% %                 amplitudeBoth = offsetMeanBoth - baselineMeanBoth;
%                 
%                 offsetAmplitudeAll = [offsetAmplitudeAll; dp u amplitudeSound amplitudeBoth];
%                 
%                 if ismember(unitData(u).neuronNumber,unitClassification.singleUnits)
%                     offsetAmplitudeSingle = [offsetAmplitudeSingle; dp u amplitudeSound amplitudeBoth];
%                 end
%                 
%                 %             %offset latency
%                 %             %temporarily flip if this is a sound-suppressed neuron
%                 %             meanData(1,:) = sign(amplitudeSound)*(meanData(1,:)-baselineMeanSound) + baselineMeanSound;
%                 %             meanData(laserOffsetIndex,:) = sign(amplitudeSound)*(meanData(laserOffsetIndex,:)-baselineMeanBoth) + baselineMeanBoth;
%                 %
%                 %             latencySound = find(meanData(1,offsetBinFirst:end)>(baselineMeanSound+thresholdSTD*baselineStdSound),1);
%                 %             latencyBoth = find(meanData(laserOffsetIndex,offsetBinFirst:end)>(baselineMeanBoth+thresholdSTD*baselineStdSound),1);
%                 %
%                 %             offsetLatencyAll = [offsetLatencyAll; latencySound latencyBoth];
%                 
%                 %             %offset duration
%                 %             aboveSound = meanData(1,offsetBinFirst:end)>(0.5*(max(meanData(1,offsetBinFirst:offsetBinLast))-baselineMeanSound)+sustainMeanSound);
%                 %             aboveSoundDiff = diff(aboveSound);
%                 %             firstSound = find(aboveSoundDiff==1,1);
%                 %             lastSound = find(aboveSoundDiff==-1,1);
%                 %             durSound = lastSound-firstSound + offsetWindow(1);
%                 %
%                 %             aboveBoth = meanData(laserOffsetIndex,offsetBinFirst:end)>(0.5*(max(meanData(laserOffsetIndex,offsetBinFirst:offsetBinLast))-baselineMeanBoth)+baselineMeanBoth);
%                 %             aboveBothDiff = diff(aboveBoth);
%                 %             firstBoth = find(aboveBothDiff==1,1);
%                 %             lastBoth = find(aboveBothDiff==-1,1);
%                 %             durBoth = lastBoth-firstBoth + offsetWindow(1);
%                 %
%                 %             offsetDurationAll = [offsetDurationAll; durSound durBoth];
%                 %
%                 %             if ismember(unitData(u).neuronNumber,unitClassification.singleUnits)
%                 %                 offsetAmplitudeSingle = [offsetAmplitudeSingle; amplitudeSound amplitudeBoth];
%                 %                 offsetLatencySingle = [offsetLatencySingle; latencySound latencyBoth];
%                 %                 offsetDurationSingle = [offsetDurationSingle; durSound durBoth];
%                 %             end
%                 %
%                 %             %return to original orientation (see offset latency section)
%                 %             meanData(1,:) = sign(amplitudeSound)*(meanData(1,:)-baselineMeanSound) + baselineMeanSound;
%                 %             meanData(laserOffsetIndex,:) = sign(amplitudeSound)*(meanData(laserOffsetIndex,:)-baselineMeanBoth) + baselineMeanBoth;
%             end
%             
        end
    end
end

f1 = figure;
subplot(2,2,1);hold on;
if size(lightModulatesOnsetAll,1)~=0
    histogram(lightModulatesOnsetAll(:,3),'BinWidth',5);
    histogram(lightModulatesOnsetAll(:,4),'BinWidth',5);
    [h p] = ttest(lightModulatesOnsetAll(:,3),lightModulatesOnsetAll(:,4));
    title(['p = ' num2str(p)]);
end
subplot(2,2,2);hold on;
if size(lightModulatesOnsetSingle,1)~=0
    histogram(lightModulatesOnsetSingle(:,3),'BinWidth',5);
    histogram(lightModulatesOnsetSingle(:,4),'BinWidth',5);
    [h p] = ttest(lightModulatesOnsetSingle(:,3),lightModulatesOnsetSingle(:,4));
    title(['p = ' num2str(p)]);
end
subplot(2,2,3);hold on;
bar([1 2],mean(lightModulatesOnsetAll(:,3:4)));
xticks([1 2]);
xticklabels({'Sound','Sound/Laser'});
er = errorbar([1 2],mean(lightModulatesOnsetAll(:,3:4)),std(lightModulatesOnsetAll(:,3:4),[],1));
er.Color = [0 0 0];
er.LineStyle = 'none';
subplot(2,2,4);hold on;
bar([1 2],mean(lightModulatesOnsetSingle(:,3:4),1));
xticks([1 2]);
xticklabels({'Sound','Sound/Laser'});
er = errorbar([1 2],mean(lightModulatesOnsetSingle(:,3:4),1),std(lightModulatesOnsetSingle(:,3:4),[],1));
er.Color = [0 0 0];
er.LineStyle = 'none';
suptitle('Onset amplitude')



f2 = figure;
subplot(2,2,1);hold on;
if size(onsetLatencyAll,1)~=0
    histogram(onsetLatencyAll(:,3),'BinWidth',5);
    histogram(onsetLatencyAll(:,4),'BinWidth',5);
    [p h] = signrank(onsetLatencyAll(:,3),onsetLatencyAll(:,4));
    title(['p = ' num2str(p)]);
end
subplot(2,2,2);hold on;
if size(onsetLatencySingle,1)~=0
    histogram(onsetLatencySingle(:,3),'BinWidth',5);
    histogram(onsetLatencySingle(:,4),'BinWidth',5);
    [p h] = signrank(onsetLatencySingle(:,3),onsetLatencySingle(:,4));
    title(['p = ' num2str(p)]);
end
subplot(2,2,3);hold on;
bar([1 2],mean(onsetLatencyAll(:,3:4)));
xticks([1 2]);
xticklabels({'Sound','Sound/Laser'});
er = errorbar([1 2],mean(onsetLatencyAll(:,3:4)),std(onsetLatencyAll(:,3:4),[],1));
er.Color = [0 0 0];
er.LineStyle = 'none';
subplot(2,2,4);hold on;
bar([1 2],mean(onsetLatencySingle(:,3:4)));
xticks([1 2]);
xticklabels({'Sound','Sound/Laser'});
er = errorbar([1 2],mean(onsetLatencySingle(:,3:4)),std(onsetLatencySingle(:,3:4),[],1));
er.Color = [0 0 0];
er.LineStyle = 'none';
suptitle('Onset response latency');

f3 = figure;
subplot(2,2,1);hold on;
if size(onsetDurationAll,1)~=0
    histogram(onsetDurationAll(:,3),'BinWidth',5);
    histogram(onsetDurationAll(:,4),'BinWidth',5);
    [p h] = signrank(onsetDurationAll(:,3),onsetDurationAll(:,4));
    title(['p = ' num2str(p)]);
end
subplot(2,2,2);hold on;
if size(onsetDurationSingle,1)~=0
    histogram(onsetDurationSingle(:,3),'BinWidth',5);
    histogram(onsetDurationSingle(:,4),'BinWidth',5);
    [p h] = signrank(onsetDurationSingle(:,3),onsetDurationSingle(:,4));
    title(['p = ' num2str(p)]);
end
subplot(2,2,3);hold on;
bar([1 2],mean(onsetDurationAll(:,3:4)));
xticks([1 2]);
xticklabels({'Sound','Sound/Laser'});
er = errorbar([1 2],mean(onsetDurationAll(:,3:4)),std(onsetDurationAll(:,3:4),[],1));
er.Color = [0 0 0];
er.LineStyle = 'none';
subplot(2,2,4);hold on;
bar([1 2],mean(onsetDurationSingle(:,3:4)));
xticks([1 2]);
xticklabels({'Sound','Sound/Laser'});
er = errorbar([1 2],mean(onsetDurationSingle(:,3:4)),std(onsetDurationSingle(:,3:4),[],1));
er.Color = [0 0 0];
er.LineStyle = 'none';
suptitle('Onset response duration');

f4 = figure;
subplot(2,2,1);hold on;
if size(lightModulatesSustainAll,1)~=0
    histogram(lightModulatesSustainAll(:,3),'BinWidth',5);
    histogram(lightModulatesSustainAll(:,4),'BinWidth',5);
    [h p] = ttest(lightModulatesSustainAll(:,3),lightModulatesSustainAll(:,4));
    title(['p = ' num2str(p)]);
end
subplot(2,2,2);hold on;
if size(lightModulatesSustainSingle,1)~=0
    histogram(lightModulatesSustainSingle(:,3),'BinWidth',5);
    histogram(lightModulatesSustainSingle(:,4),'BinWidth',5);
    [h p] = ttest(lightModulatesSustainSingle(:,3),lightModulatesSustainSingle(:,4));
    title(['p = ' num2str(p)]);
end
subplot(2,2,3);hold on;
bar([1 2],mean(lightModulatesSustainAll(:,3:4)));
xticks([1 2]);
xticklabels({'Sound','Sound/Laser'});
er = errorbar([1 2],mean(lightModulatesSustainAll(:,3:4)),std(lightModulatesSustainAll(:,3:4),[],1));
er.Color = [0 0 0];
er.LineStyle = 'none';
subplot(2,2,4);hold on;
bar([1 2],mean(lightModulatesSustainSingle(:,3:4)));
xticks([1 2]);
xticklabels({'Sound','Sound/Laser'});
er = errorbar([1 2],mean(lightModulatesSustainSingle(:,3:4)),std(lightModulatesSustainSingle(:,3:4),[],1));
er.Color = [0 0 0];
er.LineStyle = 'none';
suptitle('Sustain response')

f5 = figure;
subplot(2,2,1);hold on;
if size(lightModulatesOffsetAll,1)~=0
    histogram(lightModulatesOffsetAll(:,3),'BinWidth',5);
    histogram(lightModulatesOffsetAll(:,4),'BinWidth',5);
    [h p] = ttest(lightModulatesOffsetAll(:,3),lightModulatesOffsetAll(:,4));
    title(['p = ' num2str(p)]);
end
subplot(2,2,2);hold on;
if size(lightModulatesOffsetSingle,1)~=0
    histogram(lightModulatesOffsetSingle(:,3),'BinWidth',5);
    histogram(lightModulatesOffsetSingle(:,4),'BinWidth',5);
    [h p] = ttest(lightModulatesOffsetSingle(:,3),lightModulatesOffsetSingle(:,4));
    title(['p = ' num2str(p)]);
end
subplot(2,2,3);hold on;
bar([1 2],mean(lightModulatesOffsetAll(:,3:4)));
xticks([1 2]);
xticklabels({'Sound','Sound/Laser'});
er = errorbar([1 2],mean(lightModulatesOffsetAll(:,3:4)),std(lightModulatesOffsetAll(:,3:4),[],1));
er.Color = [0 0 0];
er.LineStyle = 'none';
subplot(2,2,4);hold on;
bar([1 2],mean(lightModulatesOffsetSingle(:,3:4)));
xticks([1 2]);
xticklabels({'Sound','Sound/Laser'});
er = errorbar([1 2],mean(lightModulatesOffsetSingle(:,3:4)),std(lightModulatesOffsetSingle(:,3:4),[],1));
er.Color = [0 0 0];
er.LineStyle = 'none';
suptitle('Offset amplitude')


analysisParams.laserOffsetIndex = laserOffsetIndex;
analysisParams.baselineWindow = baselineWindow;
analysisParams.onsetWindow = onsetWindow;
analysisParams.sustainWindow = sustainWindowAdjusted;
analysisParams.offsetWindow = offsetWindowAdjusted;
analysisParams.thresholdSTD = thresholdSTD;

responseData.onsetAmplitudeSingle = onsetAmplitudeSingle;
responseData.onsetAmplitudeAll = onsetAmplitudeAll;
responseData.onsetLatencySingle = onsetLatencySingle;
responseData.onsetLatencyAll = onsetLatencyAll;
responseData.onsetDurationSingle = onsetDurationSingle;
responseData.onsetDurationAll = onsetDurationAll;
responseData.sustainAmplitudeSingle = sustainAmplitudeSingle;
responseData.sustainAmplitudeAll = sustainAmplitudeAll;
responseData.offsetAmplitudeSingle = offsetAmplitudeSingle;
responseData.offsetAmplitudeAll = offsetAmplitudeAll;

% saveas(f1,fullfile(savePath,'Onset amplitude figure'));
% saveas(f2,fullfile(savePath,'Onset latency figure'));
% saveas(f3,fullfile(savePath,'Onset duration figure'));
% saveas(f4,fullfile(savePath,'Sustained response amplitude figure'));
% saveas(f5,fullfile(savePath,'Offset amplitude figure'));
% save(fullfile(savePath,'OptoNoiseMetaData.mat'),'dataPaths','analysisParams','responseData');


function trialsFR = extractRaster(xVector,yVector,r,window)
    edges = 0.5:1:(r+0.5);
    frScalar = 1000/(window(2)-window(1));
    
    windowSpikeIndx = find(xVector>window(1) & xVector<window(2));
    trials = yVector(windowSpikeIndx);
    trialsFR = histcounts(trials,edges)'*frScalar;
end
