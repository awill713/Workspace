
clear all;
close all;

savePath = fullfile('E:\Electrophysiology','EP001');

dataPaths{1} = fullfile('EP001','AW086','20190815');
% dataPaths{1} = fullfile('EP001','AW087','20190709');
dataPaths{2} = fullfile('EP001','AW087','20190709');
dataPaths{3} = fullfile('EP001','AW088','20190803');
dataPaths{4} = fullfile('EP001','AW089','20190717');
dataPaths{5} = fullfile('EP001','AW089','20190818');
dataPaths{6} = fullfile('EP001','AW091','20191021');
dataPaths{7} = fullfile('EP001','AW092','20191026');
dataPaths{8} = fullfile('EP001','AW093','20191029');
dataPaths{9} = fullfile('EP001','AW093','20191104-1');
dataPaths{10} = fullfile('EP001','AW093','20191104-2');

laserOffsetIndex = 3; % [-20 -10 -5 0 5 10 20] + 2

baselineWindow = [-40 -20];
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
% offsetLatencySingle = [];
% offsetLatencyAll = [];
% offsetDurationSingle = [];
% offsetDurationAll = [];


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
    
    baselineBinFirst = baselineWindow(1) - analysisParams.analysisWindow(1) - analysisParams.frBinWidth + 1;
    baselineBinLast = baselineWindow(2) - analysisParams.analysisWindow(1) - analysisParams.frBinWidth + 1;
    onsetBinFirst = onsetWindow(1) - analysisParams.analysisWindow(1) - analysisParams.frBinWidth + 1;
    onsetBinLast = onsetWindow(2) - analysisParams.analysisWindow(1) - analysisParams.frBinWidth + 1;
    sustainBinFirst = sustainWindow(1) - analysisParams.analysisWindow(1) - analysisParams.frBinWidth + 1;
    sustainBinLast = sustainWindow(2) - analysisParams.analysisWindow(1) - analysisParams.frBinWidth + stimInfo.burstDuration + 1;
    offsetBinFirst = offsetWindow(1) - analysisParams.analysisWindow(1) - analysisParams.frBinWidth + stimInfo.burstDuration + 1;
    offsetBinLast = offsetWindow(2) - analysisParams.analysisWindow(1) - analysisParams.frBinWidth + stimInfo.burstDuration + 1;
    
    %     baselineBinLast = - analysisParams.analysisWindow(1) - analysisParams.frBinWidth;
    
    %Cycle through each neuron
    for u = 1:unitClassification.totalUnits
        if unitData(u).type~=0
            
            trialData = unitData(u).frTrainTrials;
            meanData = unitData(u).frTrain;
            
            %ONSET RESPONSE
            baselineSound = mean(trialData(1,:,baselineBinFirst:baselineBinLast),3);
            onsetSound = mean(trialData(1,:,onsetBinFirst:onsetBinLast),3);
            [h p] = ttest2(baselineSound,onsetSound); h(isnan(h))=0;
            if h
                baselineMeanSound = mean(baselineSound);
                baselineStdSound = std(meanData(1,:));
                baselineMeanBoth = mean(mean(trialData(laserOffsetIndex,:,baselineBinFirst:baselineBinLast),3));
                baselineStdBoth = std(meanData(laserOffsetIndex,:));
                
                onsetMeanSound = mean(onsetSound);
                onsetMeanBoth = mean(mean(trialData(laserOffsetIndex,:,onsetBinFirst:onsetBinLast),3));
                
                %onset amplitude
                amplitudeSound = onsetMeanSound - baselineMeanSound;
                amplitudeBoth = onsetMeanBoth - baselineMeanBoth;
                
                onsetAmplitudeAll = [onsetAmplitudeAll; dp u amplitudeSound amplitudeBoth];
                
                %onset latency
                %temporarily flip if this is a sound-suppressed neuron
                meanData(1,:) = sign(amplitudeSound)*(meanData(1,:)-baselineMeanSound) + baselineMeanSound;
                meanData(laserOffsetIndex,:) = sign(amplitudeSound)*(meanData(laserOffsetIndex,:)-baselineMeanBoth) + baselineMeanBoth;
                
                latencySound = find(meanData(1,onsetBinFirst:end)>(baselineMeanSound+thresholdSTD*baselineStdSound),1);
                latencyBoth = find(meanData(laserOffsetIndex,onsetBinFirst:end)>(baselineMeanBoth+thresholdSTD*baselineStdSound),1);
                
                if isempty(latencyBoth)
                    latencyBoth=latencySound;
                end
                
                onsetLatencyAll = [onsetLatencyAll; dp u latencySound latencyBoth];
                
                %onset duration
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
                    onsetAmplitudeSingle = [onsetAmplitudeSingle; dp u amplitudeSound amplitudeBoth];
                    onsetLatencySingle = [onsetLatencySingle; dp u latencySound latencyBoth];
                    onsetDurationSingle = [onsetDurationSingle; dp u durSound durBoth];
                end
                
                %return to original orientation (see onset latency section)
                meanData(1,:) = sign(amplitudeSound)*(meanData(1,:)-baselineMeanSound) + baselineMeanSound;
                meanData(laserOffsetIndex,:) = sign(amplitudeSound)*(meanData(laserOffsetIndex,:)-baselineMeanBoth) + baselineMeanBoth;
            end
            
            %SUSTAIN RESPONSE
            sustainSound = mean(trialData(1,:,sustainBinFirst:sustainBinLast),3);
            [h p] = ttest2(baselineSound,sustainSound); h(isnan(h))=0;
            if h
                baselineMeanSound = mean(baselineSound);
                baselineMeanBoth = mean(mean(trialData(laserOffsetIndex,:,baselineBinFirst:baselineBinLast),3));
                
                sustainMeanSound = mean(sustainSound);
                sustainMeanBoth = mean(mean(trialData(laserOffsetIndex,:,sustainBinFirst:sustainBinLast),3));
                
                %sustain amplitude
                amplitudeSound = sustainMeanSound - baselineMeanBoth;
                amplitudeBoth = sustainMeanBoth - baselineMeanBoth;
                
                sustainAmplitudeAll = [sustainAmplitudeAll; dp u amplitudeSound amplitudeBoth];
                
                if ismember(unitData(u).neuronNumber,unitClassification.singleUnits)
                    sustainAmplitudeSingle = [sustainAmplitudeSingle; dp u amplitudeSound amplitudeBoth];
                end
            end
            
            %OFFSET RESPONSE
            offsetSound = mean(trialData(1,:,offsetBinFirst:offsetBinLast),3);
            [h p] = ttest2(baselineSound,offsetSound); h(isnan(h))=0;
            if h
                baselineMeanSound = mean(baselineSound);
                %             baselineStdSound = std(baselineSound);
                baselineMeanBoth = mean(mean(trialData(laserOffsetIndex,:,baselineBinFirst:baselineBinLast),3));
                %             baselineStdBoth = std(mean(trialData(laserOffsetIndex,:,baselineBinFirst:baselineBinLast),3));
                
                offsetMeanSound = mean(offsetSound);
                offsetMeanBoth = mean(mean(trialData(laserOffsetIndex,:,offsetBinFirst:offsetBinLast),3));
                
                %offset amplitude
                amplitudeSound = offsetMeanSound - baselineMeanSound;
                amplitudeBoth = offsetMeanBoth - baselineMeanBoth;
                
                offsetAmplitudeAll = [offsetAmplitudeAll; dp u amplitudeSound amplitudeBoth];
                
                if ismember(unitData(u).neuronNumber,unitClassification.singleUnits)
                    offsetAmplitudeSingle = [offsetAmplitudeSingle; dp u amplitudeSound amplitudeBoth];
                end
                
                %             %offset latency
                %             %temporarily flip if this is a sound-suppressed neuron
                %             meanData(1,:) = sign(amplitudeSound)*(meanData(1,:)-baselineMeanSound) + baselineMeanSound;
                %             meanData(laserOffsetIndex,:) = sign(amplitudeSound)*(meanData(laserOffsetIndex,:)-baselineMeanBoth) + baselineMeanBoth;
                %
                %             latencySound = find(meanData(1,offsetBinFirst:end)>(baselineMeanSound+thresholdSTD*baselineStdSound),1);
                %             latencyBoth = find(meanData(laserOffsetIndex,offsetBinFirst:end)>(baselineMeanBoth+thresholdSTD*baselineStdSound),1);
                %
                %             offsetLatencyAll = [offsetLatencyAll; latencySound latencyBoth];
                
                %             %offset duration
                %             aboveSound = meanData(1,offsetBinFirst:end)>(0.5*(max(meanData(1,offsetBinFirst:offsetBinLast))-baselineMeanSound)+sustainMeanSound);
                %             aboveSoundDiff = diff(aboveSound);
                %             firstSound = find(aboveSoundDiff==1,1);
                %             lastSound = find(aboveSoundDiff==-1,1);
                %             durSound = lastSound-firstSound + offsetWindow(1);
                %
                %             aboveBoth = meanData(laserOffsetIndex,offsetBinFirst:end)>(0.5*(max(meanData(laserOffsetIndex,offsetBinFirst:offsetBinLast))-baselineMeanBoth)+baselineMeanBoth);
                %             aboveBothDiff = diff(aboveBoth);
                %             firstBoth = find(aboveBothDiff==1,1);
                %             lastBoth = find(aboveBothDiff==-1,1);
                %             durBoth = lastBoth-firstBoth + offsetWindow(1);
                %
                %             offsetDurationAll = [offsetDurationAll; durSound durBoth];
                %
                %             if ismember(unitData(u).neuronNumber,unitClassification.singleUnits)
                %                 offsetAmplitudeSingle = [offsetAmplitudeSingle; amplitudeSound amplitudeBoth];
                %                 offsetLatencySingle = [offsetLatencySingle; latencySound latencyBoth];
                %                 offsetDurationSingle = [offsetDurationSingle; durSound durBoth];
                %             end
                %
                %             %return to original orientation (see offset latency section)
                %             meanData(1,:) = sign(amplitudeSound)*(meanData(1,:)-baselineMeanSound) + baselineMeanSound;
                %             meanData(laserOffsetIndex,:) = sign(amplitudeSound)*(meanData(laserOffsetIndex,:)-baselineMeanBoth) + baselineMeanBoth;
            end
            
        end
    end
end


f1 = figure;
subplot(1,2,1);hold on;
if size(onsetAmplitudeAll,1)~=0
    histogram(onsetAmplitudeAll(:,3));
    histogram(onsetAmplitudeAll(:,4));
    [h p] = ttest(onsetAmplitudeAll(:,3),onsetAmplitudeAll(:,4));
    title(['p = ' num2str(p)]);
end
subplot(1,2,2);hold on;
if size(onsetAmplitudeSingle,1)~=0
    histogram(onsetAmplitudeSingle(:,3));
    histogram(onsetAmplitudeSingle(:,4));
    [h p] = ttest(onsetAmplitudeSingle(:,3),onsetAmplitudeSingle(:,4));
    title(['p = ' num2str(p)]);
end
suptitle('Onset amplitude (0-25ms)')

f2 = figure;
subplot(1,2,1);hold on;
if size(onsetLatencyAll,1)~=0
    histogram(onsetLatencyAll(:,3));
    histogram(onsetLatencyAll(:,4));
    [p h] = signrank(onsetLatencyAll(:,3),onsetLatencyAll(:,4));
    title(['p = ' num2str(p)]);
end
subplot(1,2,2);hold on;
if size(onsetLatencySingle,1)~=0
    histogram(onsetLatencySingle(:,3));
    histogram(onsetLatencySingle(:,4));
    [p h] = signrank(onsetLatencySingle(:,3),onsetLatencySingle(:,4));
    title(['p = ' num2str(p)]);
end
suptitle('Onset response latency');

f3 = figure;
subplot(1,2,1);hold on;
if size(onsetDurationAll,1)~=0
    histogram(onsetDurationAll(:,3));
    histogram(onsetDurationAll(:,4));
    [p h] = signrank(onsetDurationAll(:,3),onsetDurationAll(:,4));
    title(['p = ' num2str(p)]);
end
subplot(1,2,2);hold on;
if size(onsetDurationSingle,1)~=0
    histogram(onsetDurationSingle(:,3));
    histogram(onsetDurationSingle(:,4));
    [p h] = signrank(onsetDurationSingle(:,3),onsetDurationSingle(:,4));
    title(['p = ' num2str(p)]);
end
suptitle('Onset response duration');

f4 = figure;
subplot(1,2,1);hold on;
if size(sustainAmplitudeAll,1)~=0
    histogram(sustainAmplitudeAll(:,3));
    histogram(sustainAmplitudeAll(:,4));
    [h p] = ttest(sustainAmplitudeAll(:,3),sustainAmplitudeAll(:,4));
    title(['p = ' num2str(p)]);
end
subplot(1,2,2);hold on;
if size(sustainAmplitudeSingle,1)~=0
    histogram(sustainAmplitudeSingle(:,3));
    histogram(sustainAmplitudeSingle(:,4));
    [h p] = ttest(sustainAmplitudeSingle(:,3),sustainAmplitudeSingle(:,4));
    title(['p = ' num2str(p)]);
end
suptitle('Sustain response (25ms - end of tone)')

f5 = figure;
subplot(1,2,1);hold on;
if size(offsetAmplitudeAll,1)~=0
    histogram(offsetAmplitudeAll(:,3));
    histogram(offsetAmplitudeAll(:,4));
    [h p] = ttest(offsetAmplitudeAll(:,3),offsetAmplitudeAll(:,4));
    title(['p = ' num2str(p)]);
end
subplot(1,2,2);hold on;
if size(offsetAmplitudeSingle,1)~=0
    histogram(offsetAmplitudeSingle(:,3));
    histogram(offsetAmplitudeSingle(:,4));
    [h p] = ttest(offsetAmplitudeSingle(:,3),offsetAmplitudeSingle(:,4));
    title(['p = ' num2str(p)]);
end
suptitle('Offset amplitude (0-25ms)')


analysisParams.laserOffsetIndex = laserOffsetIndex;
analysisParams.baselineWindow = baselineWindow;
analysisParams.onsetWindow = onsetWindow;
analysisParams.sustainWindow = sustainWindow;
analysisParams.offsetWindow = offsetWindow;
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

saveas(f1,fullfile(savePath,'Onset amplitude figure'));
saveas(f2,fullfile(savePath,'Onset latency figure'));
saveas(f3,fullfile(savePath,'Onset duration figure'));
saveas(f4,fullfile(savePath,'Sustained response amplitude figure'));
saveas(f5,fullfile(savePath,'Offset amplitude figure'));
save(fullfile(savePath,'OptoNoiseMetaData.mat'),'dataPaths','analysisParams','responseData');
