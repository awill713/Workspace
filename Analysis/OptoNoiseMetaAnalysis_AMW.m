
clear all;

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

baselineWindow = [-10 0];
onsetResponseWindow = [0 25]; %ms after sound onset
thresholdSTD = 1; %std above baseline to trigger latency

%column 1 is sound, column 2 is sound/laser
peakResponseSingle = [];
peakResponseAll = [];
latencySingle = [];
latencyAll = [];
responseDurationSingle = [];
responseDurationAll = [];


for dp = 1:length(dataPaths)
    dp
    %Load data
    dataFile = fullfile('E:\Electrophysiology\',dataPaths{dp},'OptoNoiseResponses','OptoNoiseResponseData.mat');
    if exist(dataFile)
        load(dataFile);
    else
        dataFile = fullfile('E:\Electrophysiology\',dataPaths{dp},'OptoNoiseResponses_25ms','OptoNoiseResponseData.mat');
        load(dataFile);
    end
    
    onsetBinFirst = onsetResponseWindow(1) - analysisParams.analysisWindow(1) - analysisParams.frBinWidth + 1;
    onsetBinLast = onsetResponseWindow(2) - analysisParams.analysisWindow(1) - analysisParams.frBinWidth + 1;
    baselineBinFirst = baselineWindow(1) - analysisParams.analysisWindow(1) - analysisParams.frBinWidth + 1;
    baselineBinLast = baselineWindow(2) - analysisParams.analysisWindow(1) - analysisParams.frBinWidth;
%     baselineBinLast = - analysisParams.analysisWindow(1) - analysisParams.frBinWidth;
    
    %Cycle through each neuron
    for u = 1:unitClassification.totalUnits
        
        if ismember(unitData(u).neuronNumber,unitClassification.soundResponsiveUnits)
            
            frData = unitData(u).frTrain;
            
            baselineMeanSound = mean(frData(1,baselineBinFirst:baselineBinLast));
            baselineStdSound = std(frData(1,baselineBinFirst:baselineBinLast));
            baselineMeanBoth = mean(frData(6,baselineBinFirst:baselineBinLast));
            baselineStdBoth = std(frData(6,baselineBinFirst:baselineBinLast));
            
            %peak response
            peakSound = mean(frData(1,onsetBinFirst:onsetBinLast)) - mean(frData(1,baselineBinFirst:baselineBinLast));
            peakBoth = mean(frData(6,onsetBinFirst:onsetBinLast)) - mean(frData(6,baselineBinFirst:baselineBinLast));
            
            peakResponseAll = [peakResponseAll; peakSound peakBoth];
            
            %response latency
            if peakSound < 0
                frData(1,:) = -1*(frData(1,:)-baselineMeanSound) + baselineMeanSound;
                frData(6,:) = -1*(frData(6,:)-baselineMeanBoth) + baselineMeanBoth;
            end
            
            latencySound = find(frData(1,baselineBinLast+1:end)>(baselineMeanSound+thresholdSTD*baselineStdSound),1);
            latencyBoth = find(frData(6,baselineBinLast+1:end)>(baselineMeanBoth+thresholdSTD*baselineStdBoth),1);
            
            latencyAll = [latencyAll; latencySound latencyBoth];
            
            %response duration
            aboveSound = frData(1,baselineBinLast+1:end)>(0.5*max(frData(1,baselineBinLast+1:end)));
%             aboveSound = frData(1,baselineBinLast+1:end)>(baselineMeanSound+thresholdSTD*baselineStdSound);
            aboveSoundDiff = diff(aboveSound);
            firstSound = find(aboveSoundDiff==1,1);
            lastSound = find(aboveSoundDiff==-1,1);
            durSound = lastSound-firstSound;
            
            aboveBoth = frData(6,baselineBinLast+1:end)>(0.5*max(frData(6,baselineBinLast+1:end)));
%             aboveBoth = frData(6,baselineBinLast+1:end)>(baselineMeanBoth+thresholdSTD*baselineStdBoth);
            aboveBothDiff = diff(aboveBoth);
            firstBoth = find(aboveBothDiff==1,1);
            lastBoth = find(aboveBothDiff==-1,1);
            durBoth = lastBoth-firstBoth;
            
            responseDurationAll = [responseDurationAll; durSound durBoth];
            
            if ismember(unitData(u).neuronNumber,unitClassification.singleUnits)
                peakResponseSingle = [peakResponseSingle; peakSound peakBoth];
                latencySingle = [latencySingle; latencySound latencyBoth];
                responseDurationSingle = [responseDurationSingle; durSound durBoth];
            end
        end
    end
end


figure;
subplot(1,2,1);hold on;
histogram(peakResponseAll(:,1));
histogram(peakResponseAll(:,2));
[h p] = ttest2(peakResponseAll(:,1),peakResponseAll(:,2));
title(['p = ' num2str(p)]);
subplot(1,2,2);hold on;
histogram(peakResponseSingle(:,1));
histogram(peakResponseSingle(:,2));
[h p] = ttest2(peakResponseSingle(:,1),peakResponseSingle(:,2));
title(['p = ' num2str(p)]);
suptitle('Onset response (0-25ms)')

figure;
subplot(1,2,1);hold on;
histogram(latencyAll(:,1));
histogram(latencyAll(:,2));
[h p] = ttest2(latencyAll(:,1),latencyAll(:,2));
title(['p = ' num2str(p)]);
subplot(1,2,2);hold on;
histogram(latencySingle(:,1));
histogram(latencySingle(:,2));
[h p] = ttest2(latencySingle(:,1),latencySingle(:,2));
title(['p = ' num2str(p)]);
suptitle('Response latency');

figure;
subplot(1,2,1);hold on;
histogram(responseDurationAll(:,1));
histogram(responseDurationAll(:,2));
[h p] = ttest2(responseDurationAll(:,1),responseDurationAll(:,2));
title(['p = ' num2str(p)]);
subplot(1,2,2);hold on;
histogram(responseDurationSingle(:,1));
histogram(responseDurationSingle(:,2));
[h p] = ttest2(responseDurationSingle(:,1),responseDurationSingle(:,2));
title(['p = ' num2str(p)]);
suptitle('Response duration');

        
        
        