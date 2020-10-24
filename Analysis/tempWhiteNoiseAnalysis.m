
clear all;
% close all;

saveDir = fullfile('D:\Electrophysiology','EP007','Muscimol analysis');
if ~exist(saveDir)
    mkdir(saveDir);
end

dataPaths{1} = fullfile('EP007','AW137','20200826-1');
dataPaths{2} = fullfile('EP007','AW137','20200826-2');
dataPaths{3} = fullfile('EP007','AW138','20200826-1');
dataPaths{4} = fullfile('EP007','AW138','20200826-2');
dataPaths{5} = fullfile('EP007','AW139','20200827-1');
dataPaths{6} = fullfile('EP007','AW139','20200827-2');

for dp = 1:length(dataPaths)
    dp
    dataFile = fullfile('D:\Electrophysiology\',dataPaths{dp},'WhiteNoiseMultidB','WhiteNoiseMultidB_Data.mat');
    load(dataFile);
    
    metaData{dp}.totalUnits = length(responsiveUnits.multiUnits) + length(responsiveUnits.singleUnits);
    metaData{dp}.responsiveUnits = length(responsiveUnits.soundResponsiveUnits);
    
    baselineFR = [];
    ampTrace = [];
    psth = [];
    for u = 1:length(unitData)
        if unitData(u).type~=0
            baselineFR = [baselineFR unitData(u).meanResponse(1,4)];
        end
        if ismember(unitData(u).neuronNumber,responsiveUnits.soundResponsiveUnits)
            ampTrace = [ampTrace; unitData(u).meanResponse(:,4)'];
            psth = [psth; unitData(u).frTrain(5,:)];
        end
    end
    
    metaData{dp}.baselineFR = [mean(baselineFR) std(baselineFR)];
    metaData{dp}.amplitudeTrace = [mean(ampTrace); std(ampTrace)];
    metaData{dp}.psth = [mean(psth); std(psth)];
end

%% July 2020 Muscimol
f1 = figure;
set(f1,'Position',[200 200 1200 500]);
subplot(1,4,1);
y = [metaData{1}.responsiveUnits (metaData{1}.totalUnits-metaData{1}.responsiveUnits);...
    metaData{2}.responsiveUnits (metaData{2}.totalUnits-metaData{2}.responsiveUnits)];
bar(y,'stacked');
title('Units detected');
xticklabels({'Pre','Musc.'});
legend({'Sound resp.','Not sound resp.'});

subplot(1,4,2); hold on;
data = [metaData{1}.baselineFR(1) metaData{2}.baselineFR(1)];
bar([1 2],data);
er = errorbar([1 2],data,[metaData{1}.baselineFR(2) metaData{2}.baselineFR(2)]);
er.Color = [0 0 0];
er.LineStyle = 'none';
ylabel('Firing rate (Hz)');
title('Baseline FR');
xticks([1 2]);
xticklabels({'Pre','Musc.'});

subplot(1,4,3); hold on;
errorbar([0 10 30 50 70 90],metaData{1}.amplitudeTrace(1,:),metaData{1}.amplitudeTrace(2,:),'Color',[0 0 0]);
errorbar([0 10 30 50 70 90],metaData{2}.amplitudeTrace(1,:),metaData{2}.amplitudeTrace(2,:),'Color',[0 0 1]);
xlabel('Sound intensity (dB)');
ylabel('Firing rate (Hz)');
legend({'Pre','Musc.'});
title('Sound response');

subplot(1,4,4); hold on;
timePoints = (analysisParams.analysisWindow(1)+analysisParams.frBinWidth):1:analysisParams.analysisWindow(2);
plot(timePoints,metaData{1}.psth(1,:),'Color',[0 0 0])
plot(timePoints,metaData{2}.psth(1,:),'Color',[0 0 1]);
xlabel('Time (ms)');
ylabel('Firing rate (Hz)');
title('PSTH (90dB)');

suptitle('July 2020 muscimol');

%% Kath's Muscimol
f2 = figure;
set(f2,'Position',[200 200 1200 500]);
subplot(1,4,1);
y = [metaData{3}.responsiveUnits (metaData{3}.totalUnits-metaData{3}.responsiveUnits);...
    metaData{4}.responsiveUnits (metaData{4}.totalUnits-metaData{4}.responsiveUnits)];
bar(y,'stacked');
title('Units detected');
xticklabels({'Pre','Musc.'});
legend({'Sound resp.','Not sound resp.'});

subplot(1,4,2); hold on;
data = [metaData{3}.baselineFR(1) metaData{4}.baselineFR(1)];
bar([1 2],data);
er = errorbar([1 2],data,[metaData{3}.baselineFR(2) metaData{4}.baselineFR(2)]);
er.Color = [0 0 0];
er.LineStyle = 'none';
ylabel('Firing rate (Hz)');
title('Baseline FR');
xticks([1 2]);
xticklabels({'Pre','Musc.'});

subplot(1,4,3); hold on;
errorbar([0 10 30 50 70 90],metaData{3}.amplitudeTrace(1,:),metaData{3}.amplitudeTrace(2,:),'Color',[0 0 0]);
errorbar([0 10 30 50 70 90],metaData{4}.amplitudeTrace(1,:),metaData{4}.amplitudeTrace(2,:),'Color',[0 0 1]);
xlabel('Sound intensity (dB)');
ylabel('Firing rate (Hz)');
legend({'Pre','Musc.'});
title('Sound response');

subplot(1,4,4); hold on;
timePoints = (analysisParams.analysisWindow(1)+analysisParams.frBinWidth):1:analysisParams.analysisWindow(2);
plot(timePoints,metaData{3}.psth(1,:),'Color',[0 0 0])
plot(timePoints,metaData{4}.psth(1,:),'Color',[0 0 1]);
xlabel('Time (ms)');
ylabel('Firing rate (Hz)');
title('PSTH (90dB)');

suptitle('Kaths muscimol');

%% Saline
f3 = figure;
set(f3,'Position',[200 200 1200 500]);
subplot(1,4,1);
y = [metaData{5}.responsiveUnits (metaData{5}.totalUnits-metaData{5}.responsiveUnits);...
    metaData{6}.responsiveUnits (metaData{6}.totalUnits-metaData{6}.responsiveUnits)];
bar(y,'stacked');
title('Units detected');
xticklabels({'Pre','Sal.'});
legend({'Sound resp.','Not sound resp.'});

subplot(1,4,2); hold on;
data = [metaData{5}.baselineFR(1) metaData{6}.baselineFR(1)];
bar([1 2],data);
er = errorbar([1 2],data,[metaData{5}.baselineFR(2) metaData{6}.baselineFR(2)]);
er.Color = [0 0 0];
er.LineStyle = 'none';
ylabel('Firing rate (Hz)');
title('Baseline FR');
xticks([1 2]);
xticklabels({'Pre','Sal.'});

subplot(1,4,3); hold on;
errorbar([0 10 30 50 70 90],metaData{5}.amplitudeTrace(1,:),metaData{5}.amplitudeTrace(2,:),'Color',[0 0 0]);
errorbar([0 10 30 50 70 90],metaData{6}.amplitudeTrace(1,:),metaData{6}.amplitudeTrace(2,:),'Color',[0 0 1]);
xlabel('Sound intensity (dB)');
ylabel('Firing rate (Hz)');
legend({'Pre','Sal.'});
title('Sound response');

subplot(1,4,4); hold on;
timePoints = (analysisParams.analysisWindow(1)+analysisParams.frBinWidth):1:analysisParams.analysisWindow(2);
plot(timePoints,metaData{5}.psth(1,:),'Color',[0 0 0])
plot(timePoints,metaData{6}.psth(1,:),'Color',[0 0 1]);
xlabel('Time (ms)');
ylabel('Firing rate (Hz)');
title('PSTH (90dB)');

suptitle('Saline');