
clear all;
% close all;

saveDir = fullfile('D:\Electrophysiology','EP008');
if ~exist(saveDir)
    mkdir(saveDir);
end

dataPaths{1} = fullfile('EP008','AW143','20200919');
dataPaths{2} = fullfile('EP008','AW144','20200919');
dataPaths{3} = fullfile('EP008','AW145','20200919');

depthOffsets(1) = 780;
depthOffsets(2) = 900;
depthOffsets(3) = 775;

laserResponsive = [];
singleSynapse = [];
allUnits = [];


for dp = 1:length(dataPaths)
    dp
    %Load data
    dataFile = fullfile('D:\Electrophysiology\',dataPaths{dp},'OptoLaserResponses','OptoResponseData.mat');
    load(dataFile);
    
    depthFile = fullfile('D:\Electrophysiology\',dataPaths{dp},'unitDepths.mat');
    load(depthFile);
    
    belowLat = responsiveUnits.latencies(find(responsiveUnits.latencies(:,2)<8),1);
    
    opto = -clusterChannelAndDepth(2,responsiveUnits.optoResponsiveUnits);
    single = -clusterChannelAndDepth(2,belowLat);
    all = -clusterChannelAndDepth(2,union(responsiveUnits.singleUnits,responsiveUnits.multiUnits));
    
%     opto = opto + depthOffsets(dp)-775;
%     all = all + depthOffsets(dp)-775;
%     single = single + depthOffsets(dp)-775;
    
    laserResponsive = [laserResponsive opto];
    singleSynapse = [singleSynapse single];
    allUnits = [allUnits all];
end

figure;hold on;
histogram(laserResponsive,0:50:1000,'Normalization','probability');
xlabel('Depth (um)');
ylabel('Number of units');
title('Laser responsive units');

figure;hold on;
histogram(singleSynapse,0:50:1000,'Normalization','probability');
xlabel('Depth (um)');
ylabel('Number of units');
title('Monosynaptic laser responsive units');


figure;hold on;
histogram(allUnits,0:50:1000,'Normalization','probability');
xlabel('Depth (um)');
ylabel('Number of units');
title('All units');
