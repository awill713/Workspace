
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

depthOffsets(1) = 1000;
depthOffsets(2) = 700;
depthOffsets(3) = 700;
depthOffsets(4) = 1100;
depthOffsets(5) = 650;
depthOffsets(6) = 550;
depthOffsets(7) = 850;
depthOffsets(8) = 850;

soundResponsive = [];
lightResponsive = [];
orientationSelective = [];
allUnits = [];


for dp = 1:length(dataPaths)
    dp
    %Load data
    dataFile = fullfile('E:\Electrophysiology\',dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise','AVMultiContrastDriftingGratingsWhiteNoiseData.mat');
    load(dataFile);
    
    depthFile = fullfile('E:\Electrophysiology\',dataPaths{dp},'unitDepths.mat');
    load(depthFile);
    
%     light =clusterChannelAndDepth(2,intersect(responsiveUnits.lightResponsiveUnits,responsiveUnits.singleUnits));
%     sound = clusterChannelAndDepth(2,intersect(responsiveUnits.soundResponsiveUnits,responsiveUnits.singleUnits));
%     orientation = clusterChannelAndDepth(2,intersect(responsiveUnits.orientationSelectiveUnits,responsiveUnits.singleUnits));
    light = -clusterChannelAndDepth(2,responsiveUnits.lightResponsiveUnits);
    sound = -clusterChannelAndDepth(2,responsiveUnits.soundResponsiveUnits);
    orientation = -clusterChannelAndDepth(2,responsiveUnits.orientationSelectiveUnits);
    all = -clusterChannelAndDepth(2,union(responsiveUnits.singleUnits,responsiveUnits.multiUnits));
    
    lightResponsive = [lightResponsive light];
    soundResponsive = [soundResponsive sound];
    orientationSelective = [orientationSelective orientation];
    allUnits = [allUnits all];
end

figure;hold on;
histogram(lightResponsive);
xlabel('Depth (um)');
ylabel('Number of units');

figure;
histogram(soundResponsive);
xlabel('Depth (um)');
ylabel('Number of units');

figure;
histogram(orientationSelective);
xlabel('Depth (um)');
ylabel('Number of units');

figure;hold on;
histogram(allUnits);
xlabel('Depth (um)');
ylabel('Number of units');