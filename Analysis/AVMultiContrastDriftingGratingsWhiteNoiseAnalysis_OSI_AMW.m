

clear all;
% close all;

saveDir = fullfile('D:\Electrophysiology','EP010','MetaData - AVMultiContrastDriftingGratingsWhiteNoise');
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
dataPaths{9} = fullfile('EP010','AW157','20201212-1');
dataPaths{10} = fullfile('EP010','AW157','20201212-2');
dataPaths{11} = fullfile('EP010','AW158','20201212-1');
dataPaths{12} = fullfile('EP010','AW158','20201212-2');

% dataPaths{13} = fullfile('EP010','AW159','20201213-1');
% dataPaths{14} = fullfile('EP010','AW159','20201213-2');
% dataPaths{15} = fullfile('EP010','AW162','20210102-1');
% dataPaths{16} = fullfile('EP010','AW162','20210102-2');
% dataPaths{17} = fullfile('EP010','AW163','20210102-1');
% dataPaths{18} = fullfile('EP010','AW163','20210102-2');
% dataPaths{19} = fullfile('EP010','AW164','20210105');
% dataPaths{20} = fullfile('EP010','AW165','20210106-1');
% dataPaths{21} = fullfile('EP010','AW165','20210106-2');

osiData = [];
dsiData = [];

for dp = 1:length(dataPaths)
    dp
    %Load data
    dataFile = fullfile('D:\Electrophysiology\',dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise_final','AVMultiContrastDriftingGratingsWhiteNoiseData_selectivity.mat');
    load(dataFile);
    
    for u = 1:length(unitData)
        
        neuronNumber = unitData(u).neuronNumber;
        
        if (ismember(neuronNumber,responsiveUnitsGLM.lightSoundInteractUnits) ||...
        (ismember(neuronNumber,responsiveUnitsGLM.soundResponsiveUnits) &&...
        ismember(neuronNumber,responsiveUnitsGLM.lightResponsiveUnits))) &&...
        ismember(neuronNumber,responsiveUnits.orientationSelectiveUnits) &&...
        ismember(neuronNumber,[responsiveUnits.singleUnits responsiveUnits.multiUnits])
    
            v = unitData(u).OSI(1);
            av = unitData(u).OSI(2);
            
            osiData = [osiData; v av dp u neuronNumber];
        end
        if (ismember(neuronNumber,responsiveUnitsGLM.lightSoundInteractUnits) ||...
        (ismember(neuronNumber,responsiveUnitsGLM.soundResponsiveUnits) &&...
        ismember(neuronNumber,responsiveUnitsGLM.lightResponsiveUnits))) &&...
        ismember(neuronNumber,responsiveUnits.directionSelectiveUnits) &&...
        ismember(neuronNumber,[responsiveUnits.singleUnits responsiveUnits.multiUnits])
            v = unitData(u).DSI(1);
            av = unitData(u).DSI(2);
            
            dsiData = [dsiData; v av dp u neuronNumber];
        end
    end
end
    