
clear all;
% close all;

saveDir = fullfile('D:\Electrophysiology','EP004','MetaData - AVMultiContrastDriftingGratingsWhiteNoise','MLE model','Individual neurons - recat');
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
% dataPaths{8} = fullfile('EP004','AW124','20200303-2');

for dp = 1:length(dataPaths)
    
    %Load data
    dataFile = fullfile('D:\Electrophysiology\',dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise_Recat','AVMultiContrastDriftingGratingsWhiteNoiseData.mat');
    load(dataFile);
    
    orientations = linspace(0,330,12);
    directionSelectiveUnits = [];
    orientationSelectiveUnits = [];
    
    for u = 1:size(unitData,2)
        neuronNumber = unitData(u).neuronNumber;
        
        if ismember(neuronNumber,responsiveUnits.lightResponsiveUnits)
            
            meanResponse = unitData(u).meanResponse;
            trialResponse = unitData(u).trialResponse;
            
            c=5;
            baseInd = (c-1)*length(orientations);
            visMeanResp = meanResponse(baseInd+1:baseInd+length(orientations),4);
            [~, maxInd] = max(visMeanResp);
            oppInd = mod(maxInd+6-1,length(orientations))+1;
            orthInd = [mod(maxInd+3-1,length(orientations)) mod(maxInd-3-1,length(orientations))]+1;
            
            
            [~, pDir] = ttest2(trialResponse(baseInd+maxInd,:),trialResponse(baseInd+oppInd,:));
            if pDir<0.05
                directionSelectiveUnits = [directionSelectiveUnits neuronNumber];
            end
            
            [~, pOrient] = ttest2(trialResponse(baseInd+maxInd,:),[trialResponse(baseInd+orthInd(1),:) trialResponse(baseInd+orthInd(2),:)]);
            if pOrient<0.05
                orientationSelectiveUnits = [orientationSelectiveUnits neuronNumber];
            end
            
        end
    end
    
    responsiveUnits2.directionSelectiveUnits = directionSelectiveUnits;
    responsiveUnits2.orientationSelectiveUnits = orientationSelectiveUnits;
    save(fullfile('D:\Electrophysiology\',dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise_Recat','AVMultiContrastDriftingGratingsWhiteNoiseData.mat'),'responsiveUnits2','-append');

end

    
    