
%% FIRING RATE

clear all;
% close all;

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

for dp = 1:length(dataPaths)
    dp
    %Load data
    dataFile = fullfile('D:\Electrophysiology\',dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise_Video_final','AVMultiContrastDriftingGratingsWhiteNoiseData.mat');
    load(dataFile);
    stimPath = dir(fullfile('D:\Electrophysiology\',dataPaths{dp},'StimInfo','*AVmultiContrastDriftingGratingsWhiteNoise_stimInfo*'));
    load(fullfile(stimPath.folder,stimPath.name));
    
    orientations = stimInfo.orientations;
    contrasts = stimInfo.contrasts;
    binEdges = analysisParams.binEdges;
    
    if ~exist('rawFiringRate','var')
        rawFiringRate = cell(1,2);
        deltaFiringRate = [];
        normFiringRate = cell(1,2);
        psth = cell(length(contrasts),2);
        linearRatio = [];
    end
    
    for u  = 1:size(unitData,2)
        neuronNumber = unitData(u).neuronNumber;
        
        if ismember(neuronNumber,responsiveUnitsGLM.lightSoundInteractUnits) ||...
                (ismember(neuronNumber,responsiveUnitsGLM.lightResponsiveUnits) &&...
                ismember(neuronNumber,responsiveUnitsGLM.soundResponsiveUnits))
            
            tempRates = zeros(2,length(contrasts));
            for c = 1:length(contrasts)
                vIndices = (c-1)*length(orientations)+1 : c*length(orientations);
                avIndices = vIndices + length(orientations)*length(contrasts);
                
                vTemp = mean(unitData(u).meanResponse(vIndices,4));
                avTemp = mean(unitData(u).meanResponse(avIndices,4));
                tempRates(:,c) = [vTemp; avTemp];
                
                vPSTH = mean(unitData(u).frTrain(vIndices,:));
                avPSTH = mean(unitData(u).frTrain(avIndices,:));
                psth{c,1} = [psth{c,1}; vPSTH];
                psth{c,2} = [psth{c,2}; avPSTH];
            end
            
            rawFiringRate{1} = [rawFiringRate{1}; tempRates(1,:)];
            rawFiringRate{2} = [rawFiringRate{2}; tempRates(2,:)];
            deltaFiringRate = [deltaFiringRate; tempRates(2,:)-tempRates(1,:)];
            normFiringRate{1} = [normFiringRate{1}; tempRates(1,:)/tempRates(1,1)];
            normFiringRate{2} = [normFiringRate{2}; tempRates(2,:)/tempRates(1,1)];
            
            linRat = (tempRates(2,:)-tempRates(1,1)) ./ (tempRates(2,1) + tempRates(1,:) - 2*tempRates(1,1));
            linearRatio = [linearRatio; linRat];
        end
    end
end
vv = mean(rawFiringRate{1});
vvSTD = std(rawFiringRate{1},[],1)./sqrt(size(rawFiringRate{1},1));
av = mean(rawFiringRate{2});
avSTD = std(rawFiringRate{2},[],1)./sqrt(size(rawFiringRate{2},1));
figure; hold on;
hA = area(contrasts,[vv-vvSTD; 2*vvSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
hA = area(contrasts,[av-avSTD; 2*avSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 1];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(contrasts,vv,'Color',[0 0 0]);
plot(contrasts,av,'Color',[0 0 1]);
xticks(contrasts);
xlabel('Contrast');
ylabel('Mean firing rate (Hz)');

dv = mean(deltaFiringRate);
dvSTD = std(deltaFiringRate,[],1)./sqrt(size(deltaFiringRate,1));
figure;hold on;
hA = area(contrasts,[dv-dvSTD; 2*dvSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(contrasts,dv,'Color',[0 0 0]);
xticks(contrasts);
xlabel('Contrast');
ylabel('\DeltaFiring rate (Hz)');

vv = mean(psth{5,1});
vvSTD = std(psth{5,1},[],1)./sqrt(size(psth{5,1},1));
av = mean(psth{5,2});
avSTD = std(psth{5,2},[],1)./sqrt(size(psth{5,2},1));
figure; hold on;
hA = area(binEdges,[vv-vvSTD; 2*vvSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
hA = area(binEdges,[av-avSTD; 2*avSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 1];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(binEdges,vv,'Color',[0 0 0]);
plot(binEdges,av,'Color',[0 0 1]);
xlabel('Time (ms)');
ylabel('Mean firing rate (Hz)');


expected = rawFiringRate{2}(:,1) + rawFiringRate{1}(:,5) - 2*rawFiringRate{1}(:,1);
actual = rawFiringRate{2}(:,5) - rawFiringRate{1}(:,1);
figure;hold on;
scatter(expected,actual);
line([-20 100],[-20 100]);
xlabel('\DeltaFR_{light} + \DeltaFR_{sound} (Hz)');
ylabel('\DeltaFR_{light+sound} (Hz)');

linearRatio([11 25 38 237 327 415 477 494 512 520],:) = [];
lr = mean(linearRatio);
lrSTD = std(linearRatio)./sqrt(size(linearRatio,1));
figure;hold on;
hA = area(contrasts,[lr-lrSTD; 2*lrSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(contrasts,lr,'Color',[0 0 0]);
xticks(contrasts);
xlabel('Contrast');
ylabel('Linearity ratio');
ylim([0.6 1.5]);



%% TIMING AND CoV

clear all;
% close all;

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

for dp = 1:length(dataPaths)
    dp
    %Load data
    dataFile = fullfile('D:\Electrophysiology\',dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise_Video_final','AVMultiContrastDriftingGratingsWhiteNoiseData.mat');
    load(dataFile);
    stimPath = dir(fullfile('D:\Electrophysiology\',dataPaths{dp},'StimInfo','*AVmultiContrastDriftingGratingsWhiteNoise_stimInfo*'));
    load(fullfile(stimPath.folder,stimPath.name));
    
    orientations = stimInfo.orientations;
    contrasts = stimInfo.contrasts;
    
    binEdges = analysisParams.binEdges;
    baselineLastBin = find(binEdges==analysisParams.baselineWindow(2));
    quantFirstBin = find(binEdges==analysisParams.quantWindow(1));
    quantLastBin = find(binEdges==analysisParams.quantWindow(2));
    %     quantLastBin = find(binEdges==200);
    
    if ~exist('latency','var')
        latency = cell(1,2);
        deltaLatency = [];
        peakLatency = cell(1,2);
        responseSlope = cell(1,2);
        deltaSlope = [];
        fwhm = cell(1,2);
        info = [];
    end
    
    for u  = 1:size(unitData,2)
        neuronNumber = unitData(u).neuronNumber;
        
        if ismember(neuronNumber,responsiveUnitsGLM.lightSoundInteractUnits) ||...
                (ismember(neuronNumber,responsiveUnitsGLM.lightResponsiveUnits) &&...
                ismember(neuronNumber,responsiveUnitsGLM.soundResponsiveUnits))
            
            tempLatency = zeros(2,length(contrasts)-1);
            tempPeak = zeros(2,length(contrasts)-1);
            tempSlope = zeros(2,length(contrasts)-1);
            tempFWHM = zeros(2,length(contrasts)-1);
            include = true;
            for c = 2:length(contrasts)
                vIndices = (c-1)*length(orientations)+1 : c*length(orientations);
                avIndices = vIndices + length(orientations)*length(contrasts);
                
                vTemp = smooth(mean(unitData(u).frTrain(vIndices,:)));
                avTemp = smooth(mean(unitData(u).frTrain(avIndices,:)));
                
                vBaseline = mean(vTemp(1:baselineLastBin)); vBaseSTD = std(vTemp(1:baselineLastBin));
                avBaseline = mean(avTemp(1:baselineLastBin)); avBaseSTD = std(avTemp(1:baselineLastBin));
                vRespWindow = vTemp(quantFirstBin:quantLastBin);
                avRespWindow = avTemp(quantFirstBin:quantLastBin);
                
                vSign = sign(mean(vTemp(quantFirstBin:quantLastBin))-mean(vTemp(1:quantFirstBin)));
                avSign = sign(mean(avTemp(quantFirstBin:quantLastBin))-mean(avTemp(1:quantFirstBin)));
                
                vBaseline = vBaseline*vSign; avBaseline = avBaseline*avSign;vTemp = vTemp*vSign;
                vRespWindow = vRespWindow*vSign; avRespWindow = avRespWindow*avSign;avTemp = avTemp*avSign;
                
                [vMaxFR vMaxBin] = max(vRespWindow);
                [avMaxFR avMaxBin] = max(avRespWindow);
                vSlope = (vMaxFR-vBaseline)/vMaxBin;
                avSlope = (avMaxFR-avBaseline)/avMaxBin;
                tempSlope(:,c-1) = [vSlope; avSlope];
                tempPeak(:,c-1) = [vMaxBin; avMaxBin];
                
                vLat = find(vRespWindow(1:vMaxBin) < (vBaseline+1*vBaseSTD),1,'last');
                avLat = find(avRespWindow(1:avMaxBin) < (avBaseline+1*avBaseSTD),1,'last');
                
                vLat = find(vRespWindow > (vBaseline+1*vBaseSTD),1,'first');
                avLat = find(avRespWindow > (avBaseline+1*avBaseSTD),1,'first');
                tempLatency(:,c-1) = [vLat; avLat];
                
                vThresh = vBaseline + 0.5*(vMaxFR-vBaseline);
                avThresh = avBaseline + 0.5*(avMaxFR-avBaseline);
                vFirst = find(fliplr(vRespWindow(1:vMaxBin)') < vThresh,1,'first');
                vLast = find(vTemp((vMaxBin+quantFirstBin):end) < vThresh,1,'first');
                vDur = vFirst + vLast;
                avFirst = find(fliplr(avRespWindow(1:avMaxBin)') < avThresh,1,'first');
                avLast = find(avTemp((quantFirstBin+avMaxBin):end) < avThresh,1,'first');
                avDur = avFirst + avLast;
                
                if length(vDur)==1 && length(avDur)==1 && vSign==1 && avSign==1
                    tempFWHM(:,c-1) = [vDur; avDur];
                else
                    include = false;
                end
            end
            
            responseSlope{1} = [responseSlope{1}; tempSlope(1,:)];
            responseSlope{2} = [responseSlope{2}; tempSlope(2,:)];
            deltaSlope = [deltaSlope; tempSlope(2,:)-tempSlope(1,:)];
            
            peakLatency{1} = [peakLatency{1}; tempPeak(1,:)];
            peakLatency{2} = [peakLatency{2}; tempPeak(2,:)];
            
            if include
                fwhm{1} = [fwhm{1}; tempFWHM(1,:)];
                fwhm{2} = [fwhm{2}; tempFWHM(2,:)];
                
                info = [info; dp u neuronNumber];
            end
            
            latency{1} = [latency{1}; tempLatency(1,:)];
            latency{2} = [latency{2}; tempLatency(2,:)];
            deltaLatency = [deltaLatency; tempLatency(2,:)-tempLatency(1,:)];
            
        end
    end
end

vv = mean(responseSlope{1});
vvSTD = std(responseSlope{1},[],1)./sqrt(size(responseSlope{1},1));
av = mean(responseSlope{2});
avSTD = std(responseSlope{2},[],1)./sqrt(size(responseSlope{2},1));
figure; hold on;
hA = area(contrasts(2:end),[vv-vvSTD; 2*vvSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
hA = area(contrasts(2:end),[av-avSTD; 2*avSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 1];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(contrasts(2:end),vv,'Color',[0 0 0]);
plot(contrasts(2:end),av,'Color',[0 0 1]);
% xticks(contrasts);
xlabel('Contrast');
ylabel('Response slope (Hz/ms)');


dv = mean(deltaSlope);
dvSTD = std(deltaSlope,[],1)./sqrt(size(deltaSlope,1));
figure;hold on;
hA = area(contrasts(2:end),[dv-dvSTD; 2*dvSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(contrasts(2:end),dv,'Color',[0 0 0]);
xlabel('Contrast');
ylabel('\Delta Response slope (Hz/ms)');

vv = mean(latency{1});
vvSTD = std(latency{1},[],1)./sqrt(size(latency{1},1));
av = mean(latency{2});
avSTD = std(latency{2},[],1)./sqrt(size(latency{2},1));
figure; hold on;
hA = area(contrasts(2:end),[vv-vvSTD; 2*vvSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
hA = area(contrasts(2:end),[av-avSTD; 2*avSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 1];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(contrasts(2:end),vv,'Color',[0 0 0]);
plot(contrasts(2:end),av,'Color',[0 0 1]);
xticks(contrasts);
xlabel('Contrast');
ylabel('Response latency (ms)');

dv = mean(deltaLatency);
dvSTD = std(deltaLatency,[],1)./sqrt(size(deltaLatency,1));
figure;hold on;
hA = area(contrasts(2:end),[dv-dvSTD; 2*dvSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(contrasts(2:end),dv,'Color',[0 0 0]);
xticks(contrasts);
xlabel('Contrast');
ylabel('\Delta response latency (ms)');

vv = mean(peakLatency{1});
vvSTD = std(peakLatency{1},[],1)./sqrt(size(peakLatency{1},1));
av = mean(peakLatency{2});
avSTD = std(peakLatency{2},[],1)./sqrt(size(peakLatency{2},1));
figure; hold on;
hA = area(contrasts(2:end),[vv-vvSTD; 2*vvSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
hA = area(contrasts(2:end),[av-avSTD; 2*avSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 1];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(contrasts(2:end),vv,'Color',[0 0 0]);
plot(contrasts(2:end),av,'Color',[0 0 1]);
xticks(contrasts);
xlabel('Contrast');
ylabel('Latency to peak response(ms)');


vv = mean(fwhm{1});
vvSTD = std(fwhm{1},[],1)./sqrt(size(fwhm{1},1));
av = mean(fwhm{2});
avSTD = std(fwhm{2},[],1)./sqrt(size(fwhm{2},1));
figure; hold on;
hA = area(contrasts(2:end),[vv-vvSTD; 2*vvSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
hA = area(contrasts(2:end),[av-avSTD; 2*avSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 1];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(contrasts(2:end),vv,'Color',[0 0 0]);
plot(contrasts(2:end),av,'Color',[0 0 1]);
xticks(contrasts);
xlabel('Contrast');
ylabel('Response FWHM (ms)');

dv = mean(fwhm{2}-fwhm{1});
dvSTD = std(fwhm{2}-fwhm{1},[],1)./sqrt(size(fwhm{2}-fwhm{1},1));
figure;hold on;
hA = area(contrasts(2:end),[dv-dvSTD; 2*dvSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(contrasts(2:end),dv,'Color',[0 0 0]);
xticks(contrasts);
xlabel('Contrast');
ylabel('\Delta FWHM (ms)');


%% Orientation shifts
clear all;
% close all;

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

orientationShuffles = 1000;
for dp = 1:length(dataPaths)
    dp
    %Load data
    dataFile = fullfile('D:\Electrophysiology\',dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise_Video_final','AVMultiContrastDriftingGratingsWhiteNoiseData.mat');
    load(dataFile);
    stimPath = dir(fullfile('D:\Electrophysiology\',dataPaths{dp},'StimInfo','*AVmultiContrastDriftingGratingsWhiteNoise_stimInfo*'));
    load(fullfile(stimPath.folder,stimPath.name));
    
    orientations = stimInfo.orientations;
    contrasts = stimInfo.contrasts;
    repeats = stimInfo.repeats;
    
    if ~exist('orientShifts','var')
        orientShifts = zeros(0,length(contrasts));
        unitOSI = cell(1,length(contrasts));
        unitDSI = cell(1,length(contrasts));
        simOrientShiffs = [];
    end
    
    for u  = 1:size(unitData,2)
        neuronNumber = unitData(u).neuronNumber;
        
        if ismember(neuronNumber,responsiveUnitsGLM.lightSoundInteractUnits) ||...
                (ismember(neuronNumber,responsiveUnitsGLM.lightResponsiveUnits) &&...
                ismember(neuronNumber,responsiveUnitsGLM.soundResponsiveUnits)) &&...
                (ismember(neuronNumber,responsiveUnits.orientationSelectiveUnits) ||...
                ismember(neuronNumber,responsiveUnits.directionSelectiveUnits))
            
            tempShifts = zeros(1,length(contrasts));
            for c = 1:length(contrasts)
                vIndices = (c-1)*length(orientations)+1 : c*length(orientations);
                avIndices = vIndices + length(orientations)*length(contrasts);
                
                vResp = unitData(u).meanResponse(vIndices,4);
                avResp = unitData(u).meanResponse(avIndices,4);
                
                [valV indV] = max(vResp);
                [valAV indAV] = max(avResp);
                orientDiff = wrapTo180(orientations(indAV) - orientations(indV));
                if orientDiff == 180
                    orientDiff = -180;
                end
                tempShifts(c) = orientDiff;
                
                
                oppIndV = mod(indV+6-1,length(orientations))+1;
                orthIndV = [mod(indV+3-1,length(orientations)) mod(indV-3-1,length(orientations))]+1;
                vOSI = (vResp(indV) - mean(vResp(orthIndV))) / (vResp(indV) + mean(vResp(orthIndV)));
                vDSI = (vResp(indV) - vResp(oppIndV)) / (vResp(indV) + vResp(oppIndV));
                
                oppIndAV = mod(indAV+6-1,length(orientations))+1;
                orthIndAV = [mod(indAV+3-1,length(orientations)) mod(indAV-3-1,length(orientations))]+1;
                avOSI = (avResp(indAV) - mean(avResp(orthIndAV))) / (avResp(indAV) + mean(avResp(orthIndAV)));
                avDSI = (avResp(indAV) - avResp(oppIndAV)) / (avResp(indAV) + avResp(oppIndV));
                if ismember(neuronNumber,responsiveUnits.orientationSelectiveUnits)
                    unitOSI{c} = [unitOSI{c}; vOSI avOSI dp u neuronNumber];
                end
                if ismember(neuronNumber,responsiveUnits.directionSelectiveUnits)
                    unitDSI{c} = [unitDSI{c}; vDSI avDSI dp u neuronNumber];
                end
                
                if c == length(contrasts)
                    
                    shuffledDiffs = zeros(1,orientationShuffles);
                    theVaverages = unitData(u).meanResponse(vIndices,4);
                    theVstds = unitData(u).meanResponse(vIndices,5)*sqrt(repeats);
                    
                    indShift = indV-indAV;
                    theAVaverages = circshift(unitData(u).meanResponse(avIndices,4),indShift);
                    theAVstds = circshift(unitData(u).meanResponse(avIndices,5)*sqrt(repeats),indShift);
                    for shuff = 1:orientationShuffles
                        vShuffMeans = zeros(1,length(orientations));
                        avShuffMeans = zeros(1,length(orientations));
                        for oo = 1:length(orientations)
                            poissLambda = vResp(oo)./analysisParams.quantScalar;
                            randPoissTrials = poissrnd(poissLambda,[1 repeats]);
                            vShuffMeans(oo) = mean(randPoissTrials);
                            
                            poissLambda = avResp(oo)./analysisParams.quantScalar;
                            randPoissTrials = poissrnd(poissLambda,[1 repeats]);
                            avShuffMeans(oo) = mean(randPoissTrials);
                            
                            dirVal = mean(normrnd(theVaverages(oo),theVstds(oo),1,repeats));
%                             shuffMeans(oo) = dirVal;
                        end
                        [~, vShuffPeakInd] = max(vShuffMeans);
                        [~, avShuffPeakInd] = max(avShuffMeans);
                        shuffDiff = wrapTo180(orientations(avShuffPeakInd) - orientations(vShuffPeakInd));
                        if shuffDiff == 180
                            shuffDiff = -180;
                        end
                        shuffledDiffs(shuff) = shuffDiff;
                    end
                    simOrientShiffs = [simOrientShiffs; shuffledDiffs];
                end
            end
            orientShifts = [orientShifts; tempShifts];
        end
    end
end

figure;
for c = 1:length(contrasts)
    subplot(1,5,c);
    histogram(orientShifts(:,c),-180:30:180,'Normalization','probability')
end
figure;
histogram(orientShifts(:,length(contrasts)),-180:30:180,'Normalization','probability')

shuffPercent = zeros(orientationShuffles,length(orientations)+1);
shuffVector = zeros(1,orientationShuffles);
for sh = 1:orientationShuffles
    percent = histcounts(simOrientShiffs(:,sh),(-180-orientations(2)/2):orientations(2):(180+orientations(2)/2),'Normalization','probability');
    percent(end) = percent(1);
    shuffPercent(sh,:) = percent;
    meanVector = mean(abs(simOrientShiffs(:,sh)));
    shuffVector(sh) = meanVector;
end
percentMean = mean(shuffPercent);
percentStd = std(shuffPercent);
actualPercent = histcounts(orientShifts(:,length(contrasts)),(-180-orientations(2)/2):orientations(2):(180+orientations(2)/2),'Normalization','probability');
actualPercent(end) = actualPercent(1);
actualVector = mean(abs(orientShifts(:,length(contrasts))));
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
hx = histogram(shuffVector,'Normalization','probability');
line([actualVector actualVector],[0 0.16]);
xlabel('Mean \Delta orientation');
ylabel('Randomization instances');

figure;
for c = 1:length(contrasts)
    subplot(1,length(contrasts),c);hold on;
    scatter(unitOSI{c}(:,1),unitOSI{c}(:,2));
    line([0 1],[0 1]);
    [h p] = ttest(unitOSI{c}(:,1),unitOSI{c}(:,2));
    title(['p = ' num2str(p)]);
    
    meanOSI(1,c) = mean(unitOSI{c}(:,1));
    meanOSI(2,c) = mean(unitOSI{c}(:,2));
    meanOSI(3,c) = mean(unitOSI{c}(:,2))-mean(unitOSI{c}(:,1));
    stdOSI(1,c) = std(unitOSI{c}(:,1)) ./ sqrt(length(unitOSI{c}(:,1)));
    stdOSI(2,c) = std(unitOSI{c}(:,2)) ./ sqrt(length(unitOSI{c}(:,2)));
    stdOSI(3,c) = std(unitOSI{c}(:,2))-mean(unitOSI{c}(:,1)) ./ sqrt(length(unitOSI{c}(:,1)));
end
suptitle('OSI');


figure; hold on;
hA = area(contrasts,[meanOSI(1,:)-stdOSI(1,:); 2*stdOSI(1,:)]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
hA = area(contrasts,[meanOSI(2,:)-stdOSI(2,:); 2*stdOSI(2,:)]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 1];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(contrasts,meanOSI(1,:),'Color',[0 0 0]);
plot(contrasts,meanOSI(2,:),'Color',[0 0 1]);
xticks(contrasts);
xlabel('Contrast');
ylabel('Orientation selectivity index');


figure;hold on;
hA = area(contrasts,[meanOSI(3,:)-stdOSI(3,:); 2*stdOSI(3,:)]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(contrasts,meanOSI(3,:),'Color',[0 0 0]);
xticks(contrasts);
xlabel('Contrast');
ylabel('\Delta Orientation selectivity index');

figure;
for c = 1:length(contrasts)
    subplot(1,length(contrasts),c);hold on;
    scatter(unitDSI{c}(:,1),unitDSI{c}(:,2));
    line([0 1],[0 1]);
    [h p] = ttest(unitDSI{c}(:,1),unitDSI{c}(:,2));
    title(['p = ' num2str(p)]);
    
    meanDSI(1,c) = mean(unitDSI{c}(:,1));
    meanDSI(2,c) = mean(unitDSI{c}(:,2));
    meanDSI(3,c) = mean(unitDSI{c}(:,2))-mean(unitDSI{c}(:,1));
    stdDSI(1,c) = std(unitDSI{c}(:,1)) ./ sqrt(length(unitDSI{c}(:,1)));
    stdDSI(2,c) = std(unitDSI{c}(:,2)) ./ sqrt(length(unitDSI{c}(:,2)));
    stdDSI(3,c) = std(unitDSI{c}(:,2))-mean(unitDSI{c}(:,1)) ./ sqrt(length(unitDSI{c}(:,1)));
end
suptitle('DSI');


figure; hold on;
hA = area(contrasts,[meanDSI(1,:)-stdDSI(1,:); 2*stdDSI(1,:)]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
hA = area(contrasts,[meanDSI(2,:)-stdDSI(2,:); 2*stdDSI(2,:)]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 1];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(contrasts,meanDSI(1,:),'Color',[0 0 0]);
plot(contrasts,meanDSI(2,:),'Color',[0 0 1]);
xticks(contrasts);
xlabel('Contrast');
ylabel('Direction selectivity index');


figure;hold on;
hA = area(contrasts,[meanDSI(3,:)-stdDSI(3,:); 2*stdDSI(3,:)]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(contrasts,meanDSI(3,:),'Color',[0 0 0]);
xticks(contrasts);
xlabel('Contrast');
ylabel('\Delta Direction selectivity index');
            
%% Coefficient of variation or Fano Factor

clear all;
% close all;

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

for dp = 1:length(dataPaths)
    dp
    %Load data
    dataFile = fullfile('D:\Electrophysiology\',dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise_Video_final','AVMultiContrastDriftingGratingsWhiteNoiseData.mat');
    load(dataFile);
    stimPath = dir(fullfile('D:\Electrophysiology\',dataPaths{dp},'StimInfo','*AVmultiContrastDriftingGratingsWhiteNoise_stimInfo*'));
    load(fullfile(stimPath.folder,stimPath.name));
    
    orientations = stimInfo.orientations;
    contrasts = stimInfo.contrasts;
    uniqueEvents = length(stimInfo.index);
    
    %     quantLastBin = find(binEdges==200);
    
    if ~exist('coefVar','var')
        cvData = cell(1,length(contrasts)); 
    end
    
    for u  = 1:size(unitData,2)
        neuronNumber = unitData(u).neuronNumber;
        
        if ismember(neuronNumber,responsiveUnitsGLM.lightSoundInteractUnits) ||...
                (ismember(neuronNumber,responsiveUnitsGLM.lightResponsiveUnits) &&...
                ismember(neuronNumber,responsiveUnitsGLM.soundResponsiveUnits))
                
            for c = 1:length(contrasts)
                vIndices = (c-1)*length(orientations)+1 : c*length(orientations);
                avIndices = vIndices + length(orientations)*length(contrasts);
                
                vMeans = mean(unitData(u).trialResponse(vIndices,:),2);
                vSTDs = std(unitData(u).trialResponse(vIndices,:),[],2);
%                 vMeans = mean(reshape(unitData(u).trialResponse(vIndices,:),[1 120]));
%                 vSTDs = std(reshape(unitData(u).trialResponse(vIndices,:),[1 120]));
                vCV = nanmean(vSTDs./vMeans);
                vFF = nanmean((vSTDs.^2)./vMeans);
                
                avMeans = mean(unitData(u).trialResponse(avIndices,:),2);
                avSTDs = std(unitData(u).trialResponse(avIndices,:),[],2);
%                 avMeans = mean(reshape(unitData(u).trialResponse(avIndices,:),[1 120]));
%                 avSTDs = std(reshape(unitData(u).trialResponse(avIndices,:),[1 120]));
                avCV = nanmean(avSTDs./avMeans);
                avFF = nanmean((avSTDs.^2)./avMeans);
                
                cvData{c} = [cvData{c}; vCV vFF avCV avFF dp u neuronNumber];
            end
                
                data = zeros(uniqueEvents,5); %mean, std, var, cv, fano
                for e = 1:uniqueEvents
                    avg = mean(unitData(u).trialResponse(e,:)/analysisParams.quantScalar);
                    stdd = std(unitData(u).trialResponse(e,:)/analysisParams.quantScalar);
                    
%                     avg = mean(unitData(u).trialResponse(e,:));
%                     stdd = std(unitData(u).trialResponse(e,:));
                    
                    data(e,1) = avg;
                    data(e,2) = stdd;
                    data(e,3) = stdd^2;
                    data(e,4) = stdd / avg;
                    data(e,5) = stdd^2 / avg;
                end
                figure;hold on;
                scatter(data(1:60,1),data(1:60,4),'MarkerFaceColor',[0 0 0]);
                scatter(data(61:120,1),data(61:120,4),'MarkerFaceColor',[0 0 1]);
                xs = min(data(:,1)):.1:max(data(:,1));
                ys = sqrt(nanmean(data(:,5)))./sqrt(xs);
%                 ys = 1./sqrt(xs);
                plot(xs,ys);
        end
    end
end

for c = 1:length(contrasts)
    vCVMean(c) = mean(cvData{c}(:,1));
    vCVSTD(c) = std(cvData{c}(:,1))./sqrt(length(cvData{c}(:,1)));
    avCVMean(c) = mean(cvData{c}(:,3));
    avCVSTD(c) = std(cvData{c}(:,3))./sqrt(length(cvData{c}(:,3)));
    dCVMean(c) = mean(cvData{c}(:,3) - cvData{c}(:,1));
    dCVSTD(c) = std(cvData{c}(:,3) - cvData{c}(:,1))./sqrt(length(cvData{c}(:,1)));
    
    vFFMean(c) = mean(cvData{c}(:,2));
    vFFSTD(c) = std(cvData{c}(:,2))./sqrt(length(cvData{c}(:,2)));
    avFFMean(c) = mean(cvData{c}(:,4));
    avFFSTD(c) = std(cvData{c}(:,4))./sqrt(length(cvData{c}(:,4)));
    dFFMean(c) = mean(cvData{c}(:,4) - cvData{c}(:,2));
    dFFSTD(c) = std(cvData{c}(:,4) - cvData{c}(:,2))./sqrt(length(cvData{c}(:,1)));
end
figure; hold on;
hA = area(contrasts,[vCVMean-vCVSTD; 2*vCVSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
hA = area(contrasts,[avCVMean-avCVSTD; 2*avCVSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 1];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(contrasts,vCVMean,'Color',[0 0 0]);
plot(contrasts,avCVMean,'Color',[0 0 1]);
xticks(contrasts);
xlabel('Contrast');
ylabel('Coefficient of variation');


figure;hold on;
hA = area(contrasts,[dCVMean-dCVSTD; 2*dCVSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(contrasts,dCVMean,'Color',[0 0 0]);
xticks(contrasts);
xlabel('Contrast');
ylabel('\Delta Coefficient of variation');


figure; hold on;
hA = area(contrasts,[vFFMean-vFFSTD; 2*vFFSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
hA = area(contrasts,[avFFMean-avFFSTD; 2*avFFSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 1];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(contrasts,vFFMean,'Color',[0 0 0]);
plot(contrasts,avFFMean,'Color',[0 0 1]);
xticks(contrasts);
xlabel('Contrast');
ylabel('Fano factor (Hz)');


figure;hold on;
hA = area(contrasts,[dFFMean-dFFSTD; 2*dFFSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(contrasts,dFFMean,'Color',[0 0 0]);
xticks(contrasts);
xlabel('Contrast');
ylabel('\Delta Fano factor (Hz)');

