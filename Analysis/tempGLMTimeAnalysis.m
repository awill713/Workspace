clear 

saveDir = fullfile('/Volumes/AARON DATA/','Electrophysiology');

dataPaths{1} = fullfile('EP010','AW159','20201213-1');
dataPaths{2} = fullfile('EP010','AW159','20201213-2');
dataPaths{3} = fullfile('EP010','AW162','20210102-1');
dataPaths{4} = fullfile('EP010','AW162','20210102-2');
dataPaths{5} = fullfile('EP010','AW163','20210102-1');
dataPaths{6} = fullfile('EP010','AW163','20210102-2');
dataPaths{7} = fullfile('EP010','AW164','20210105');
dataPaths{8} = fullfile('EP010','AW165','20210106-1');
dataPaths{9} = fullfile('EP010','AW165','20210106-2');



for dp = 1:length(dataPaths)
    dp
%     dataFile = fullfile('D:\Electrophysiology\',dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise_Video_final','AVMultiContrastDriftingGratingsWhiteNoiseData.mat');
    dataFile = fullfile(saveDir,dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise_Video_final','AVMultiContrastDriftingGratingsWhiteNoiseData.mat');
    
    load(dataFile);
%     glmPath = dir(fullfile('D:\Electrophysiology\',dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise_GLM_TimeCourse','glmTimeCourseData_10msBin.mat'));
    glmPath = dir(fullfile(saveDir,dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise_GLM_TimeCourse','glmTimeCourseData_10msBin.mat'));
    
    load(fullfile(glmPath(end).folder,glmPath(end).name));
    
    binCount = glmAnalysisParams.binCount;
    binEdges = glmAnalysisParams.binEdges;
    totalUnits = length(unitData);
    frBinWidth = glmAnalysisParams.frBinWidth;
    
    if ~exist('lightTrain','var')
        lightTrain = zeros(0,binCount);
        soundTrain = zeros(0,binCount);
        lightSoundTrain = zeros(0,binCount);
        lightMoveTrain = zeros(0,binCount);
        soundMoveTrain = zeros(0,binCount);
        lightSoundMoveTrain = zeros(0,binCount);
        lightSoundMoveMeanTrain = zeros(0,binCount);
        nothingTrain = zeros(0,binCount);
        infooo = [];
        rSquares = [];
        rmse = [];
        theEstimates = [];
        coeffMat = [];
    end
    
    
    for n = 1:totalUnits
        neuronNumber = unitData(n).neuronNumber;
%                 if ismember(neuronNumber,responsiveUnitsGLM.lightSoundInteractUnits) ||...
%                         (ismember(neuronNumber,responsiveUnitsGLM.lightResponsiveUnits) &&...
%                         ismember(neuronNumber,responsiveUnitsGLM.soundResponsiveUnits))
        if (ismember(neuronNumber,responsiveUnitsGLM.lightMoveInteractUnits) ||...
                (ismember(neuronNumber,responsiveUnitsGLM.lightResponsiveUnits) &&...
                ismember(neuronNumber,responsiveUnitsGLM.movementResponsiveUnits))) &&...
                ~isempty(unitGLMTimeCourse(n).coefficients)
%         if ismember(neuronNumber,responsiveUnitsGLM.lightMoveInteractUnits) &&...
%                 ~isempty(unitGLMTimeCourse(n).coefficients)
            
            coeff = unitGLMTimeCourse(n).coefficients;
            coeffMat = cat(3,coeffMat,coeff);
            
            lTrain = 1000/frBinWidth * exp(coeff(1,:) + coeff(2,:)*1 + coeff(3,:)*boolean(0) + coeff(4,:)*1 + coeff(5,:)*1*boolean(0) + coeff(6,:)*1*1 + coeff(7,:)*boolean(0)*1);
            sTrain = 1000/frBinWidth * exp(coeff(1,:) + coeff(2,:)*0 + coeff(3,:)*boolean(1) + coeff(4,:)*1 + coeff(5,:)*0*boolean(1) + coeff(6,:)*0*1 + coeff(7,:)*boolean(1)*1);
            lsTrain = 1000/frBinWidth * exp(coeff(1,:) + coeff(2,:)*1 + coeff(3,:)*boolean(1) + coeff(4,:)*1 + coeff(5,:)*1*boolean(1) + coeff(6,:)*1*1 + coeff(7,:)*boolean(1)*1);
            lmTrain = 1000/frBinWidth * exp(coeff(1,:) + coeff(2,:)*1 + coeff(3,:)*boolean(0) + coeff(4,:)*4 + coeff(5,:)*1*boolean(0) + coeff(6,:)*1*4 + coeff(7,:)*boolean(0)*4);
            lsmTrain = 1000/frBinWidth * exp(coeff(1,:) + coeff(2,:)*1 + coeff(3,:)*boolean(1) + coeff(4,:)*4 + coeff(5,:)*1*boolean(1) + coeff(6,:)*1*4 + coeff(7,:)*boolean(1)*4);
            nothing = 1000/frBinWidth * exp(coeff(1,:) + coeff(2,:)*0 + coeff(3,:)*boolean(0) + coeff(4,:)*1 + coeff(5,:)*0*boolean(0) + coeff(6,:)*0*1 + coeff(7,:)*boolean(0)*1);
            smTrain = 1000/frBinWidth * exp(coeff(1,:) + coeff(2,:)*0 + coeff(3,:)*boolean(1) + coeff(4,:)*2 + coeff(5,:)*0*boolean(1) + coeff(6,:)*0*2 + coeff(7,:)*boolean(1)*2);
            lsmTrainMean = 1000/frBinWidth * exp(coeff(1,:) + coeff(2,:)*1 + coeff(3,:)*boolean(1) + coeff(4,:)*2 + coeff(5,:)*1*boolean(1) + coeff(6,:)*1*2 + coeff(7,:)*boolean(1)*2);
            
%             baselineBins = [1 find(binEdges==-1)];
%             lBase = mean(lTrain(baselineBins(1):baselineBins(2))); lSTD = std(lTrain(baselineBins(1):baselineBins(2)));
%             lsBase = mean(lsTrain(baselineBins(1):baselineBins(2))); lsSTD = std(lsTrain(baselineBins(1):baselineBins(2)));
%             lmBase = mean(lmTrain(baselineBins(1):baselineBins(2))); lmSTD = std(lmTrain(baselineBins(1):baselineBins(2)));
%             lsmBase = mean(lsmTrain(baselineBins(1):baselineBins(2))); lsmSTD = std(lsmTrain(baselineBins(1):baselineBins(2)));
%             nothingBase = mean(nothing(baselineBins(1):baselineBins(2))); nothingSTD = std(nothing(baselineBins(1):baselineBins(2)));
%             
%             lNorm = (lTrain-lBase)./lSTD;
%             lsNorm = (lsTrain-lsBase)./lsSTD;
%             lmNorm = (lmTrain-lmBase)./lmSTD;
%             lsmNorm = (lsmTrain-lsmBase)./lsmSTD;
%             nothingNorm = (nothing-nothingBase)./nothingSTD;
            
            lightTrain = [lightTrain; lTrain];
            soundTrain = [soundTrain; sTrain];
            soundMoveTrain = [soundMoveTrain; smTrain];
            lightSoundTrain = [lightSoundTrain; lsTrain];
            lightMoveTrain = [lightMoveTrain; lmTrain];
            lightSoundMoveTrain = [lightSoundMoveTrain; lsmTrain];
            lightSoundMoveMeanTrain = [lightSoundMoveMeanTrain; lsmTrainMean];
            nothingTrain = [nothingTrain; nothing];
            
            infooo = [infooo; dp n neuronNumber];
            
            psth1 = mean(unitData(n).frTrain(49:60,:));
            psth2 = mean(unitData(n).frTrain(61:72,:));
            psth3 = mean(unitData(n).frTrain(109:120,:));
            
            lTrainShort = lTrain(101:end);
            smTrainShort = smTrain(101:end);
            lsmTrainShort =lsmTrainMean(101:end);
            
            lOnset = [mean(psth1(91:301)) mean(lTrainShort(91:301))];
            smOnset = [mean(psth2(91:301)) mean(smTrainShort(91:301))];
            lsmOnset = [mean(psth3(91:301)) mean(lsmTrainShort(91:301))];
            
            if max([lOnset smOnset lsmOnset])<10000
                theEstimates = [theEstimates; lOnset smOnset lsmOnset];
            end
            
%             f1 = figure;
%             set(f1,'Position',[200 200 1200 500]);
%             subplot(1,3,1);hold on;
%             plot(analysisParams.binEdges,psth1);
%             plot(glmAnalysisParams.binEdges,lTrain);
            xc = corrcoef(psth1,lTrainShort);r2_l = xc(1,2)^2;
            rmse_l = sqrt(mean((lTrainShort-psth1).^2));
%             legend({'Light','GLM'});
%             title(['r^2 = ' num2str(r2_l)]);
%             
%             subplot(1,3,2);hold on;
%             plot(analysisParams.binEdges,psth2);
%             plot(glmAnalysisParams.binEdges,smTrain);
            xc = corrcoef(psth2,smTrainShort);r2_sm = xc(1,2)^2;
            rmse_sm = sqrt(mean((smTrainShort-psth2).^2));
%             legend({'Sound','GLM'});
%             title(['r^2 = ' num2str(r2_sm)]);
%             
%             subplot(1,3,3);hold on;
%             plot(analysisParams.binEdges,psth3);
%             plot(glmAnalysisParams.binEdges,lsmTrainMean);
            xc = corrcoef(psth3,lsmTrainShort);r2_lsm = xc(1,2)^2;
            rmse_lsm = sqrt(mean((lsmTrainShort-psth3).^2));
%             legend({'Audiovisual','GLM'});
%             title(['r^2 = ' num2str(r2_lsm)]);
%             
%             suptitle(['DP = ' num2str(dp) ', n = ' num2str(n) ', neuronNumber = ' num2str(neuronNumber)]);
%         
%             saveas(f1,fullfile('D:\Electrophysiology\',dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise_GLM_TimeCourse',['Unit ' num2str(neuronNumber) ' actual and GLM PSTH.fig']));
%             saveas(f1,fullfile('D:\Electrophysiology\',dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise_GLM_TimeCourse',['Unit ' num2str(neuronNumber) ' actual and GLM PSTH.jpg']));
%             close(f1);
            
            rSquares = [rSquares; r2_l r2_sm r2_lsm dp n neuronNumber];
            rmse = [rmse; rmse_l rmse_sm rmse_lsm dp n neuronNumber];
        end
    end
    
%     ex1 = isinf(max(lightTrain,[],2));
%     [row col] = find(isnan(lightTrain)); ex2 = unique(row);
%     ex3 = find(max(lightTrain,[],2)>10000 | max(soundTrain,[],2)>10000);
%     exclusionRows = union(union(ex1,ex2),ex3);
    
%     exclusionRows = find(isinf(max(lightTrain,[],2)) | isnan(max(lightTrain,[],2)));
end
    ex1 = find(isinf(max(lightTrain,[],2)));
    [row col] = find(isnan(lightTrain)); ex2 = unique(row);
    ex3 = find(max(lightTrain,[],2)>10000 | max(soundTrain,[],2)>10000 | max(lightSoundMoveTrain,[],2)>10000);
    exclusionRows = union(union(ex1,ex2),ex3);
% exclusionRows = find(isinf(max(lightTrain,[],2)) | isnan(max(lightTrain,[],2)));

lightTrain(exclusionRows,:) = [];
soundTrain(exclusionRows,:) = [];
lightSoundTrain(exclusionRows,:) = [];
lightMoveTrain(exclusionRows,:) = [];
soundMoveTrain(exclusionRows,:) = [];
lightSoundMoveTrain(exclusionRows,:) = [];
lightSoundMoveMeanTrain(exclusionRows,:) = [];
nothingTrain(exclusionRows,:) = [];
rmse(exclusionRows,:) = [];

f1 = figure;hold on;
nt = mean(nothingTrain);lt = mean(lightTrain);
ntSTD = std(nothingTrain,[],1)./sqrt(size(nothingTrain,1));
ltSTD = std(lightTrain,[],1)./sqrt(size(lightTrain,1));
hA = area(binEdges./1000,[nt-ntSTD; 2*ntSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.3;
hA(2).EdgeColor = [1 1 1];
hA = area(binEdges./1000,[lt-ltSTD; 2*ltSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(binEdges./1000,nt,'Color',[0 0 0],'LineStyle',':');
plot(binEdges./1000,lt,'Color',[0 0 0]);
xlabel('Time (ms)');
ylabel('GLM-estimated firing rate (Hz)');


f2 = figure;hold on;
lt = mean(lightTrain);ls = mean(lightSoundTrain);
ltSTD = std(lightTrain,[],1)./sqrt(size(lightTrain,1));
lsSTD = std(lightSoundTrain,[],1)./sqrt(size(lightSoundTrain,1));
hA = area(binEdges./1000,[lt-ltSTD; 2*ltSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
hA = area(binEdges./1000,[ls-lsSTD; 2*lsSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 1];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(binEdges./1000,lt,'Color',[0 0 0]);
plot(binEdges./1000,ls,'Color',[0 0 1]);
xlabel('Time (s)');
ylabel('GLM-estimated firing rate (Hz)');


f3 = figure;hold on;
lt = mean(lightTrain);ls = mean(lightSoundTrain);
lm = mean(lightMoveTrain); lsm = mean(lightSoundMoveTrain);
ltSTD = std(lightTrain,[],1)./sqrt(size(lightTrain,1));
lsSTD = std(lightSoundTrain,[],1)./sqrt(size(lightSoundTrain,1));
lmSTD = std(lightMoveTrain,[],1)./sqrt(size(lightMoveTrain,1));
lsmSTD = std(lightSoundMoveTrain,[],1)./sqrt(size(lightSoundMoveTrain,1));
hA = area(binEdges./1000,[lt-ltSTD; 2*ltSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
hA = area(binEdges./1000,[ls-lsSTD; 2*lsSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 1];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
hA = area(binEdges./1000,[lm-lmSTD; 2*lmSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [1 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
hA = area(binEdges./1000,[lsm-lsmSTD; 2*lsmSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [1 0 1];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(binEdges./1000,lt,'Color',[0 0 0]);
plot(binEdges./1000,ls,'Color',[0 0 1]);
plot(binEdges./1000,lm,'Color',[1 0 0]);
plot(binEdges./1000,lsm,'Color',[1 0 1]);
xlabel('Time (s)');
ylabel('GLM-estimated firing rate (Hz)');

f4 = figure;hold on;
lt = mean(lightTrain);lm = mean(lightMoveTrain);
ltSTD = std(lightTrain,[],1)./sqrt(size(lightTrain,1));
lmSTD = std(lightMoveTrain,[],1)./sqrt(size(lightMoveTrain,1));
hA = area(binEdges./1000,[lt-ltSTD; 2*ltSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
hA = area(binEdges./1000,[lm-lmSTD; 2*lmSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [1 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(binEdges./1000,lt,'Color',[0 0 0]);
plot(binEdges./1000,lm,'Color',[1 0 0]);
xlabel('Time (s)');
ylabel('GLM-estimated firing rate (Hz)');


% figure;
% subplot(1,3,1);
% histogram(rSquares(:,1));
% title('Light-evoked PSTH R^2');
% subplot(1,3,2);
% histogram(rSquares(:,2));
% title('Sound-evoked PSTH R^2');
% subplot(1,3,3);
% histogram(rSquares(:,3));
% title('Light/Sound-evoked PSTH R^2');