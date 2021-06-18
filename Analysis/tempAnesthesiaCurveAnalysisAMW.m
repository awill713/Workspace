

% dataPaths{1} = fullfile('EP013','AW182','20210428');
% dataPaths{2} = fullfile('EP013','AW183','20210428');
% dataPaths{3} = fullfile('EP013','AW184','20210429');


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


contrasts = [0 0.25 0.5 0.75 1];
orientations = linspace(0,330,12);

normalizedResponseCurvesAWAKE = cell(5,2);
% normalizedResponseCurvesANESTH = cell(5,2);

for dp = 1: length(dataPaths)
    dataFile = fullfile('D:\Electrophysiology\',dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise_final','AVMultiContrastDriftingGratingsWhiteNoiseData_selectivity.mat');
    load(dataFile);
    
%     anesthFile = fullfile('D:\Electrophysiology\',dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise_ANESTHETIZED','AVMultiContrastDriftingGratingsWhiteNoiseData_selectivity.mat');
%     anesth = load(anesthFile);
    
    for u =1:length(unitData)
        
        neuronNumber = unitData(u).neuronNumber;
        
        if (ismember(neuronNumber,responsiveUnitsGLM.lightSoundInteractUnits) ||...
                (ismember(neuronNumber,responsiveUnitsGLM.soundResponsiveUnits) &&...
                ismember(neuronNumber,responsiveUnitsGLM.lightResponsiveUnits))) &&...
                (ismember(neuronNumber,responsiveUnits.orientationSelectiveUnitsAND) ||...
                ismember(neuronNumber,responsiveUnits.directionSelectiveUnitsAND)) &&...
                ismember(neuronNumber,[responsiveUnits.singleUnits responsiveUnits.multiUnits])
            [dp neuronNumber]
            allResponses = unitData(u).meanResponse(:,4);
            highContVis = allResponses(((length(contrasts)-1)*length(orientations))+1:((length(contrasts)-1)*length(orientations)+length(orientations)));
            [peak ind] = max(highContVis);
            maxVal = max(allResponses(:));
            minVal = min(allResponses(:));
            
            for c = 1:length(contrasts)
                
                responseVis = allResponses(((c-1)*length(orientations))+1:((c-1)*length(orientations)+length(orientations)));
                responseVisAud = allResponses(((c-1)*length(orientations)+1+length(contrasts)*length(orientations)):((c-1)*length(orientations)+length(contrasts)*length(orientations)+length(orientations)));
                
                %             [peak ind] = max(responseVis);
                targetIndex = floor(length(orientations)/2) + 1;
                
                realignedVis = circshift(responseVis,targetIndex - ind);
                realignedVis = [realignedVis; realignedVis(1)];
                realignedVisAud = circshift(responseVisAud,targetIndex - ind);
                realignedVisAud = [realignedVisAud; realignedVisAud(1)];
                
                
                %             normVis = (realignedVis - minVal) / (maxVal - minVal);
                %             normVisAud = (realignedVisAud - minVal) / (maxVal - minVal);
                normVis = realignedVis ./ peak;
                normVisAud = realignedVisAud ./ peak;
                
                normalizedResponseCurvesAWAKE{c,1} = [normalizedResponseCurvesAWAKE{c,1}; normVis'];
                normalizedResponseCurvesAWAKE{c,2} = [normalizedResponseCurvesAWAKE{c,2}; normVisAud'];
            end
            
            
%             allResponses = anesth.unitData(u).meanResponse(:,4);
%             highContVis = allResponses(((length(contrasts)-1)*length(orientations))+1:((length(contrasts)-1)*length(orientations)+length(orientations)));
%             [peak ind] = max(highContVis);
%             maxVal = max(allResponses(:));
%             minVal = min(allResponses(:));
%             
%             for c = 1:length(contrasts)
%                 
%                 responseVis = allResponses(((c-1)*length(orientations))+1:((c-1)*length(orientations)+length(orientations)));
%                 responseVisAud = allResponses(((c-1)*length(orientations)+1+length(contrasts)*length(orientations)):((c-1)*length(orientations)+length(contrasts)*length(orientations)+length(orientations)));
%                 
%                 %             [peak ind] = max(responseVis);
%                 targetIndex = floor(length(orientations)/2) + 1;
%                 
%                 realignedVis = circshift(responseVis,targetIndex - ind);
%                 realignedVis = [realignedVis; realignedVis(1)];
%                 realignedVisAud = circshift(responseVisAud,targetIndex - ind);
%                 realignedVisAud = [realignedVisAud; realignedVisAud(1)];
%                 
%                 
%                 %             normVis = (realignedVis - minVal) / (maxVal - minVal);
%                 %             normVisAud = (realignedVisAud - minVal) / (maxVal - minVal);
%                 normVis = realignedVis ./ peak;
%                 normVisAud = realignedVisAud ./ peak;
%                 
%                 normalizedResponseCurvesANESTH{c,1} = [normalizedResponseCurvesANESTH{c,1}; normVis'];
%                 normalizedResponseCurvesANESTH{c,2} = [normalizedResponseCurvesANESTH{c,2}; normVisAud'];
%             end
            
        end
    end
end


figure;
for c = 1:5
    vMean = mean(normalizedResponseCurvesAWAKE{c,1});
    avMean = mean(normalizedResponseCurvesAWAKE{c,2});
    vSTD = std(normalizedResponseCurvesAWAKE{c,1})/sqrt(size(normalizedResponseCurvesAWAKE{c,1},1));
    avSTD = std(normalizedResponseCurvesAWAKE{c,2})/sqrt(size(normalizedResponseCurvesAWAKE{c,2},1));
    
    subplot(1,5,c);
    hold on;
    %     plot(-180:30:180,mean(normalizedResponseCurves{c,1}),'Color',[0 0 0]);
    %     plot(-180:30:180,mean(normalizedResponseCurves{c,2}),'Color',[0 0 1]);
    
    
    hA = area(-180:30:180,[vMean-vSTD; 2*vSTD]');
    hA(1).FaceAlpha = 0;
    hA(1).EdgeColor = [1 1 1];
    hA(2).FaceColor = [0 0 0];
    hA(2).FaceAlpha = 0.5;
    hA(2).EdgeColor = [1 1 1];
    hA = area(-180:30:180,[avMean-avSTD; 2*avSTD]');
    hA(1).FaceAlpha = 0;
    hA(1).EdgeColor = [1 1 1];
    hA(2).FaceColor = [0 0 1];
    hA(2).FaceAlpha = 0.5;
    hA(2).EdgeColor = [1 1 1];
    
    plot(-180:30:180,vMean,'Color',[0 0 0]);
    plot(-180:30:180,avMean,'Color',[0 0 1]);
    
    ylim([0 1.2]);
    title(['Contrast = ' num2str(contrasts(c))]);
    xlabel('\Delta Orientation (deg)');
    ylabel('Normalized FR');
end

% 
% figure;
% for c = 1:5
%     v = normalizedResponseCurvesANESTH{c,1};
%     av = normalizedResponseCurvesANESTH{c,2};
%     
%     [row col] = find(isnan(v) | isinf(v)); row = unique(row);
%     v(row,:) = [];
%     [row col] = find(isnan(av) | isinf(av)); row = unique(row);
%     av(row,:) = [];
%     
%     vMean = mean(v);
%     avMean = mean(av);
%     vSTD = std(v)/sqrt(size(v,1));
%     avSTD = std(av)/sqrt(size(av,1));
%     
%     subplot(1,5,c);
%     hold on;
%     %     plot(-180:30:180,mean(normalizedResponseCurves{c,1}),'Color',[0 0 0]);
%     %     plot(-180:30:180,mean(normalizedResponseCurves{c,2}),'Color',[0 0 1]);
%     
%     
%     hA = area(-180:30:180,[vMean-vSTD; 2*vSTD]');
%     hA(1).FaceAlpha = 0;
%     hA(1).EdgeColor = [1 1 1];
%     hA(2).FaceColor = [0 0 0];
%     hA(2).FaceAlpha = 0.5;
%     hA(2).EdgeColor = [1 1 1];
%     hA = area(-180:30:180,[avMean-avSTD; 2*avSTD]');
%     hA(1).FaceAlpha = 0;
%     hA(1).EdgeColor = [1 1 1];
%     hA(2).FaceColor = [0 0 1];
%     hA(2).FaceAlpha = 0.5;
%     hA(2).EdgeColor = [1 1 1];
%     
%     plot(-180:30:180,vMean,'Color',[0 0 0]);
%     plot(-180:30:180,avMean,'Color',[0 0 1]);
%     
%     ylim([0 2]);
%     title(['Contrast = ' num2str(contrasts(c))]);
%     xlabel('\Delta Orientation (deg)');
%     ylabel('Normalized FR');
% end