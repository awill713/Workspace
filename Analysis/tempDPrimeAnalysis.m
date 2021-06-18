


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

dPrimeCurves = cell(5,2);

for dp = 1: length(dataPaths)
    dataFile = fullfile('D:\Electrophysiology\',dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise_final','AVMultiContrastDriftingGratingsWhiteNoiseData_selectivity.mat');
    load(dataFile);
    

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
            [peak vInd] = max(highContVis);
            highContAudVis = allResponses(((length(contrasts)-1)*length(orientations))+61:((length(contrasts)-1)*length(orientations)+length(orientations)+60));
            [peak avInd] = max(highContAudVis);
            
            vIndices = circshift(1:12,7-vInd);
            avIndices = circshift(1:12,7-avInd);
            
            for c = 1:length(contrasts)
                
                vBaseInd = (c-1)*length(orientations);
                avBaseInd = vBaseInd + 60;
                
                tempCurves = zeros(2,12);
                for dd = 1:12
                    theInd = vIndices(dd);
                    bestInd = vIndices(7);
                    
                    v1Data = unitData(u).trialResponse(vBaseInd + theInd,:);
                    v2Data = unitData(u).trialResponse(vBaseInd + bestInd,:);
                    v1Mean = mean(v1Data);
                    v1Std = std(v1Data);
                    v2Mean = mean(v2Data);
                    v2Std = std(v2Data);
                    
                    theInd = avIndices(dd);
                    bestInd = avIndices(7);
                    av1Data = unitData(u).trialResponse(avBaseInd + theInd,:);
                    av2Data = unitData(u).trialResponse(avBaseInd + bestInd,:);
                    av1Mean = mean(av1Data);
                    av1Std = std(av1Data);
                    av2Mean = mean(av2Data);
                    av2Std = std(av2Data);
                    
                    vPrime = (v2Mean - v1Mean) / sqrt(0.5* (v2Std^2 + v1Std^2));
                    avPrime = (av2Mean - av1Mean) / sqrt(0.5* (av2Std^2 + av1Std^2));
                    
                    tempCurves(:,dd) = [vPrime avPrime]';
                end
                tempCurves = [tempCurves tempCurves(:,1)];
                
                dPrimeCurves{c,1} = [dPrimeCurves{c,1}; tempCurves(1,:)];
                dPrimeCurves{c,2} = [dPrimeCurves{c,2}; tempCurves(2,:)];
            end
        end
    end
end


figure;
for c = 1:5
    vMean = mean(dPrimeCurves{c,1});
    avMean = mean(dPrimeCurves{c,2});
    vSTD = std(dPrimeCurves{c,1})/sqrt(size(dPrimeCurves{c,1},1));
    avSTD = std(dPrimeCurves{c,2})/sqrt(size(dPrimeCurves{c,2},1));
    
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
    
    ylim([-0.3 1.6]);
    title(['Contrast = ' num2str(contrasts(c))]);
    xlabel('\Delta Orientation (deg)');
    ylabel('D-prime selectivity index');
end
                