
totalUnits = length(unitData);


normPupil = zeros(120,10,121);
pupilComp = zeros(120,10);
for ii = 1:size(pupilData,1)
    for rr = 1:size(pupilData,2)
        trace = pupilData(ii,rr,:);
        theMean = nanmean(trace(1:30));
        tempTrace = (trace-theMean)./theMean;
        normPupil(ii,rr,:) = tempTrace;
        pupilComp(ii,rr) = nanmean(tempTrace(41:121));
    end
end

popData = zeros(0,120,10);
for u = 1:totalUnits
    neuronNumber = unitData(u).neuronNumber;
    
    if (ismember(neuronNumber,responsiveUnitsGLM.lightSoundInteractUnits) ||...
            (ismember(neuronNumber,responsiveUnitsGLM.soundResponsiveUnits) &&...
            ismember(neuronNumber,responsiveUnitsGLM.lightResponsiveUnits))) &&...
            ismember(neuronNumber,[responsiveUnits.singleUnits responsiveUnits.multiUnits])
        
        responses = unitData(u).trialResponse;
        
        linResponses = reshape(responses,[1 size(responses,1)*size(responses,2)]);
        avgResp = mean(linResponses);
        stdResp = std(linResponses);
        
        zscoredResp = zeros([1 size(responses)]);
        for t = 1:size(responses,1)
            for r = 1:size(responses,2)
                temp = responses(t,r);
                z = (temp-avgResp)./stdResp;
                zscoredResp(1,t,r) = z;
            end
        end
        
        popData = cat(1,popData,zscoredResp);
    end  
end

meanZ = squeeze(mean(popData,1));
linMeanZ = reshape(meanZ,[1 1200]);
linPupil = reshape(pupilComp,[1 1200]);
figure;scatter(linPupil,linMeanZ)
exclude = find(isnan(linPupil));
linPupil(exclude) = []; linMeanZ(exclude) = [];
[CC pp] = corrcoef(linPupil,linMeanZ);
CC = CC(1,2)
pp = pp(1,2)
title(['CC = ' num2str(CC) ', p-value = ' num2str(pp)])

figure;
for c = 1:5
    vInd = (c-1)*12+1 : c*12;
    avInd = vInd + 60;
    
    neuralV = reshape(meanZ(vInd,:),[1 120]);
    pupV = reshape(pupilComp(vInd,:),[1 120]);
    neuralAV = reshape(meanZ(avInd,:),[1 120]);
    pupAV = reshape(pupilComp(avInd,:),[1 120]);
    subplot(1,5,c);hold on;
    scatter(pupV,neuralV,[],[0 0 0]);
    scatter(pupAV,neuralAV,[],[0 0 1]);
    neuralBoth = [neuralV neuralAV];
    pupBoth = [pupV pupAV];
    exclude = find(isnan(pupBoth));
    pupBoth(exclude) = []; neuralBoth(exclude) = [];
    [CC pp] = corrcoef(pupBoth,neuralBoth);
    CC = CC(1,2);
    pp = pp(1,2);
    xlabel('Mean pupil growth (%baseline)');
    ylabel('Z-scored population neural response');
    title(['CC = ' num2str(CC) ', p = ' num2str(pp)])
%     p = polyfit(pupBoth,neuralBoth,1);
%     minx = min(pupBoth);maxx = max(pupBoth);
%     miny = minx*p(1) + p(2);
%     maxy = maxx*p(1) + p(2);
%     line([minx maxx],[miny maxy]);
end
