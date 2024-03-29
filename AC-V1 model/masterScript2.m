
totalUnits = 10;

repeats = 5;

contrasts = [0 0.25 0.5 0.75 1];
vStim = [zeros(1,1000) ones(1,1000) zeros(1,1000)];
aStim = [zeros(1,1000) ones(1,1000) zeros(1,1000)];

totalTime = length(vStim);

meanFR = zeros(2,length(contrasts));
stdFR = zeros(2,length(contrasts));


%%
percentOrientation = 0.3;
orientationSelectiveUnits = round(totalUnits*percentOrientation);

orientations = 0:30:330;
defaultTuningCurve = abs(sin(orientations*2*pi/360));
prefOrientations = zeros(1,totalUnits);
tuningCurves = zeros(orientationSelectiveUnits,length(orientations));

for n = 1:orientationSelectiveUnits
    shift = randsample(0:length(orientations)-1,1);
    temp = circshift(defaultTuningCurve,shift);
    prefOrientations(n) = orientations(mod(4+shift-1,length(orientations))+1);
    tuningCurves(n,:) = temp;
end

%%
stimResponses = cell(2,length(contrasts));
for c = 1:length(contrasts)
    vResponses = zeros(totalUnits,length(orientations),repeats);
    avResponses = zeros(totalUnits,length(orientations),repeats);
    
    for o = 1:length(orientations)
        for r = 1:repeats
            acNeurons = zeros(totalUnits,totalTime);
            v1Neurons = zeros(totalUnits,totalTime);
            av1Neurons = zeros(totalUnits,totalTime);
            for t = 1:totalTime
                [c o r t]
                for a = 1:size(acNeurons,1)
                    in = aStim(t);
                    prob = 0.2*exp(in*10-5)./(exp(in*10-5)+1);
                    acNeurons(a,t) = randsample([0 1],1,true,[1-prob prob]);
                end
                
                for n = 1:size(v1Neurons,1)
                    if n<=orientationSelectiveUnits
                        in = 0.5*vStim(t)*contrasts(c)*tuningCurves(n,o);
                    else
                        in = 0.5*vStim(t)*contrasts(c);
                    end
                    prob = 0.2*exp(in*10-5)./(exp(in*10-5)+1);
                    
                    firstBin = max([1 t-500]);
                    recentSpikes = sum(v1Neurons(n,firstBin:t))/20;
                    refract = 0.75*exp(-10*(recentSpikes-1))./(exp(-10*(recentSpikes-1))+1)+0.25;
                    
                    prob = prob*refract;
                    v1Neurons(n,t) = randsample([0 1],1,true,[1-prob prob]);
                end
                
                for n = 1:size(av1Neurons,1)
                    recentBin = max([1 t-10]);
                    acActivity = mean(sum(acNeurons(:,recentBin:t),2));
                    acInput = acActivity/20;
                    
                    if n<=orientationSelectiveUnits
                        in = 0.5*vStim(t)*contrasts(c)*tuningCurves(n,o) + acInput;
                    else
                        in = 0.5*vStim(t)*contrasts(c) + acInput;
                    end
                    prob = 0.2*exp(in*10-5)./(exp(in*10-5)+1);
                    
                    firstBin = max([1 t-500]);
                    recentSpikes = sum(av1Neurons(n,firstBin:t))/20;
                    refract = 0.75*exp(-10*(recentSpikes-1))./(exp(-10*(recentSpikes-1))+1)+0.25;
                    
                    prob = prob*refract;
                    av1Neurons(n,t) = randsample([0 1],1,true,[1-prob prob]);
                end
            end
            
            v1FR = zeros(size(v1Neurons));
            for t = 1:totalTime
                for n = 1:size(v1Neurons,1)
                    firstBin = max([1 t-9]);
                    spikes = sum(v1Neurons(n,firstBin:t));
                    fr = spikes*100;
                    v1FR(n,t) = fr;
                end
            end
            meanFR(1,c) = mean(mean(v1FR(:,1000:1300),2));
            stdFR(1,c) = std(mean(v1FR(:,1000:1300),2));
            varFR(1,c) = std(mean(v1FR(:,1000:1300),2)).^2;
            
            vResponses(:,o,r) = mean(v1FR(:,1000:1300),2);
            
            
            av1FR = zeros(size(av1Neurons));
            for t = 1:totalTime
                for n = 1:size(av1Neurons,1)
                    firstBin = max([1 t-9]);
                    spikes = sum(av1Neurons(n,firstBin:t));
                    fr = spikes*100;
                    av1FR(n,t) = fr;
                end
            end
            meanFR(2,c) = mean(mean(av1FR(:,1000:1300),2));
            stdFR(2,c) = std(mean(av1FR(:,1000:1300),2));
            varFR(2,c) = std(mean(av1FR(:,1000:1300),2)).^2;
            
            avResponses(:,o,r) = mean(av1FR(:,1000:1300),2);
        end
    end
    
    stimResponses{1,c} = vResponses;
    stimResponses{2,c} = avResponses;
end

%%
figure;
tc = mean(vResponses(1,:,:),3);
polarplot(orientations*2*pi/360,tc);

%%
choiceMap = cell(2,length(contrasts)); %v and av, by contrasts
accuracy = cell(2,length(contrasts));
orientationMap = cell(1,length(contrasts));

for n = 1:orientationSelectiveUnits
    
    for c = 1:length(contrasts)
        
        %%Overall accuracy
        tempChoiceMap = zeros(2,length(orientations),length(orientations));
        for testDir = 1:length(orientations)
            [n c testDir]
            
            mleClassification = zeros(2,length(orientations));
            
            
            for trial = 1:repeats
                if exist('neuronStats','var')
                    clear neuronStats
                end
                
                for orient = 1:length(orientations)
                    
                    %                 vInd = (c-1)*length(orientations) + orient;
                    %                 avInd = (c-1)*length(orientations) + orient + length(orientations)*length(contrasts);
                    
                    trialsIncluded = 1:repeats;
                    if orient==testDir
                        trialsIncluded  = setdiff(trialsIncluded,trial);
                        
                        probeTrialV = squeeze(stimResponses{1,c}(n,orient,trial));
                        probeTrialAV = squeeze(stimResponses{2,c}(n,orient,trial));
                        probeData = [probeTrialV; probeTrialAV];
                        
                        %                     probeTrialV = unitData(u).trialResponse(vInd,trial);
                        %                     probeTrialAV = unitData(u).trialResponse(avInd,trial);
                        %                     probeData = [probeTrialV; probeTrialAV];
                    end
                    
                    vData = squeeze(stimResponses{1,c}(n,orient,trialsIncluded));
                    neuronStats(orient,1) = fitdist(vData,'Normal');
                    
                    avData = squeeze(stimResponses{2,c}(n,orient,trialsIncluded));
                    neuronStats(orient,2) = fitdist(avData,'Normal');
                    
                    %                 vData = unitData(u).trialResponse(vInd,trialsIncluded);
                    %                 avData = unitData(u).trialResponse(avInd,trialsIncluded);
                    %
                    %                 neuronStats(orient,1) = fitdist(vData','Normal');
                    %                 neuronStats(orient,2) = fitdist(avData','Normal');
                end
                
                vEstimate = maximumLikelihoodFunctionSingle(neuronStats(:,1),probeData(1));
                mleClassification(1,vEstimate) = mleClassification(1,vEstimate)+1;
                
                avEstimate = maximumLikelihoodFunctionSingle(neuronStats(:,2),probeData(2));
                mleClassification(2,avEstimate) = mleClassification(2,avEstimate) +1;
            end
            
            tempChoiceMap(1,testDir,:) = mleClassification(1,:) / sum(mleClassification(1,:));
            tempChoiceMap(2,testDir,:) = mleClassification(2,:) / sum(mleClassification(2,:));
        end
        
        choiceMap{1,c} = cat(3,choiceMap{1,c},squeeze(tempChoiceMap(1,:,:)));
        choiceMap{2,c} = cat(3,choiceMap{2,c},squeeze(tempChoiceMap(2,:,:)));
        
        accuracy{1,c} = [accuracy{1,c}; mean(diag(squeeze(tempChoiceMap(1,:,:)))) n];
        accuracy{2,c} = [accuracy{2,c}; mean(diag(squeeze(tempChoiceMap(2,:,:)))) n];
        
        
        %%Orientation accuracy
        %     vIndBase = (c-1)*length(orientations) + 1;
        %     avIndBase = (c-1)*length(orientations) + 1 + length(orientations)*length(contrasts);
        %
        %     vResp = unitData(u).meanResponse(vIndBase:vIndBase+length(orientations)-1,4);
        %     avResp = unitData(u).meanResponse(avIndBase:avIndBase+length(orientations)-1,4);
        %     [~, prefOrientV] = max(vResp);
        %     [~, prefOrientAV] = max(avResp);
        %
        %     orthIndV = [mod(prefOrientV+3-1,length(orientations)) mod(prefOrientV-3-1,length(orientations))]+1;
        %     orthIndAV = [mod(prefOrientAV+3-1,length(orientations)) mod(prefOrientAV-3-1,length(orientations))]+1;
        %
        prefInd = find(orientations==prefOrientations(n));
        orthInd = [mod(prefInd+3-1,length(orientations)) mod(prefInd-3-1,length(orientations))]+1;
        
        orientClassification = zeros(2,2);
        for trial = 1:repeats
            if exist('neuronStats','var')
                clear neuronStats
            end
            
            probeTrialV = stimResponses{1,c}(n,prefInd,trial);
            probeTrialAV = stimResponses{2,c}(n,prefInd,trial);
            probeData = [probeTrialV; probeTrialAV];
            
            prefDataV = squeeze(stimResponses{1,c}(n,prefInd,setdiff(1:repeats,trial)));
            prefDataAV = squeeze(stimResponses{2,c}(n,prefInd,setdiff(1:repeats,trial)));
            neuronStats(1,1) = fitdist(prefDataV,'Normal');
            neuronStats(1,2) = fitdist(prefDataAV,'Normal');
            
            orthDataV = reshape(stimResponses{1,c}(n,orthInd,:),[2*repeats 1]);
            orthDataAV = reshape(stimResponses{2,c}(n,orthInd,:),[2*repeats 1]);
            neuronStats(2,1) = fitdist(orthDataV,'Normal');
            neuronStats(2,2) = fitdist(orthDataAV,'Normal');
            
            
            %         probeTrialV = unitData(u).trialResponse(vIndBase + prefOrientV,trial);
            %         probeTrialAV = unitData(u).trialResponse(avIndBase +prefOrientAV,trial);
            %         probeData = [probeTrialV; probeTrialAV];
            
            %         prefDataV = unitData(u).trialResponse(vIndBase+prefOrientV,setdiff(1:repeats,trial));
            %         prefDataAV = unitData(u).trialResponse(avIndBase+prefOrientAV,setdiff(1:repeats,trial));
            %         neuronStats(1,1) = fitdist(prefDataV','Normal');
            %         neuronStats(1,2) = fitdist(prefDataAV','Normal');
            
            %         orthDataV = reshape(unitData(u).trialResponse(vIndBase+orthIndV,:),[1 2*repeats]);
            %         orthDataAV = reshape(unitData(u).trialResponse(avIndBase+orthIndAV,:),[1 2*repeats]);
            %         neuronStats(2,1) = fitdist(orthDataV','Normal');
            %         neuronStats(2,2) = fitdist(orthDataAV','Normal');
            
            vEstimate = maximumLikelihoodFunctionSingle(neuronStats(:,1),probeData(1));
            orientClassification(1,vEstimate) = orientClassification(1,vEstimate)+1;
            
            avEstimate = maximumLikelihoodFunctionSingle(neuronStats(:,2),probeData(2));
            orientClassification(2,avEstimate) = orientClassification(2,avEstimate)+1;
        end
        
        orientationMap{1,c} = [orientationMap{1,c}; orientClassification(:,1)'./repeats];
    end
end


for c = 1:length(contrasts)
    vAcc(c) = mean(accuracy{1,c}(:,1));
    vAccStd(c) = std(accuracy{1,c}(:,1))/sqrt(size(accuracy{1,c},1));
    avAcc(c) = mean(accuracy{2,c}(:,1));
    avAccStd(c) = std(accuracy{2,c}(:,1))/sqrt(size(accuracy{2,c},1));
    
    dAcc(c) = mean(accuracy{2,c}(:,1)-accuracy{1,c}(:,1));
    dAccStd(c) = std(accuracy{2,c}(:,1)-accuracy{1,c}(:,1))/sqrt(size(accuracy{2,c},1));
end
figure;hold on;
hA = area(contrasts,[vAcc-vAccStd; 2*vAccStd]')
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
hA = area(contrasts,[avAcc-avAccStd; 2*avAccStd]')
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 1];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];

plot(contrasts,vAcc,'Color',[0 0 0]);
plot(contrasts,avAcc,'Color',[0 0 1]);

figure;hold on;
hA = area(contrasts,[dAcc-dAccStd; 2*dAccStd]')
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];

plot(contrasts,dAcc,'Color',[0 0 0]);




figure;
subplot(1,2,1);
imagesc(mean(choiceMap{1,5},3));
subplot(1,2,2);
imagesc(mean(choiceMap{2,5},3));






    
%     av1FR = zeros(size(av1Neurons));
%     for t = 1:totalTime
%         for n = 1:size(av1Neurons,1)
%             firstBin = max([1 t-9]);
%             spikes = sum(av1Neurons(n,firstBin:t));
%             fr = spikes*100;
%             av1FR(n,t) = fr;
%         end
%     end
%     meanFR(2,c) = mean(mean(av1FR(:,1000:1300),2));
%     stdFR(2,c) = std(mean(av1FR(:,1000:1300),2));
%     varFR(2,c) = std(mean(av1FR(:,1000:1300),2)).^2;

%%


figure;hold on;
plot(contrasts,meanFR(1,:),'Color',[0 0 0]);
plot(contrasts,stdFR(1,:),'Color',[0 0 0]);
plot(contrasts,meanFR(2,:),'Color',[0 0 1]);
plot(contrasts,stdFR(2,:),'Color',[0 0 1]);

figure;hold on;
plot(contrasts,stdFR(1,:)./meanFR(1,:),'Color',[0 0 0]);
plot(contrasts,stdFR(2,:)./meanFR(2,:),'Color',[0 0 1]);
title('Coefficient of Variation');

figure;hold on;
plot(contrasts,varFR(1,:)./meanFR(1,:),'Color',[0 0 0]);
plot(contrasts,varFR(2,:)./meanFR(2,:),'Color',[0 0 1]);
title('Fano factor');