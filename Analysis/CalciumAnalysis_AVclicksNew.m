
%Analysis for 4-indices (sound, light, both, neither) and 2-way ANOVA

clear;
[file folder] = uigetfile('/Users/Aaron/Documents/'); %looking for *processedData.mat
if file==0
    return;
end
load([folder '/' file]);

mouse = exptInfo.mouse;
date = exptInfo.recDate;

totalNeurons = length(calcium.n);

% totalTones = length(stimInfo.order);
repeats = stimInfo.repeats;

frameRate = exptInfo.fr;
preToneTime = 1000; %ms
preToneFrames = round(frameRate*preToneTime/1000);
postToneTime = stimInfo.interClickInterval*1000; %ms
postToneFrames = ceil(frameRate*postToneTime/1000);
totalFrames = postToneFrames + preToneFrames + 1;

responses = zeros(totalNeurons,4,totalFrames);
timeForAverage = 1000; %ms;
framesForAverage = ceil(frameRate*timeForAverage/1000);

stimTypes = unique(stimInfo.order);

mutualInformation = zeros(totalNeurons,2);

slNeurons = [];
sbNeurons = [];
lbNeurons = [];


stimulusMatrix = zeros(2,length(stimInfo.order));
for e = 1:length(stimInfo.order)
    if stimInfo.order(e)==1
        stimulusMatrix(:,e) = [1;0];
    elseif stimInfo.order(e)==2
        stimulusMatrix(:,e) = [0;1];
    elseif stimInfo.order(e)==3
        stimulusMatrix(:,e) = [1;1];
    elseif stimInfo.order(e)==4
        stimulusMatrix(:,e) = [0;0];
    end
end

respMeanVar = zeros(totalNeurons,4,2);
for n = 1:totalNeurons
    responseMatrix = zeros(1,length(stimInfo.order));
%     respForStats = zeros(repeats,uniqueTones);
%     respForStats2 = zeros(repeats,uniqueTones);
    tempResp = zeros(4,repeats,totalFrames);
    
    for type = 1:length(stimTypes)
        eventsOfType = find(stimInfo.order == type);
        for r = 1:repeats
            ev = eventsOfType(r);
            firstFrame = events.eventsOn(ev) - preToneFrames;
            lastFrame = events.eventsOn(ev) + postToneFrames;
            
            pretoneAverage = mean(calcium.npilSubTraces(n,firstFrame:events.eventsOn(ev)-1));
            pretoneSTD = std(calcium.npilSubTraces(n,firstFrame:events.eventsOn(ev)-1));
            
            trace = calcium.npilSubTraces(n,firstFrame:lastFrame) - pretoneAverage;
            trace = trace ./ pretoneSTD;
            
            tempResp(type,r,:) = trace;
            
            responseMatrix(1,ev) = mean(trace(preToneFrames+1:preToneFrames+1+framesForAverage));
        end
    end
    
%     tempResp = tempResp ./ repeats;
    
    responses(n,:,:) = squeeze(mean(tempResp,2));
    
    soundPre = mean(tempResp(1,:,1:preToneFrames),3);
    lightPre = mean(tempResp(2,:,1:preToneFrames),3);
    bothPre = mean(tempResp(3,:,1:preToneFrames),3);
    nonePre = mean(tempResp(4,:,1:preToneFrames),3);
    
    soundPost = mean(tempResp(1,:,preToneFrames+1:preToneFrames+1+framesForAverage),3);
    lightPost = mean(tempResp(2,:,preToneFrames+1:preToneFrames+1+framesForAverage),3);
    bothPost = mean(tempResp(3,:,preToneFrames+1:preToneFrames+1+framesForAverage),3);
    nonePost = mean(tempResp(4,:,preToneFrames+1:preToneFrames+1+framesForAverage),3);

    [p table stats] = anova2([nonePost' soundPost'; lightPost' bothPost'],repeats,'off');
    tuningStats(n).statistics = p;
    if p(1)<0.05
        tuningStats(n).soundSig=1;
        tuningStats(n).soundResponseSign = sign(mean(soundPost)-mean(soundPre));
    else
        tuningStats(n).soundSig=0;
    end
    if p(2)<0.05
        tuningStats(n).lightSig=1;
        tuningStats(n).lightResponseSign = sign(mean(lightPost)-mean(lightPre));
    else
        tuningStats(n).lightSig=0;
    end
    if p(3)<0.05
        tuningStats(n).AVinteractSig = 1;
        tuningStats(n).AVinteractSign = sign(mean(bothPost) - mean(soundPost) - mean(lightPost));
    else
        tuningStats(n).AVinteractSig=0;
    end
    
    respMeanVar(n,1,1) = mean(soundPost); respMeanVar(n,1,2) = std(soundPost);
    respMeanVar(n,2,1) = mean(lightPost); respMeanVar(n,2,2) = std(lightPost);
    respMeanVar(n,3,1) = mean(bothPost); respMeanVar(n,3,2) = std(bothPost);
    respMeanVar(n,4,1) = mean(nonePost); respMeanVar(n,4,2) = std(nonePost);
        
    
    audioMean = mean(soundPost);
    audioSTD = std(soundPost);
    bothMean = mean(bothPost);
    bothSTD = std(bothPost);
    meanInput = [audioMean bothMean];
    stdInput = [audioSTD bothSTD];
    [miValues miTotal] = MutualInformation2(meanInput,stdInput); %miTotal is nonsense
    mutualInformation(n,:) = miValues;
    
    glm = fitglm(stimulusMatrix',responseMatrix','interactions','distr','normal');
    tuningStatsGLM(n).table = glm.Coefficients;
    tuningStatsGLM(n).estimates = glm.Coefficients{:,1};
    tuningStatsGLM(n).significance = glm.Coefficients{:,4};
    if tuningStatsGLM(n).significance(2) <0.05
        tuningStatsGLM(n).soundSig = 1;
        tuningStatsGLM(n).soundResponseSign = sign(tuningStatsGLM(n).estimates(2));
    else
        tuningStatsGLM(n).soundSig = 0;
    end
    if tuningStatsGLM(n).significance(3) <0.05
        tuningStatsGLM(n).lightSig = 1;
        tuningStatsGLM(n).lightResponseSign = sign(tuningStatsGLM(n).estimates(3));
    else
        tuningStatsGLM(n).lightSig = 0;
    end
    if tuningStatsGLM(n).significance(4) <0.05
        tuningStatsGLM(n).AVinteractSig = 1;
        tuningStatsGLM(n).AVinteractSign = sign(tuningStatsGLM(n).estimates(4));
    else
        tuningStatsGLM(n).AVinteractSig = 0;
    end

    
    [TesttuningStats(n).soundSig TesttuningStats(n).soundpValue] = ttest2(soundPost,nonePost);
    [TesttuningStats(n).lightSig TesttuningStats(n).lightpValue] = ttest2(lightPost,nonePost);
    [TesttuningStats(n).bothSig TesttuningStats(n).bothpValue] = ttest2(bothPost,nonePost);
    if TesttuningStats(n).soundSig
        TesttuningStats(n).soundResponseSign = sign(mean(soundPost)-mean(soundPre));
    end
    if TesttuningStats(n).lightSig
        TesttuningStats(n).lightResponseSign = sign(mean(lightPost)-mean(lightPre));
    end
    if TesttuningStats(n).bothSig
        TesttuningStats(n).bothResponseSign = sign(mean(bothPost)-mean(bothPre));
    end
    
end

soundResponsiveNeurons = find([tuningStats.soundSig]==1);
lightResponsiveNeurons = find([tuningStats.lightSig]==1);
interactNeurons = find([tuningStats.AVinteractSig]==1);
soundUp = find([tuningStats(soundResponsiveNeurons).soundResponseSign] == 1);
soundDown = find([tuningStats(soundResponsiveNeurons).soundResponseSign] == -1);
soundUp = soundResponsiveNeurons(soundUp);
soundDown = soundResponsiveNeurons(soundDown);

soundResponsiveNeuronsGLM = find([tuningStatsGLM.soundSig]==1);
lightResponsiveNeuronsGLM = find([tuningStatsGLM.lightSig]==1);
interactNeuronsGLM = find([tuningStatsGLM.AVinteractSig]==1);
soundUpGLM = find([tuningStatsGLM(soundResponsiveNeuronsGLM).soundResponseSign] == 1);
soundDownGLM = find([tuningStatsGLM(soundResponsiveNeuronsGLM).soundResponseSign] == -1);
soundUpGLM = soundResponsiveNeuronsGLM(soundUpGLM);
soundDownGLM = soundResponsiveNeuronsGLM(soundDownGLM);

TestsoundResponsiveNeurons = find([TesttuningStats.soundSig]==1);
TestlightResponsiveNeurons = find([TesttuningStats.lightSig]==1);
TestsoundUp = find([TesttuningStats(TestsoundResponsiveNeurons).soundResponseSign] == 1);
TestsoundDown = find([TesttuningStats(TestsoundResponsiveNeurons).soundResponseSign] == -1);
TestsoundUp = TestsoundResponsiveNeurons(TestsoundUp);
TestsoundDown = TestsoundResponsiveNeurons(TestsoundDown);


soundUp = find([tuningStats(soundResponsiveNeurons).soundResponseSign] == 1);
soundDown = find([tuningStats(soundResponsiveNeurons).soundResponseSign] == -1);
soundUp = soundResponsiveNeurons(soundUp);
soundDown = soundResponsiveNeurons(soundDown);
f1 = figure;hold on;
plot(squeeze(mean(responses(soundUp,1,:),1)))
plot(squeeze(mean(responses(soundDown,1,:),1)))
legend({['n = ' num2str(length(soundUp))], ['n = ' num2str(length(soundDown))]});
title('Sound responsive neurons');
xticks(1:round(frameRate):totalFrames);
xticklabels({'-1','0','1','2','3','4'});
xlabel('Time relative to stimulus onset');
ylabel('\Delta F / std(F)');

lightUp = find([tuningStats(lightResponsiveNeurons).lightResponseSign] == 1);
lightDown = find([tuningStats(lightResponsiveNeurons).lightResponseSign] == -1);
lightUp = lightResponsiveNeurons(lightUp);
lightDown = lightResponsiveNeurons(lightDown);
f2 = figure;hold on;
plot(squeeze(mean(responses(lightUp,2,:),1)))
plot(squeeze(mean(responses(lightDown,2,:),1)))
legend({['n = ' num2str(length(lightUp))], ['n = ' num2str(length(lightDown))]});
title('Light responsive neurons');
xticks(1:round(frameRate):totalFrames);
xticklabels({'-1','0','1','2','3'});
xlabel('Time relative to stimulus onset');
ylabel('\Delta F / std(F)');

avInteractUp = find([tuningStats(interactNeurons).AVinteractSign] == 1);
avInteractDown = find([tuningStats(interactNeurons).AVinteractSign] == -1);
avInteractUp = interactNeurons(avInteractUp);
avInteractDown = interactNeurons(avInteractDown);
avUpExpected = squeeze(mean(responses(avInteractUp,1,:) + responses(avInteractUp,2,:),1));
avDownExpected = squeeze(mean(responses(avInteractDown,1,:) + responses(avInteractDown,2,:),1));
f3 = figure;hold on;
plot(squeeze(mean(responses(avInteractUp,3,:),1)))
plot(avUpExpected);
title('Enhanced AV interaction');
xticks(1:round(frameRate):totalFrames);
xticklabels({'-1','0','1','2','3'});
xlabel('Time relative to stimulus onset');
ylabel('\Delta F / std(F)');
legend({['n = ' num2str(length(avInteractUp))]});
f4 = figure;hold on;
plot(squeeze(mean(responses(avInteractDown,3,:),1)))
plot(avDownExpected);
title('Diminished AV interaction');
xticks(1:round(frameRate):totalFrames);
xticklabels({'-1','0','1','2','3'});
xlabel('Time relative to stimulus onset');
ylabel('\Delta F / std(F)');
legend({['n = ' num2str(length(avInteractDown))]});

miResponsive = mutualInformation(soundResponsiveNeurons,2)./mutualInformation(soundResponsiveNeurons,1);
miAll = mutualInformation(:,2)./mutualInformation(:,1);
f5 = figure;histogram(miResponsive,10);
xlabel('MI ratio (audiovisual:audio)');
title(['Sound responsive neurons (n = ' num2str(length(soundResponsiveNeurons)) ')'])
f6 = figure;histogram(miAll,10);
xlabel('MI ratio (audiovisual:audio)');
title(['All neurons (n = ' num2str(totalNeurons) ')'])
% 
% save([folder '/' mouse '_' date '_AVclickResponses'],'responses','tuningStats','soundResponsiveNeurons','lightResponsiveNeurons','interactNeurons','mutualInformation');
% saveas(f1,fullfile(folder,[mouse '_' date '_AVnoise2P_soundResponseFig.fig']));
% saveas(f2,fullfile(folder,[mouse '_' date '_AVnoise2P_lightResponseFig.fig']));
% saveas(f3,fullfile(folder,[mouse '_' date '_AVnoise2P_AVinteractUpFig.fig']));
% saveas(f4,fullfile(folder,[mouse '_' date '_AVnoise2P_AVinteractDownFig.fig']));
% saveas(f5,fullfile(folder,[mouse '_' date '_AVnoise2P_mutualInfoSoundResponsive.fig']));
% saveas(f6,fullfile(folder,[mouse '_' date '_AVnoise2P_mutualInfoAll.fig']));

% 
% barResponses = zeros(52,3);
% 
% for n = 1:52
%     for e = 1:3
%         tri = zeros(1,20);
%          
%         barResponses(n,e) = squeeze(mean(responses(n,e,31:60),3));
%     end
% end

xOffset = spatialInfo.xrange(1)-1;
yOffset = spatialInfo.yrange(1)-1;

tempImage = zeros(1,size(spatialInfo.im,1)*size(spatialInfo.im,2));
roiLocations = cell(1,length(spatialInfo.ROIs));
for n = 1:length(spatialInfo.ROIs)
    tempImage(spatialInfo.ipix{1,n})=1;
    
    roiRow = mod(spatialInfo.ipix{1,n},size(spatialInfo.im,1));
    roiRow(find(roiRow==0)) = size(spatialInfo.im,1);
    roiColumn = floor(spatialInfo.ipix{1,n}/size(spatialInfo.im,1))+1;

    roiRow = roiRow + yOffset;
    roiColumn = roiColumn + xOffset;
    roiLocations{1,n} = [roiRow roiColumn];
    
end

rgbField = zeros(512,512,3);
rgbField(:,:,:) = 0.7;
soundRGB = [0 0 1]; %blue
lightRGB = [1 0 0]; %red
bothRGB =  [0 1 0]; %green
% slRGB =  
% sbRGB = 
% lbRGB = 
% slbRGB = [1 1 1];

for n = 1:totalNeurons
    neuronColor = [0 0 0];
    if ismember(n,soundResponsiveNeurons)
        neuronColor = neuronColor + soundRGB;
    end
    if ismember(n,lightResponsiveNeurons)
        neuronColor = neuronColor + lightRGB;
    end
    if ismember(n,interactNeurons)
        neuronColor = neuronColor + bothRGB;
    end
    for p = 1:length(roiLocations{1,n})
        rgbField(roiLocations{1,n}(p,1),roiLocations{1,n}(p,2),:) = neuronColor;
    end
end
f6 = figure;imshow(rgbField);
title('Blue = sound, red = light, green = interact');
% saveas(f6,fullfile(folder,[mouse '_' date '_AVnoise2P_neuronsColorCoded.fig']));

% 
% for n = 1:totalNeurons
%     tempImage = zeros(512,512,3);
%     for p = 1:length(roiLocations{1,n})
%         tempImage(roiLocations{1,n}(p,1),roiLocations{1,n}(p,2),:) = [1 1 1];
%     end
%     figure;imshow(tempImage);
%     title(num2str(n))
% end
% 
% repetitions=1000;
% meanVar = squeeze(mean(respMeanVar,1));
% sigsCategories = zeros(repetitions,3);
% var1 = fitdist(respMeanVar(:,1,2),'InverseGaussian');
% var2 = fitdist(respMeanVar(:,2,2),'InverseGaussian');
% var3 = fitdist(respMeanVar(:,3,2),'InverseGaussian');
% var4 = fitdist(respMeanVar(:,4,2),'InverseGaussian');
% for rep = 1:repetitions
%     rep
%     for a = 1:totalNeurons
%         soundTemp = var1.random*randn(repeats,1) + meanVar(1,1);
%         lightTemp = var2.random*randn(repeats,1) + meanVar(2,1);
%         bothTemp = var3.random*randn(repeats,1) + meanVar(3,1);
%         noneTemp = var4.random*randn(repeats,1) + meanVar(4,1);
%         
% %         soundTemp = meanVar(1,2)*randn(repeats,1) + meanVar(1,1);
% %         lightTemp = meanVar(2,2)*randn(repeats,1) + meanVar(2,1);
% %         bothTemp = meanVar(3,2)*randn(repeats,1) + meanVar(3,1);
% %         noneTemp = meanVar(4,2)*randn(repeats,1) + meanVar(4,1);
%         
%         [p table stats] = anova2([noneTemp soundTemp; lightTemp bothTemp],repeats,'off');
%         if p(1)<0.05
%             sigsCategories(rep,1) = sigsCategories(rep,1) + 1;
%         end
%         if p(2)<0.05
%             sigsCategories(rep,2) = sigsCategories(rep,2) + 1;
%         end
%         if p(3)<0.05
%             sigsCategories(rep,3) = sigsCategories(rep,3) + 1;
%         end
%     end
% end
% mean(sigsCategories)
% figure;histogram(sigsCategories(:,1));title('Sound');
% figure;histogram(sigsCategories(:,2));title('Light');
% figure;histogram(sigsCategories(:,3));title('Interact');