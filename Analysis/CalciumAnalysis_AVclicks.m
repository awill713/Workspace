
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

responses = zeros(totalNeurons,3,totalFrames);
timeForAverage = 1000; %ms;
framesForAverage = ceil(frameRate*timeForAverage/1000);

stimTypes = unique(stimInfo.order);

mutualInformation = zeros(totalNeurons,2);

slNeurons = [];
sbNeurons = [];
lbNeurons = [];

testSound = [];
testLight = [];
testInteract = [];

groh = [];
grohUp = [];
grohDown = [];

for n = 1:totalNeurons
    
%     respForStats = zeros(repeats,uniqueTones);
%     respForStats2 = zeros(repeats,uniqueTones);
    tempResp = zeros(3,repeats,totalFrames);
    
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
        end
    end
    
%     tempResp = tempResp ./ repeats;
    
    responses(n,:,:) = squeeze(mean(tempResp,2));
    
    soundPre = mean(tempResp(1,:,1:preToneFrames),3);
    lightPre = mean(tempResp(2,:,1:preToneFrames),3);
    bothPre = mean(tempResp(3,:,1:preToneFrames),3);
    
    soundPost = mean(tempResp(1,:,preToneFrames+1:preToneFrames+1+framesForAverage),3);
    lightPost = mean(tempResp(2,:,preToneFrames+1:preToneFrames+1+framesForAverage),3);
    bothPost = mean(tempResp(3,:,preToneFrames+1:preToneFrames+1+framesForAverage),3);
    
    pre = mean(tempResp(:,1:preToneFrames),2);
    post = mean(tempResp(:,preToneFrames+1:preToneFrames+1+framesForAverage),2);
    
    
    [tuningStats(n).soundSig tuningStats(n).soundpValue] = ttest(soundPost);
    [tuningStats(n).lightSig tuningStats(n).lightpValue] = ttest(lightPost);
    [tuningStats(n).bothSig tuningStats(n).bothpValue] = ttest(bothPost);
    
    [tuningStats(n).significant tuningStats(n).pValue] = ttest(post);
    if tuningStats(n).soundSig
        tuningStats(n).soundResponseSign = sign(mean(soundPost)-mean(soundPre));
    end
    if tuningStats(n).lightSig
        tuningStats(n).lightResponseSign = sign(mean(lightPost)-mean(lightPre));
    end
    if tuningStats(n).bothSig
        tuningStats(n).bothResponseSign = sign(mean(bothPost)-mean(bothPre));
    end
    
    [tuningStats(n).AVinteractpValue table stats] = anova1([soundPost' lightPost' bothPost'],[],'off');
    tuningStats(n).AVinteractSig = tuningStats(n).AVinteractpValue < 0.05;
    if tuningStats(n).AVinteractSig
        comparison = multcompare(stats,'display','off');
        tuningStats(n).AVcomparison = comparison(:,6); %row 1 is sound/light, row 2 is sound/both, row 3 is light/both
        tuningStats(n).AVcomparisonSig = comparison(:,6)<0.05; %row 1 is sound/light, row 2 is sound/both, row 3 is light/both
        
        if tuningStats(n).AVcomparisonSig(1)==1
            slNeurons = [slNeurons n];
        end
        if tuningStats(n).AVcomparisonSig(2)==1
            sbNeurons = [sbNeurons n];
        end
        if tuningStats(n).AVcomparisonSig(3)==1
            lbNeurons = [lbNeurons n];
        end
    end
        
    
    audioMean = mean(soundPost);
    audioSTD = std(soundPost);
    bothMean = mean(bothPost);
    bothSTD = std(bothPost);
    meanInput = [audioMean bothMean];
    stdInput = [audioSTD bothSTD];
    [miValues miTotal] = MutualInformation2(meanInput,stdInput); %miTotal is nonsense
    mutualInformation(n,:) = miValues;
    
%     nonePost = zeros(1,repeats);
%     [p table stats] = anova2([nonePost' soundPost'; lightPost' bothPost'],repeats,'off');
%     TesttuningStats(n).statistics = p;
%     if p(1)<0.05
%         TesttuningStats(n).soundSig=1;
%         TesttuningStats(n).soundResponseSign = sign(mean(soundPost)-mean(soundPre));
%     else
%         TesttuningStats(n).soundSig=0;
%     end
%     if p(2)<0.05
%         TesttuningStats(n).lightSig=1;
%         TesttuningStats(n).lightResponseSign = sign(mean(lightPost)-mean(lightPre));
%     else
%         TesttuningStats(n).lightSig=0;
%     end
%     if p(3)<0.05
%         TesttuningStats(n).AVinteractSig = 1;
%         TesttuningStats(n).AVinteractSign = sign(mean(bothPost) - mean(soundPost) - mean(lightPost));
%     else
%         TesttuningStats(n).AVinteractSig=0;
%     end

    [p table stats] = anova2([soundPre' soundPost'; lightPost' bothPost'],repeats,'off');
    if p(1)<0.05
        testSound = [testSound n];
%         if mean(soundPost) > mean(soundPre)
%             testSoundUp = [testSoundUp n];
%         else
%             testSoundDown = [testSoundDown n];
%         end
    end
    if p(2)<0.05
        testLight = [testLight n];
%         if mean(lightPost) > mean(soundPre)
%             testLightUp = [testLightUp n];
%         else
%             testLightDown = [testLightDown n];
%         end
    end
    if p(3)<0.05
        testInteract = [testInteract n];
    end
    
    [h p] = ttest2(soundPost,bothPost);
    if p<0.05
        groh = [groh n];
        if mean(bothPost)>mean(soundPost)
            grohUp = [grohUp n];
        else
            grohDown = [grohDown n];
        end
    end
    
end

soundResponsiveNeurons = find([tuningStats.soundSig]==1);
lightResponsiveNeurons = find([tuningStats.lightSig]==1);
bothResponsiveNeurons = find([tuningStats.bothSig]==1);

% TestsoundResponsiveNeurons = find([TesttuningStats.soundSig]==1);
% TestlightResponsiveNeurons = find([TesttuningStats.lightSig]==1);
% TestitneractNeurons = find([TesttuningStats.AVinteractSig]==1);

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
xticklabels({'-1','0','1','2','3','4','5','6','7'});
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
xticklabels({'-1','0','1','2','3','4','5','6','7'});
xlabel('Time relative to stimulus onset');
ylabel('\Delta F / std(F)');

bothUp = find([tuningStats(bothResponsiveNeurons).bothResponseSign] == 1);
bothDown = find([tuningStats(bothResponsiveNeurons).bothResponseSign] == -1);
bothUp = bothResponsiveNeurons(bothUp);
bothDown = bothResponsiveNeurons(bothDown);
f3 = figure;hold on;
plot(squeeze(mean(responses(bothUp,3,:),1)))
plot(squeeze(mean(responses(bothDown,3,:),1)))
legend({['n = ' num2str(length(bothUp))], ['n = ' num2str(length(bothDown))]});
title('Both responsive neurons');
xticks(1:round(frameRate):totalFrames);
xticklabels({'-1','0','1','2','3','4','5','6','7'});
xlabel('Time relative to stimulus onset');
ylabel('\Delta F / std(F)');

miResponsive = mutualInformation(soundResponsiveNeurons,2)./mutualInformation(soundResponsiveNeurons,1);
miAll = mutualInformation(:,2)./mutualInformation(:,1);
f4 = figure;histogram(miResponsive,10);
xlabel('MI ratio (audiovisual:audio)');
title(['Sound responsive neurons (n = ' num2str(length(soundResponsiveNeurons)) ')'])
f5 = figure;histogram(miAll,10);
xlabel('MI ratio (audiovisual:audio)');
title(['Sound responsive neurons (n = ' num2str(totalNeurons) ')'])
% 
% save([folder '/' mouse '_' date '_AVclickResponses'],'responses','tuningStats','soundResponsiveNeurons','lightResponsiveNeurons','bothResponsiveNeurons','mutualInformation');
% saveas(f1,fullfile(folder,[mouse '_' date '_AVnoise2P_soundResponseFig.fig']));
% saveas(f2,fullfile(folder,[mouse '_' date '_AVnoise2P_lightResponseFig.fig']));
% saveas(f3,fullfile(folder,[mouse '_' date '_AVnoise2P_bothResponseFig.fig']));
% saveas(f4,fullfile(folder,[mouse '_' date '_AVnoise2P_mutualInfoSoundResponsive.fig']));
% saveas(f5,fullfile(folder,[mouse '_' date '_AVnoise2P_mutualInfoAll.fig']));

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
soundRGB = [0 1 0]; %blue
lightRGB = [1 0 0]; %red
bothRGB =  [0 0 1]; %green
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
    if ismember(n,bothResponsiveNeurons)
        neuronColor = neuronColor + bothRGB;
    end
    for p = 1:length(roiLocations{1,n})
        rgbField(roiLocations{1,n}(p,1),roiLocations{1,n}(p,2),:) = neuronColor;
    end
end
f6 = figure;imshow(rgbField);
title('Blue = sound, red = light, green = both');
% saveas(f6,fullfile(folder,[mouse '_' date '_AVnoise2P_neuronsColorCoded.fig']));


%% Color labels superimposed on actual image
theImage = imread('/Volumes/AARON FILES/Two photon/2PM 004/AW069/AV 2P 20190407/20190407_AW069_AVnoise-001_01_AVG.tif');
grayImage = mat2gray(theImage);
rgbImage = ind2rgb(round(grayImage*255),gray);
figure;imshow(rgbImage);

xOffset = spatialInfo.xrange(1)-1;
yOffset = spatialInfo.yrange(1)-1;


roiLocations = cell(1,length(spatialInfo.ROIs));
roiBoundaries = cell(1,length(spatialInfo.ROIs));
for n = 1:length(spatialInfo.ROIs)

    roiRow = mod(spatialInfo.ipix{1,n},size(spatialInfo.im,1));
    roiRow(find(roiRow==0)) = size(spatialInfo.im,1);
    roiColumn = floor(spatialInfo.ipix{1,n}/size(spatialInfo.im,1))+1;

    roiRow = roiRow + yOffset;
    roiColumn = roiColumn + xOffset;
    roiLocations{1,n} = [roiRow roiColumn];
    
    tempImage = zeros(512,512);
    for p = 1:length(roiLocations{1,n})
        tempImage(roiLocations{1,n}(p,1),roiLocations{1,n}(p,2)) = 1;
    end
    boundaries = bwboundaries(tempImage);
    roiBoundaries{1,n} = boundaries{1,1};
end

soundRGB = [0 0 1]; %blue
lightRGB = [1 0 0]; %red
interactRGB =  [0 1 0]; %green

for n = 1:totalNeurons
    neuronColor = [0 0 0];
    if ismember(n,testSound)
        neuronColor = neuronColor + soundRGB;
    end
    if ismember(n,testLight)
        neuronColor = neuronColor + lightRGB;
    end
    if ismember(n,testInteract)
        neuronColor = neuronColor + interactRGB;
    end
    for p = 1:length(roiBoundaries{1,n})
        rgbImage(roiBoundaries{1,n}(p,1),roiBoundaries{1,n}(p,2),:) = neuronColor;
    end
end
figure;imshow(rgbImage);
title('Blue = sound responsive, Red = light responsive, Green = sound/light interaction');