%% Stimulus-induced locomotion

clear;

saveDir = fullfile('/Volumes/AARON DATA/','Electrophysiology','EP010');

dataPaths{1} = fullfile('AW159','20201213-1');
dataPaths{2} = fullfile('AW159','20201213-2');
dataPaths{3} = fullfile('AW162','20210102-1');
dataPaths{4} = fullfile('AW162','20210102-2');
dataPaths{5} = fullfile('AW163','20210102-1');
dataPaths{6} = fullfile('AW163','20210102-2');
dataPaths{7} = fullfile('AW164','20210105');
dataPaths{8} = fullfile('AW165','20210106-1');
dataPaths{9} = fullfile('AW165','20210106-2');

s5 = [];
l5 = [];
for dp = 1:length(dataPaths)
    dp
    
    load(fullfile(saveDir,dataPaths{dp},'Video data','movementData.mat'));
    file = dir(fullfile(saveDir,dataPaths{dp},'StimInfo','*stimInfo.mat'));
    load(fullfile(file(end).folder,file(end).name));

    
    for i = 1:5
        lightEOI{i} = find(ismember(stimInfo.order,12*(i-1) + (1:12)));
        soundEOI{i} = find(ismember(stimInfo.order,12*(i-1+5) + (1:12)));
        lightLocomotion{i} = [];
        soundLocomotion{i} = [];
    end
    for e = 1:length(lightEOI{1})-1
        for i = 1:5
            ev = lightEOI{i}(e);
            frameRange = (-15:1:45) + movementData.eventFrames(ev);
            trial = movementData.frameMovement(frameRange);
            
            baseline = mean(trial(1:15));
            baseSTD = std(trial(1:15));
            normed = (trial-baseline)./baseSTD;
            lightLocomotion{i} = [lightLocomotion{i};normed];
            
            ev = soundEOI{i}(e);
            frameRange = (-15:1:45) + movementData.eventFrames(ev);
            trial = movementData.frameMovement(frameRange);
            
            baseline = mean(trial(1:15));
            baseSTD = std(trial(1:15));
            normed = (trial-baseline)./baseSTD;
            soundLocomotion{i} = [soundLocomotion{i};normed];
        end
    end
    
    l5 = [l5; mean(lightLocomotion{5})];
    s5 = [s5; mean(soundLocomotion{5})];
    
end
time = (-15:1:45)./movementData.frameRate;
figure;hold on;
plot(time,mean(l5));
plot(time,mean(s5));

l5x = [l5(1:5,:); l5(7:8,:)];
s5x = [s5(1:5,:); s5(7:8,:)];

figure;hold on;
vv = mean(l5x); av = mean(s5x);
vvSTD = std(l5x)./sqrt(size(l5x,1));
avSTD = std(s5x)./sqrt(size(s5x,1));
hA = area(time,[vv-vvSTD; 2*vvSTD]')
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
hA = area(time,[av-avSTD; 2*avSTD]')
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 1];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(time,vv,'Color',[0 0 0]);
plot(time,av,'Color',[0 0 1]);



%% Z-scored movement trials histograms

clear;

saveDir = fullfile('/Volumes/AARON DATA/','Electrophysiology','EP010');

dataPaths{1} = fullfile('AW159','20201213-1');
dataPaths{2} = fullfile('AW159','20201213-2');
dataPaths{3} = fullfile('AW162','20210102-1');
dataPaths{4} = fullfile('AW162','20210102-2');
dataPaths{5} = fullfile('AW163','20210102-1');
dataPaths{6} = fullfile('AW163','20210102-2');
dataPaths{7} = fullfile('AW164','20210105');
dataPaths{8} = fullfile('AW165','20210106-1');
dataPaths{9} = fullfile('AW165','20210106-2');

s5 = [];
l5 = [];
for dp = 1:length(dataPaths)
    dp
    
    load(fullfile(saveDir,dataPaths{dp},'Video data','movementData.mat'));
    movementTrials(1,1)
    l = reshape(movementTrials(49:60,:),[1 120]);
    s = reshape(movementTrials(109:120,:),[1 120]);
    
    lNorm = (l - mean(l)) ./ std(l);
    sNorm = (s - mean(l)) ./ std(s);
    
    l5 = [l5; lNorm];
    s5 = [s5; sNorm];
end
figure;
histogram(l5);
hold on;
histogram(s5);

figure;
for r = 1:9
    subplot(3,3,r);
    histogram(l5(r,:));
    hold on;
    histogram(s5(r,:));
end


%% Get numbers for Venn diagram
clear;

saveDir = fullfile('/Volumes/AARON DATA/','Electrophysiology','EP010');

dataPaths{1} = fullfile('AW159','20201213-1');
dataPaths{2} = fullfile('AW159','20201213-2');
dataPaths{3} = fullfile('AW162','20210102-1');
dataPaths{4} = fullfile('AW162','20210102-2');
dataPaths{5} = fullfile('AW163','20210102-1');
dataPaths{6} = fullfile('AW163','20210102-2');
dataPaths{7} = fullfile('AW164','20210105');
dataPaths{8} = fullfile('AW165','20210106-1');
dataPaths{9} = fullfile('AW165','20210106-2');

totalUnits = 0;
light = 0;
soundLight = 0;
motionLight = 0;
soundMotionLight = 0;
sInteract = 0;
mInteract = 0;

for dp = 1:length(dataPaths)
    dp
    
    load(fullfile(saveDir,dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise_Video_final','AVMultiContrastDriftingGratingsWhiteNoiseData.mat'));
    
    for n = 1:length(unitData)
        neuronNumber = unitData(n).neuronNumber;
        
        if ismember(neuronNumber,[responsiveUnits.singleUnits responsiveUnits.multiUnits])
            totalUnits = totalUnits+1;
            
            if ismember(neuronNumber,responsiveUnitsGLM.lightResponsiveUnits)
                light = light+1;
                
                if ismember(neuronNumber,responsiveUnitsGLM.soundResponsiveUnits)
                    soundLight = soundLight+1;
                    if ismember(neuronNumber,responsiveUnitsGLM.lightSoundInteractUnits)
                        sInteract = sInteract+1;
                    end
                end
                if ismember(neuronNumber,responsiveUnitsGLM.movementResponsiveUnits)
                    motionLight = motionLight+1;
                    if ismember(neuronNumber,responsiveUnitsGLM.lightMoveInteractUnits)
                        mInteract = mInteract+1;
                    end
                end
                if ismember(neuronNumber,responsiveUnitsGLM.soundResponsiveUnits) &&...
                        ismember(neuronNumber,responsiveUnitsGLM.movementResponsiveUnits)
                    soundMotionLight = soundMotionLight+1;
                end
            end
        end
    end
end



%% P-value of GLM predictor variables

clear;

saveDir = fullfile('/Volumes/AARON DATA/','Electrophysiology','EP010');

dataPaths{1} = fullfile('AW159','20201213-1');
dataPaths{2} = fullfile('AW159','20201213-2');
dataPaths{3} = fullfile('AW162','20210102-1');
dataPaths{4} = fullfile('AW162','20210102-2');
dataPaths{5} = fullfile('AW163','20210102-1');
dataPaths{6} = fullfile('AW163','20210102-2');
dataPaths{7} = fullfile('AW164','20210105');
dataPaths{8} = fullfile('AW165','20210106-1');
dataPaths{9} = fullfile('AW165','20210106-2');

s = [];
l = [];
m = [];
for dp = 1:length(dataPaths)
    dp
    
    load(fullfile(saveDir,dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise_Video_NEW','AVMultiContrastDriftingGratingsWhiteNoiseData.mat'));
    
    for u = 1:length(unitData)
        neuronNumber = unitData(u).neuronNumber;
        
        if ismember(neuronNumber,responsiveUnitsGLM.lightResponsiveUnits)
            l = [l; unitData(u).glmPval(1)];
        end
        if ismember(neuronNumber,responsiveUnitsGLM.soundResponsiveUnits)
            s = [s; unitData(u).glmPval(2)];
        end
        if ismember(neuronNumber,responsiveUnitsGLM.movementResponsiveUnits)
            m = [m; unitData(u).glmPval(3)];
        end
    end
end
figure;subplot(1,3,1);histogram(floor(log10(l)));
subplot(1,3,2);histogram(floor(log10(s)));
subplot(1,3,3);histogram(floor(log10(m)));
subplot(1,3,2);xlabel('Log10(pval)')
subplot(1,3,1);ylabel('Number of units')


%%  Residuals with and without sound
clear;

saveDir = fullfile('/Volumes/AARON DATA/','Electrophysiology','EP010');

dataPaths{1} = fullfile('AW159','20201213-1');
dataPaths{2} = fullfile('AW159','20201213-2');
dataPaths{3} = fullfile('AW162','20210102-1');
dataPaths{4} = fullfile('AW162','20210102-2');
dataPaths{5} = fullfile('AW163','20210102-1');
dataPaths{6} = fullfile('AW163','20210102-2');
dataPaths{7} = fullfile('AW164','20210105');
dataPaths{8} = fullfile('AW165','20210106-1');
dataPaths{9} = fullfile('AW165','20210106-2');


soundResid = [];
noSoundResid = [];
soundResidS = [];
noSoundResidS = [];
soundResidM = [];
noSoundResidM = [];
soundResidSM = [];
noSoundResidSM = [];
soundResidAll = [];
noSoundResidAll = [];
shuffSoundAll = [];
soundRsquareAll = [];
noSoundRsquareAll = [];
soundRsquareInteract = [];
soundRsquareNoInteract = [];
noSoundRsquareInteract = [];
noSoundRsquareNoInteract = [];
for dp = 1:length(dataPaths)
    dp
    
    load(fullfile(saveDir,dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise_Video_final','AVMultiContrastDriftingGratingsWhiteNoiseData.mat'));
    load(fullfile(saveDir,dataPaths{dp},'Video data','movementData.mat'));
    
    for u = 1:length(unitData)
        neuronNumber = unitData(u).neuronNumber;
        
        if (ismember(neuronNumber,responsiveUnits.singleUnits) ||...
                ismember(neuronNumber,responsiveUnits.multiUnits)) &&...
                    ~(dp==3 && u==59)
            
%             [u neuronNumber]
            glmLight = [zeros(12,10);ones(12,10)*0.25;ones(12,10)*0.5;ones(12,10)*0.75;ones(12,10);...
                zeros(12,10);ones(12,10)*0.25;ones(12,10)*0.5;ones(12,10)*0.75;ones(12,10)];
            glmSound = [zeros(60,10);ones(60,10)];
            
            quantScalar = 1000/(analysisParams.quantWindow(2)-analysisParams.quantWindow(1));
            totalTrials = 1200;
            lightCol = reshape(glmLight,[totalTrials 1]);
            soundCol = reshape(glmSound,[totalTrials 1]);
            locoCol = reshape(movementTrials,[totalTrials 1]);
            responseCol = reshape(unitData(u).trialResponse/quantScalar,[totalTrials 1]);
            
            glmOutput = fitglm([lightCol soundCol locoCol],responseCol,'interactions','distr','poisson','link','log');
%             glmPval = glmOutput.Coefficients.pValue(2:4);
            sResid = mean(abs(glmOutput.Residuals.Raw));
            sResid = glmOutput.Rsquared.Adjusted;
            sRsquare = [glmOutput.Rsquared.Ordinary glmOutput.Rsquared.Adjusted glmOutput.Rsquared.AdjGeneralized];
            soundResidAll = [soundResidAll; sResid];
            soundRsquareAll = [soundRsquareAll; sRsquare];
            
            glmOutput = fitglm([lightCol locoCol],responseCol,'interactions','distr','poisson','link','log');
%             glmPval = glmOutput.Coefficients.pValue(2:3);
            nsResid = mean(abs(glmOutput.Residuals.Raw));
            nsResid = glmOutput.Rsquared.Adjusted;
            nsRsquare = [glmOutput.Rsquared.Ordinary glmOutput.Rsquared.Adjusted glmOutput.Rsquared.AdjGeneralized];
            noSoundResidAll = [noSoundResidAll; nsResid];
            noSoundRsquareAll = [noSoundRsquareAll; nsRsquare];
            
            glmOutput = fitglm([lightCol soundCol(randperm(length(soundCol))) locoCol],responseCol,'interactions','distr','poisson','link','log');
%             glmPval = glmOutput.Coefficients.pValue(2:4);
            xResid = mean(abs(glmOutput.Residuals.Raw));
            shuffSoundAll = [shuffSoundAll; xResid];
            
            if ismember(neuronNumber,responsiveUnitsGLM.lightSoundInteractUnits)
                soundRsquareInteract = [soundRsquareInteract; sResid];
                noSoundRsquareInteract = [noSoundRsquareInteract; nsResid];
            elseif ismember(neuronNumber,responsiveUnitsGLM.soundResponsiveUnits)
                soundRsquareNoInteract = [soundRsquareNoInteract; sResid];
                noSoundRsquareNoInteract = [noSoundRsquareNoInteract; nsResid];
            end
            
            if ismember(neuronNumber,responsiveUnitsGLM.soundResponsiveUnits) &&...
                ismember(neuronNumber,responsiveUnitsGLM.movementResponsiveUnits)
                soundResidSM = [soundResidSM; sResid];
                noSoundResidSM = [noSoundResidSM; nsResid];
            elseif ismember(neuronNumber,responsiveUnitsGLM.soundResponsiveUnits)
                soundResidS = [soundResidS; sResid];
                noSoundResidS = [noSoundResidS; nsResid];
            elseif ismember(neuronNumber,responsiveUnitsGLM.movementResponsiveUnits)
                soundResidM = [soundResidM; sResid];
                noSoundResidM = [noSoundResidM; nsResid];
            else
                soundResid = [soundResid; sResid];
                noSoundResid = [noSoundResid; nsResid];
            end
        end
    end
end
% figure;scatter(soundResidAll,noSoundResidAll)
% xlabel('Sound induded GLM residuals');
% ylabel('Sound excluded GLM residuals');
% 
figure;hold on;
scatter(noSoundResidS,soundResidS);
scatter(noSoundResidSM,soundResidSM);
scatter(noSoundResidM,soundResidM);
ylabel('GLM R-squared (sound included)');
xlabel('GLM R-squared (sound excluded)');
legend({'Sound resp','Sound and motion resp','Motion resp'})

figure;hold on;
histogram(soundResidS - noSoundResidS,-.05:0.05:0.6,'Normalization','probability');
histogram(soundResidSM - noSoundResidSM,-.05:0.05:0.6,'Normalization','probability');
histogram(soundResidM - noSoundResidM,-.05:0.05:0.6,'Normalization','probability');
xlabel('\DeltaR-squared (Sound included - sound excluded)');
ylabel('Number of units');
legend({'Sound resp','Sound and motion resp','Motion resp'})

figure;
subplot(1,2,1);
scatter(soundRsquareAll(:,2),noSoundRsquareAll(:,2));
xlabel('Sound induded GLM R-squared');
ylabel('Sound excluded GLM R-squared');

subplot(1,2,2);
histogram(soundRsquareAll(:,2) - noSoundRsquareAll(:,2))
xlabel('\DeltaR-squared (Sound included - sound excluded)');
ylabel('Number of units');

% figure;
% subplot(1,2,1);hold on;
% scatter(soundRsquareInteract, noSoundRsquareInteract);
% scatter(soundRsquareNoInteract,noSoundRsquareNoInteract);
% xlabel('Sound induded GLM R-squared');
% ylabel('Sound excluded GLM R-squared');
% legend({'Sound-interact resp','Sound alone resp'})
% 
% subplot(1,2,2);hold on;
% histogram(soundRsquareInteract - noSoundRsquareInteract,-.05:0.05:0.6,'Normalization','probability');
% histogram(soundRsquareNoInteract - noSoundRsquareNoInteract,-.05:0.05:0.6,'Normalization','probability');
% xlabel('\DeltaR-squared (Sound included - sound excluded)');
% ylabel('Number of units');
% legend({'Sound-interact resp','Sound alone resp'})

% figure;scatter(soundResidAll,shuffSoundAll);
    
%% GLM-estimated activity with and without sound 
% showing sound is net positive on neuron's firing rate

clear;

saveDir = fullfile('/Volumes/AARON DATA/','Electrophysiology','EP010');

dataPaths{1} = fullfile('AW159','20201213-1');
dataPaths{2} = fullfile('AW159','20201213-2');
dataPaths{3} = fullfile('AW162','20210102-1');
dataPaths{4} = fullfile('AW162','20210102-2');
dataPaths{5} = fullfile('AW163','20210102-1');
dataPaths{6} = fullfile('AW163','20210102-2');
dataPaths{7} = fullfile('AW164','20210105');
dataPaths{8} = fullfile('AW165','20210106-1');
dataPaths{9} = fullfile('AW165','20210106-2');


with = [];
without = [];
for dp = 1:length(dataPaths)
    dp
    
    load(fullfile(saveDir,dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise_Video','AVMultiContrastDriftingGratingsWhiteNoiseData.mat'));
    load(fullfile(saveDir,dataPaths{dp},'Video data','movementData.mat'));
    
    for u = 1:length(unitData)
        [dp u]
        neuronNumber = unitData(u).neuronNumber;
        
        if ismember(neuronNumber,responsiveUnitsGLM.soundResponsiveUnits)
            
            
            glmLight = [zeros(12,10);ones(12,10)*0.25;ones(12,10)*0.5;ones(12,10)*0.75;ones(12,10);...
                zeros(12,10);ones(12,10)*0.25;ones(12,10)*0.5;ones(12,10)*0.75;ones(12,10)];
            glmSound = [zeros(60,10);ones(60,10)];
            
            totalTrials = 1200;
            lightCol = reshape(glmLight,[totalTrials 1]);
            soundCol = reshape(glmSound,[totalTrials 1]);
            locoCol = reshape(movementTrials,[totalTrials 1]);
            responseCol = reshape(unitData(u).trialResponse,[totalTrials 1]);
            
            glmOutput = fitglm([lightCol soundCol locoCol],responseCol,'distr','poisson','link','log');
            
            a = glmOutput.Coefficients.Estimate(2);
            b = glmOutput.Coefficients.Estimate(3);
            c = glmOutput.Coefficients.Estimate(4);
            const = glmOutput.Coefficients.Estimate(1);
            
            thisUnitWith = zeros(1,1200);
            thisUnitWithout = zeros(1,1200);
            for t = 1:1200
                x1 = glmOutput.Variables.x1(t);
                x2 = glmOutput.Variables.x2(t);
                x3 = glmOutput.Variables.x3(t);
                
                lp = const + (a*x1 + b*x2 + c*x3);
                y = exp(lp);
                
                thisUnitWith(t) = y;
                
                lp = const + (a*x1 + c*x3);
                y = exp(lp);
                
                thisUnitWithout(t) = y;
            end
            with = [with; mean(thisUnitWith)];
            without = [without; mean(thisUnitWithout)];
        end
    end
end
figure;scatter(with,without);


%% Sign (positive or negative) of Sound and Sound-light and sound-movement coefficients

clear;

saveDir = fullfile('/Volumes/AARON DATA/','Electrophysiology','EP010');

dataPaths{1} = fullfile('AW159','20201213-1');
dataPaths{2} = fullfile('AW159','20201213-2');
dataPaths{3} = fullfile('AW162','20210102-1');
dataPaths{4} = fullfile('AW162','20210102-2');
dataPaths{5} = fullfile('AW163','20210102-1');
dataPaths{6} = fullfile('AW163','20210102-2');
dataPaths{7} = fullfile('AW164','20210105');
dataPaths{8} = fullfile('AW165','20210106-1');
dataPaths{9} = fullfile('AW165','20210106-2');


soundCoef = [];
soundLightCoef = [];
soundMotionCoef = [];

for dp = 1:length(dataPaths)
    dp
    
    load(fullfile(saveDir,dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise_Video_NEW','AVMultiContrastDriftingGratingsWhiteNoiseData.mat'));
    load(fullfile(saveDir,dataPaths{dp},'Video data','movementData.mat'));
    
    for u = 1:length(unitData)
        neuronNumber = unitData(u).neuronNumber;
        
        if (ismember(neuronNumber,responsiveUnitsGLM.soundResponsiveUnits) ||...
                ismember(neuronNumber,responsiveUnitsGLM.lightSoundInteractUnits) ||...
                ismember(neuronNumber,responsiveUnitsGLM.soundMoveInteractUnits)) &&...
                    ~(dp==3 && u==59)
            
            coeff = unitData(u).glmOutput.Coefficients.Estimate(2:end);
            
            if ismember(neuronNumber,responsiveUnitsGLM.soundResponsiveUnits)
                soundCoef = [soundCoef; coeff(2)];
            end
            
            if ismember(neuronNumber,responsiveUnitsGLM.lightSoundInteractUnits)
                soundLightCoef = [soundLightCoef; coeff(4)];
            end
            
            if ismember(neuronNumber,responsiveUnitsGLM.soundMoveInteractUnits)
                soundMotionCoef = [soundMotionCoef; coeff(6)];
            end
        end
    end
end

figure;
subplot(1,3,1)
histogram(soundCoef);
xlabel('Coefficient of Sound');

subplot(1,3,2);
histogram(soundLightCoef);
xlabel('Coefficient of sound*light');

subplot(1,3,3);
histogram(soundMotionCoef);
xlabel('Coefficient of sound*movement');



%% FR with and without sound

clear;

saveDir = fullfile('/Volumes/AARON DATA/','Electrophysiology','EP010');

dataPaths{1} = fullfile('AW159','20201213-1');
dataPaths{2} = fullfile('AW159','20201213-2');
dataPaths{3} = fullfile('AW162','20210102-1');
dataPaths{4} = fullfile('AW162','20210102-2');
dataPaths{5} = fullfile('AW163','20210102-1');
dataPaths{6} = fullfile('AW163','20210102-2');
dataPaths{7} = fullfile('AW164','20210105');
dataPaths{8} = fullfile('AW165','20210106-1');
dataPaths{9} = fullfile('AW165','20210106-2');

v = [];
ll = [];
a = [];
for dp = 1:length(dataPaths)
    dp
    
    load(fullfile(saveDir,dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise_Video_final','AVMultiContrastDriftingGratingsWhiteNoiseData.mat'));
    load(fullfile(saveDir,dataPaths{dp},'Video data','movementData.mat'));
    
    for u = 1:length(unitData)
        neuronNumber = unitData(u).neuronNumber;
        
        if ismember(neuronNumber,responsiveUnitsGLM.soundResponsiveUnits) &&...
                ismember(neuronNumber,responsiveUnitsGLM.lightResponsiveUnits)
            
            temp = zeros(3,5);
            for c = 1:5
                indices = 12*(c-1) + (1:12);
                li = mean(unitData(u).meanResponse(indices,4));
                so = mean(unitData(u).meanResponse(indices+60,4));
                
                temp(1,c) = li;
                temp(2,c) = so;
                temp(3,c) = so - li;
            end
            v = [v; temp(1,:)];
            a = [a; temp(2,:)];
            ll = [ll; temp(3,:)];
        end
    end
end

figure;plot(0:0.25:1,mean(ll))

figure;hold on;
vv = mean(ll);
vvSTD = std(ll)./sqrt(size(ll,1));
hA = area(0:0.25:1,[vv-vvSTD; 2*vvSTD]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(0:0.25:1,vv,'Color',[0 0 0]);


figure;hold on;
vv = mean(v); av = mean(a);
vvSTD = std(v)./sqrt(size(v,1));
avSTD = std(a)./sqrt(size(a,1));
hA = area(0:0.25:1,[vv-vvSTD; 2*vvSTD]')
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
hA = area(0:0.25:1,[av-avSTD; 2*avSTD]')
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 1];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(0:0.25:1,vv,'Color',[0 0 0]);
plot(0:0.25:1,av,'Color',[0 0 1]);



%% Light and light/sound FR binned by z-scored movement


clear;

saveDir = fullfile('/Volumes/AARON DATA/','Electrophysiology','EP010');

dataPaths{1} = fullfile('AW159','20201213-1');
dataPaths{2} = fullfile('AW159','20201213-2');
dataPaths{3} = fullfile('AW162','20210102-1');
dataPaths{4} = fullfile('AW162','20210102-2');
dataPaths{5} = fullfile('AW163','20210102-1');
dataPaths{6} = fullfile('AW163','20210102-2');
dataPaths{7} = fullfile('AW164','20210105');
dataPaths{8} = fullfile('AW165','20210106-1');
dataPaths{9} = fullfile('AW165','20210106-2');

range = linspace(-1.5,5,14);
ll = cell(1,13);
ss = cell(1,13);
for dp = 1:length(dataPaths)
    dp
    
    load(fullfile(saveDir,dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise_Video_final','AVMultiContrastDriftingGratingsWhiteNoiseData.mat'));
    load(fullfile(saveDir,dataPaths{dp},'Video data','movementData.mat'));
    
    for u = 1:length(unitData)
        neuronNumber = unitData(u).neuronNumber;
        
        if ismember(neuronNumber,responsiveUnitsGLM.lightSoundInteractUnits) &&...
                ismember(neuronNumber,responsiveUnitsGLM.lightResponsiveUnits)
            
            moves = reshape(movementTrials(49:60,:),[1 120]);
            avmoves = reshape(movementTrials(109:120,:),[1 120]);
            both = [moves avmoves];
            Lmoves = (moves-mean(moves)) ./ std(moves);
%             Lmoves = (moves-mean(both)) ./ std(both);
            
            li = reshape(unitData(u).trialResponse(49:60,:),[1 120]);
            for i = 1:13
                TOI = find(Lmoves>range(i) & Lmoves<range(i+1));
                temp = mean(li(TOI));
                if ~isnan(temp)
                    ll{i} = [ll{i}; mean(li(TOI))];
                end
            end
            
            avmoves = reshape(movementTrials(109:120,:),[1 120]);
            Smoves = (avmoves-mean(moves)) ./ std(moves);
%             Smoves = (avmoves-mean(both)) ./ std(both);
            
            so = reshape(unitData(u).trialResponse(109:120,:),[1 120]);
            for i = 1:13
                TOI = find(Smoves>range(i) & Smoves<range(i+1));
                temp = mean(so(TOI));
                if ~isnan(temp)
                    ss{i} = [ss{i}; mean(so(TOI))];
                end
            end
        end
    end
end
for i = 2:13
    lll(i-1) = mean(ll{i});
    lstd(i-1) = std(ll{i})/sqrt(length(ll{i}));
    sss(i-1) = mean(ss{i});
    sstd(i-1) = std(ss{i})/sqrt(length(ss{i}));
end
figure;hold on;
plot(range(2:13)+0.5,lll)
plot(range(2:13)+0.5,sss)

figure;hold on;
hA = area(range(2:13)+0.5,[lll-lstd; 2*lstd]')
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
hA = area(range(2:13)+0.5,[sss-sstd; 2*sstd]')
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 1];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(range(2:13)+0.5,lll,'Color',[0 0 0]);
plot(range(2:13)+0.5,sss,'Color',[0 0 1]);


%% Light and light/sound FR bar plot with stationary, low and high movement


clear;

saveDir = fullfile('/Volumes/AARON DATA/','Electrophysiology','EP010');

dataPaths{1} = fullfile('AW159','20201213-1');
dataPaths{2} = fullfile('AW159','20201213-2');
dataPaths{3} = fullfile('AW162','20210102-1');
dataPaths{4} = fullfile('AW162','20210102-2');
dataPaths{5} = fullfile('AW163','20210102-1');
dataPaths{6} = fullfile('AW163','20210102-2');
dataPaths{7} = fullfile('AW164','20210105');
dataPaths{8} = fullfile('AW165','20210106-1');
dataPaths{9} = fullfile('AW165','20210106-2');

range = [-10 -0.5 1.5 10];
ll = cell(1,3);
ss = cell(1,3);
binCounts = zeros(2,3);
both = [];
for dp = 1:length(dataPaths)
    load(fullfile(saveDir,dataPaths{dp},'Video data','movementData.mat'));
    
    vMoves = reshape(movementTrials(49:60,:),[1 120]);
    avMoves = reshape(movementTrials(109:120,:),[1 120]);
    bothMoves = [vMoves avMoves];
%     zVMoves = (vMoves-mean(vMoves)) ./ std(vMoves);
    zVMoves = (vMoves-mean(bothMoves)) ./ std(bothMoves);
    binCounts(1,1) = binCounts(1,1) + length(find(zVMoves>range(1) & zVMoves<range(2)));
    binCounts(1,2) = binCounts(1,2) + length(find(zVMoves>range(2) & zVMoves<range(3)));
    binCounts(1,3) = binCounts(1,3) + length(find(zVMoves>range(3) & zVMoves<range(4)));

    avMoves = reshape(movementTrials(109:120,:),[1 120]);
%     zAVMoves = (avMoves-mean(avMoves)) ./ std(avMoves);
    zAVMoves = (avMoves-mean(bothMoves)) ./ std(bothMoves);
    binCounts(2,1) = binCounts(2,1) + length(find(zAVMoves>range(1) & zAVMoves<range(2)));
    binCounts(2,2) = binCounts(2,2) + length(find(zAVMoves>range(2) & zAVMoves<range(3)));
    binCounts(2,3) = binCounts(2,3) + length(find(zAVMoves>range(3) & zAVMoves<range(4)));
    both = [both; zVMoves zAVMoves];
end
    
for dp = 1:length(dataPaths)
    dp
    
    load(fullfile(saveDir,dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise_Video_final','AVMultiContrastDriftingGratingsWhiteNoiseData.mat'));
    load(fullfile(saveDir,dataPaths{dp},'Video data','movementData.mat'));
    
    for u = 1:length(unitData)
        neuronNumber = unitData(u).neuronNumber;
        
        if ismember(neuronNumber,responsiveUnitsGLM.soundResponsiveUnits) &&...
                ismember(neuronNumber,responsiveUnitsGLM.lightResponsiveUnits)
%         if ismember(neuronNumber,responsiveUnitsGLM.lightSoundInteractUnits) ||...
%             (ismember(neuronNumber,responsiveUnitsGLM.lightResponsiveUnits) &&...
%             ismember(neuronNumber,responsiveUnitsGLM.soundResponsiveUnits))
        
            moves = reshape(movementTrials(49:60,:),[1 120]);
            avmoves = reshape(movementTrials(109:120,:),[1 120]);
            bothmoves = [moves avmoves];
%             Lmoves = (moves-mean(moves)) ./ std(moves);
            Lmoves = (moves-mean(bothmoves)) ./ std(bothmoves);
            
            li = reshape(unitData(u).trialResponse(49:60,:),[1 120]);
            for i = 1:3
                TOI = find(Lmoves>range(i) & Lmoves<range(i+1));
                temp = mean(li(TOI));
                if ~isnan(temp)
                    ll{i} = [ll{i}; mean(li(TOI))];
                end
            end
            
            avmoves = reshape(movementTrials(109:120,:),[1 120]);
%             Smoves = (avmoves-mean(moves)) ./ std(moves);
            Smoves = (avmoves-mean(bothmoves)) ./ std(bothmoves);
            
            so = reshape(unitData(u).trialResponse(109:120,:),[1 120]);
            for i = 1:3
                TOI = find(Smoves>range(i) & Smoves<range(i+1));
                temp = mean(so(TOI));
                if ~isnan(temp)
                    ss{i} = [ss{i}; mean(so(TOI))];
                end
            end
        end
    end
end
for i = 1:3
    lll(i) = mean(ll{i});
    lstd(i) = std(ll{i})/sqrt(length(ll{i}));
    sss(i) = mean(ss{i});
    sstd(i) = std(ss{i})/sqrt(length(ss{i}));
end
figure;
bar([lll;sss]')
hold on;
ngroups = 3;
nbars = 2;
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
y = [lll;sss]';
err = [lstd;sstd]';
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, y(:,i), err(:,i), '.');
end
hold off


%% Reclassify neurons with interaction terms

clear;

saveDir = fullfile('/Volumes/AARON DATA/','Electrophysiology','EP010');

dataPaths{1} = fullfile('AW159','20201213-1');
dataPaths{2} = fullfile('AW159','20201213-2');
dataPaths{3} = fullfile('AW162','20210102-1');
dataPaths{4} = fullfile('AW162','20210102-2');
dataPaths{5} = fullfile('AW163','20210102-1');
dataPaths{6} = fullfile('AW163','20210102-2');
dataPaths{7} = fullfile('AW164','20210105');
dataPaths{8} = fullfile('AW165','20210106-1');
dataPaths{9} = fullfile('AW165','20210106-2');


for dp = 1:length(dataPaths)
    dp
    
    load(fullfile(saveDir,dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise_Video','AVMultiContrastDriftingGratingsWhiteNoiseData.mat'));
    load(fullfile(saveDir,dataPaths{dp},'Video data','movementData.mat'));
    
    for u = 1:length(unitData)
        [dp u]
        neuronNumber = unitData(u).neuronNumber;
        
        if ismember(neuronNumber,responsiveUnitsGLM.soundResponsiveUnits)
            
            
            glmLight = [zeros(12,10);ones(12,10)*0.25;ones(12,10)*0.5;ones(12,10)*0.75;ones(12,10);...
                zeros(12,10);ones(12,10)*0.25;ones(12,10)*0.5;ones(12,10)*0.75;ones(12,10)];
            glmSound = [zeros(60,10);ones(60,10)];
            
            totalTrials = 1200;
            lightCol = reshape(glmLight,[totalTrials 1]);
            soundCol = reshape(glmSound,[totalTrials 1]);
            locoCol = reshape(movementTrials,[totalTrials 1]);
            responseCol = reshape(unitData(u).trialResponse,[totalTrials 1]);
            
            glmOutput = fitglm([lightCol soundCol locoCol],responseCol,'distr','poisson','link','log');
            
            a = glmOutput.Coefficients.Estimate(2);
            b = glmOutput.Coefficients.Estimate(3);
            c = glmOutput.Coefficients.Estimate(4);
            const = glmOutput.Coefficients.Estimate(1);
            
            thisUnitWith = zeros(1,1200);
            thisUnitWithout = zeros(1,1200);
            for t = 1:1200
                x1 = glmOutput.Variables.x1(t);
                x2 = glmOutput.Variables.x2(t);
                x3 = glmOutput.Variables.x3(t);
                
                lp = const + (a*x1 + b*x2 + c*x3);
                y = exp(lp);
                
                thisUnitWith(t) = y;
                
                lp = const + (a*x1 + c*x3);
                y = exp(lp);
                
                thisUnitWithout(t) = y;
            end
            with = [with; mean(thisUnitWith)];
            without = [without; mean(thisUnitWithout)];
        end
    end
end
figure;scatter(with,without);

