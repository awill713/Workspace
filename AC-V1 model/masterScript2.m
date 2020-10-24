
totalUnits = 100;

repeats = 10;

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
tuningCurves = zeros(totalUnits,length(orientations));

for n = 1:orientationSelectiveUnits
    shift = randsample(0:length(orientations)-1,1);
    temp = circshift(defaultTuningCurve,shift);
    prefOrientations(n) = orientations(mod(4+shift-1,length(orientations))+1);
    tuningCurves(n,:) = temp;
end

%%
stimResponses = cell(2,length(contrasts));
for c = 1:length(contrasts)
    responses = zeros(totalUnits,length(orientations),repeats);
    
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
                    recentBin = max([1 t-10]);
                    acActivity = mean(sum(acNeurons(:,recentBin:t),2));
                    acInput = acActivity/50;
                    
                    in = 0.5*vStim(t)*contrasts(c)*tuningCurves(n,o);
                    prob = 0.2*exp(in*10-5)./(exp(in*10-5)+1);
                    
                    firstBin = max([1 t-500]);
                    recentSpikes = sum(v1Neurons(n,firstBin:t))/20;
                    refract = 0.75*exp(-10*(recentSpikes-1))./(exp(-10*(recentSpikes-1))+1)+0.25;
                    
                    prob = prob*refract;
                    v1Neurons(n,t) = randsample([0 1],1,true,[1-prob prob]);
                end
                
                %         for n = 1:size(av1Neurons,1)
                %             recentBin = max([1 t-10]);
                %             acActivity = mean(sum(acNeurons(:,recentBin:t),2));
                %             acInput = acActivity/20;
                %
                %             in = 0.5*vStim(t)*contrasts(c) + acInput;
                %             prob = 0.2*exp(in*10-5)./(exp(in*10-5)+1);
                %
                %             firstBin = max([1 t-500]);
                %             recentSpikes = sum(av1Neurons(n,firstBin:t))/20;
                %             refract = 0.75*exp(-10*(recentSpikes-1))./(exp(-10*(recentSpikes-1))+1)+0.25;
                %
                %             prob = prob*refract;
                %             av1Neurons(n,t) = randsample([0 1],1,true,[1-prob prob]);
                %         end
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
            
            responses(:,o,r) = mean(v1FR(:,1000:1300),2);
        end
    end
end

    
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