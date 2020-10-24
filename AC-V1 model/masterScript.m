
close all;
clear;

dt = 1; %ms

stimDuration = 1000; %ms
totalTime = stimDuration*2;

v1NeuronCount = 1000;
v1TunedPercent = 0.4;

baselineFR = 10; %Hz

v1population = zeros(v1NeuronCount,totalTime);

vStim = [zeros(1,500) ones(1,1000) zeros(1,500)];

% Input to FR prob
x = 0:0.01:1;
y = 0.2*exp(x*10-5)./(exp(x*10-5)+1);
figure;
plot(x,y);
ylim([0 1])

% FR decay
x = 0:0.01:1;
y = exp(-10*x+0.5)./(exp(-10*x+5)+1);
figure;
plot(x,y);


x = 0:0.01:1;
y = 0.75*exp(-10*(x-1))./(exp(-10*(x-1))+1)+0.25;
figure;
plot(x,y);

trials = zeros(1,100);
burst = zeros(1,100);
for tt = 1:length(trials)
    v1Neurons = zeros(1,2000);
    for t = 1:2000
        in = vStim(t);
        prob = 0.2*exp(in*10-5)./(exp(in*10-5)+1);
        v1Neurons(1,t) = randsample([0 1],1,true,[1-prob prob]);
    end
    trials(tt) = sum(v1Neurons(1:500));
    burst(tt) = sum(v1Neurons(500:1000));
end
std(trials) / mean(trials)
std(burst) / mean(burst)

%%
contrasts = [0 0.25 0.5 0.75 1];
vStim = [zeros(1,1000) ones(1,1000) zeros(1,1000)];
aStim = [zeros(1,1000) ones(1,1000) zeros(1,1000)];

totalTime = length(vStim);

meanFR = zeros(2,length(contrasts));
stdFR = zeros(2,length(contrasts));
for c = 1:length(contrasts)
    acNeurons = zeros(1000,totalTime);
    v1Neurons = zeros(1000,totalTime);
    av1Neurons = zeros(1000,totalTime);
    for t = 1:totalTime
        [c t]
        for a = 1:size(acNeurons,1)
            in = aStim(t);
            prob = 0.2*exp(in*10-5)./(exp(in*10-5)+1);
            acNeurons(a,t) = randsample([0 1],1,true,[1-prob prob]);
        end
        
        for n = 1:size(v1Neurons,1)
            recentBin = max([1 t-10]);
            acActivity = mean(sum(acNeurons(:,recentBin:t),2));
            acInput = acActivity/50;
            
            in = 0.5*vStim(t)*contrasts(c);
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
            
            in = 0.5*vStim(t)*contrasts(c) + acInput;
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
end
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


% acFR = zeros(size(acNeurons));
% for t = 1:totalTime
%     for a = 1:size(acNeurons,1)
%         firstBin = max([1 t-9]);
%         spikes = sum(acNeurons(a,firstBin:t));
%         fr = spikes*100;
%         acFR(a,t) = fr;
%     end
% end
% figure;
% plot(mean(acFR));
% 
% v1FR = zeros(size(v1Neurons));
% for t = 1:totalTime
%     for n = 1:size(v1Neurons,1)
%         firstBin = max([1 t-9]);
%         spikes = sum(v1Neurons(n,firstBin:t));
%         fr = spikes*100;
%         v1FR(n,t) = fr;
%     end
% end
% figure;
% plot(mean(v1FR));
% std(mean(v1FR(:,1000:1300),2))/mean(mean(v1FR(:,1000:1300),2))
% std(mean(v1FR(:,3000:3300),2))/mean(mean(v1FR(:,3000:3300),2))




%%

x = linspace(0,1,1000);
m = zeros(size(x));
for t = 1:length(x)
    outcome = randsample([0 1],100,true,[1-x(t) x(t)]);
    m(t) = mean(outcome);
    s(t) = std(outcome);
    v(t) = std(outcome)^2;
end
figure; hold on;
plot(x,m);
plot(x,s);
plot(x,v);