
sampleRate = 40e3;

repeats = 5;

toneDuration = 1000;
itiDuration = 1000;

fundFreq = 440;
freqMax = 20000;
numFreq = 10;
overtoneAttn = 0.5;

toneSamples = toneDuration*sampleRate/1000;
ramp(1,:) = logspace(0,-2,toneSamples);
ramp(2,:) = ones(1,toneSamples);
ramp(3,:) = logspace(-2,0,toneSamples);
% linRamp = linspace(0.5,1.5,toneSamples);


itiSamples = itiDuration*sampleRate/1000;
trialSamples = toneSamples+itiSamples;

freqList = logspace(log10(fundFreq),log10(freqMax),numFreq);
% freqList = fundFreq:fundFreq:freqMax;

order = repmat([-1 0 1],[1 repeats]);
order = order(randperm(length(order)));

stimulus = zeros(1,trialSamples*repeats);

for e = 1:length(order)
    type = order(e);
    index = trialSamples*(e-1)+1;
    
    tempStim = zeros(1,toneSamples);
    for f=1:length(freqList)
        coef = freqList(f)*2*pi/sampleRate;
        period = sampleRate/freqList(f);
        phaseOffset = randi(round(period),1);
%         phaseOffset=0;

        freqStim = cos(coef*((1+phaseOffset):(toneSamples+phaseOffset))) * overtoneAttn^(f-1);
%         freqStim = sawtooth(coef*((1+phaseOffset):(toneSamples+phaseOffset))) * overtoneAttn^(f-1);
        
        %     tempStim = cos(coef*(1:toneSamples)) * overtoneAttn^(f-1);
        tempStim = tempStim + freqStim;
    end
    tempStim = rand(1,toneSamples)*2-1;
    
    tempStim = tempStim.*ramp(type+2,:);
    tempStim = tempStim./max(abs(tempStim));
    
    
    stimulus(index:index+toneSamples-1) = tempStim;
end

upTrials = find(order==1);
downTrials = find(order==-1);
staticTrials = find(order==0);
for r = 1:repeats
    upRange(r,:) = trialSamples*(upTrials(r)-1)+1:trialSamples*(upTrials(r)-1)+toneSamples;
    downRange(r,:) = trialSamples*(downTrials(r)-1)+1:trialSamples*(downTrials(r)-1)+toneSamples;
    staticRange(r,:) = trialSamples*(staticTrials(r)-1)+1:trialSamples*(staticTrials(r)-1)+toneSamples;
end
upStim = mean(stimulus(upRange));
downStim = mean(stimulus(downRange));
staticStim = mean(stimulus(staticRange));

figure;
subplot(2,3,1);plot(upStim);subplot(2,3,4);periodogram(upStim,[],[],sampleRate);
subplot(2,3,2);plot(downStim);subplot(2,3,5);periodogram(downStim,[],[],sampleRate);
subplot(2,3,3);plot(staticStim);subplot(2,3,6);periodogram(staticStim,[],[],sampleRate);

sound(stimulus,sampleRate);
