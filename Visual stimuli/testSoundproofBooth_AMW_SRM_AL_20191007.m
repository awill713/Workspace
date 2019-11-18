
close all;
clear all;

daqreset;

%% Design stuff
dbSteps = [-20 -15 -10 -5 0 5 10 15 20];
inputChannel = 1;
outputChannel = 0;
sampleRate = 400e3;
duration = 5; %seconds

ref_PA = 20e-6;
volts_per_PA = .316;
[fb, fa] = butter(5, 2*300 / sampleRate, 'high');
%load the filter
filterPath = 'D:\GitHub\filters\20191008_ephysBoothSpkr_300-70k_fs400k_TESTACTUAL';
load(filterPath);

%% Set up nidaq stuff
session = daq.createSession('ni');
session.Rate = sampleRate;
addAnalogOutputChannel(session,'dev3',outputChannel,'Voltage');
addAnalogInputChannel(session,'dev3',inputChannel,'Voltage');
session.Channels(2).InputType = 'SingleEnded';

%% Make noise and filter it
for db = 1:length(dbSteps)
%     intensity = 70+dbSteps(db); %WRONG, changed intensity to attenuation
    attenuation = dbSteps(db);
    noise(db,:) = randn(duration*fs,1).*10^(attenuation/20);
%     noise(db,:) = noise(db,:) - mean(noise(db,:));
    noise(db,:) = conv(noise(db,:), FILT, 'same');
    noise(db,:) = envelopeKCW(noise(db,:),5,sampleRate);
%     figure;plot(noise(db,:));
    temp = noise(db,:);
    temp(find(temp>=10))=9.9999;
    temp(find(temp<=-10)) = -9.9999;
    noise(db,:) = temp;
    noise(db,:) = zeros(duration*fs,1);
%     figure;plot(noise(db,:))
end

%% Cycle through dbSteps and play the noise from speaker, recording via the microphone
for db = 1:length(dbSteps)
    dbSteps(db)
    session.queueOutputData(noise(db,:)');
    [S(db,:), time(db,:)] = session.startForeground();
%     [toneS toneT] = session.startForeground();
    [fb, fa] = butter(5, 2*300 / sampleRate, 'high');
    S(db,:) = filter(fb,fa,S(db,:));
    [power(db,:), freq(db,:)] = pwelch(S(db,:)/ref_PA/volts_per_PA, 1024,120,[],sampleRate,'onesided');
end

%% Plot power spectrum
figure;hold on;
for db = 1:length(dbSteps)
    plot(freq(db,:),10*log10(power(db,:)));
    
    p = power(db,:); f = freq(db,:);
    totalVolume(db) = 10*log10(mean(p) * (f(end)-f(1)));
end
xlabel('Frequency')
ylabel('Power (dB)');
legend(num2str((dbSteps+70)'))
title('Power spectrum (speaker outside of booth, microphone inside, 10/8/2019)');

figure;
scatter(dbSteps+70,totalVolume);
xlabel('SPL of Gaussian noise from speaker (dB)');
ylabel('SPL picked up by microphone (dB)');
title('Total volume (speaker outside of booth, microphone inside, 10/8/2019)');



%% Test tones
% disp('Acquiring 3s of background noise:');
% stim = zeros(1,3*sampleRate);
% b = getResponse_sess(stim,1,session);
% b = filtfilt(fb, fa, b) / ref_PA / volts_per_PA;
% b = b - mean(b);
% noise_ms = mean(b.^2);
% 
% tone = conv(cos(15000*(1:duration*sampleRate)*2*pi/sampleRate),FILT,'same');
% tone = envelopeKCW(tone,5,sampleRate);
% 
% session.queueOutputData(tone');
% [S time] = session.startForeground();
% pause(3)
% session.queueOutputData(tone');
% [S time] = session.startForeground();
% S = S - mean(S);
% 
% butterS = filtfilt(fb,fa,S);
% butterS = butterS(fs*0.005 : end-(fs*0.005));
% [tonePower toneFreq] = pwelch(butterS/ref_PA/volts_per_PA, 1024,120,[],sampleRate,'onesided');
% 
% figure;
% plot(toneFreq, 10*log10(tonePower));
% xlabel('Frequency')
% ylabel('Power');
% 
% 10*log10(mean(tonePower) * (toneFreq(end)-toneFreq(1)))
% 20*log10(rms(butterS/ref_PA/volts_per_PA))
% 
% 20*log10(sqrt(mean((butterS/ref_PA/volts_per_PA).^2 ) - noise_ms))