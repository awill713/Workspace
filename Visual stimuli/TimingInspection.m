

clear;
close all;

tic
soundThresh = 1e-7;
soundFiltThresh = 0.01;
timeThresh = [30,30,40,50,20];
lightThresh = 1.55e-3;


master.v0.voltage = load('D:\Two photon\2PM 003\AW000\AV 2P 20181005\v0data.mat');
master.v0.info = load('D:\Two photon\2PM 003\AW000\AV 2P 20181005\20181005_AW000_AVtwophoton_exptInfo_v0.mat');

master.v1.voltage = load('D:\Two photon\2PM 003\AW000\AV 2P 20181005\v1data.mat');
master.v1.info = load('D:\Two photon\2PM 003\AW000\AV 2P 20181005\20181005_AW000_AVtwophoton_exptInfo_v1.mat');

master.v2.voltage = load('D:\Two photon\2PM 003\AW000\AV 2P 20181005\v2data.mat');
master.v2.info = load('D:\Two photon\2PM 003\AW000\AV 2P 20181005\20181005_AW000_AVtwophoton_exptInfo_v2.mat');

master.v3.voltage = load('D:\Two photon\2PM 003\AW000\AV 2P 20181005\v3data.mat');
master.v3.info = load('D:\Two photon\2PM 003\AW000\AV 2P 20181005\20181005_AW000_AVtwophoton_exptInfo_v3.mat');

master.v5.voltage = load('D:\Two photon\2PM 003\AW000\AV 2P 20181005\v5data.mat');
master.v5.info = load('D:\Two photon\2PM 003\AW000\AV 2P 20181005\20181005_AW000_AVtwophoton_exptInfo_v5.mat');

fields = fieldnames(master);

for i = 1:numel(fields)
    i
    data = master.(fields{i}).voltage.data;
    stimInfo = master.(fields{i}).info.exptInfo;
    
%     [soundSpec freq t] = pspectrum(data(:,1),200e3,'spectrogram','TimeResolution',0.01,'FrequencyLimits',[0 10e3]);
%     
%     freqTime = soundSpec(round(end/2),:);
%     a = freqTime>soundThresh;
%     da = diff(a);
%     up = intersect(find(da==1)+1,find(t>timeThresh(i)));
%     down = intersect(find(da==-1)+1,find(t>timeThresh(i)));
%     timeUp = t(up)';
%     timeDown = t(down)';
%     master.(fields{i}).soundUpSpec = timeUp;
%     master.(fields{i}).soundDownSpec = timeDown;
%     master.(fields{i}).soundDurationSpec = timeDown - timeUp;
    
    sound = data(:,1);
    soundFilt = highpass(sound,1000,200e3);
    a = soundFilt>soundFiltThresh;
    da = diff(a);
    up = find(da==1)+1;
    dUp = diff(up);
    firstUp = find(dUp>(4*200e3))+1;
    frameUp = up(firstUp)';
    timeUp = [up(1) frameUp]/200e3;
    timeUp = timeUp(find(timeUp>timeThresh(i)));
    
    down = find(da==-1)+1;
    dDown = diff(down);
    lastDown = find(dDown>(4*200e3));
    frameDown = down(lastDown)';
    timeDown = [frameDown down(end)]/200e3;
    timeDown = timeDown(find(timeDown>timeThresh(i)));
    
    master.(fields{i}).soundUpRaw = timeUp;
    master.(fields{i}).soundDownRaw = timeDown;
    master.(fields{i}).soundDurationRaw = timeDown - timeUp;
    
    
    light = data(:,2);
    lightFilt = lowpass(light,0.001,200e3);
    
    a = lightFilt>lightThresh;
    da = diff(a);
    up = find(da==1)+1;
    dUp = diff(up);
    firstUp = find(dUp>(4*200e3))+1;
    frameUp = up(firstUp)';
    timeUp = [up(1) frameUp]/200e3;
    timeUp = timeUp(find(timeUp>timeThresh(i)));
    
    down = find(da==-1)+1;
    dDown = diff(down);
    lastDown = find(dDown>(4*200e3));
    frameDown = down(lastDown)';
    timeDown = [frameDown down(end)]/200e3;
    timeDown = timeDown(find(timeDown>timeThresh(i)));
    
    master.(fields{i}).lightUp = timeUp;
    master.(fields{i}).lightDown = timeDown;
    master.(fields{i}).lightDuration = timeDown- timeUp;
end

disp('Done');
toc

master.v0.soundDurationRaw = master.v0.soundDurationRaw(1:end-1);
master.v0.soundUpRaw = master.v0.soundUpRaw(1:end-1);
master.v0.soundDownRaw = master.v0.soundDownRaw(1:end-1);

%% Plot

for i = 1:numel(fields)
    figure;
    
    subplot(3,2,1);
    hist(master.(fields{i}).soundDurationRaw*1000);
    theRange = range(master.(fields{i}).soundDurationRaw*1000);
    avg = mean(master.(fields{i}).soundDurationRaw*1000);
    xlabel('Duration (ms)');
    title([fields{i} ' sound duration, mean =  ' num2str(avg) 'ms, range = ' num2str(theRange) 'ms']);
    
    subplot(3,2,2);
    hist(master.(fields{i}).lightDuration(1:end-1)*1000);
    theRange = range(master.(fields{i}).lightDuration(1:end-1)*1000);
    avg = mean(master.(fields{i}).lightDuration(1:end-1)*1000);
    xlabel('Duration (ms)');
    title([fields{i} ' light duration, mean =  ' num2str(avg) 'ms, range = ' num2str(theRange) 'ms']);
    
    both = master.(fields{i}).info.exptInfo.order;
    soundBoth = find(both(union(find(both==1),find(both==3)))==3);
    lightBoth = find(both(union(find(both==2),find(both==3)))==3);
    
    onDelay = master.(fields{i}).lightUp(lightBoth) - master.(fields{i}).soundUpRaw(soundBoth);
    subplot(3,2,3);
    hist(onDelay*1000);
    theRange = range(onDelay*1000);
    avg = mean(onDelay*1000);
    xlabel('Delay (ms)');
    title([fields{i} ' light onset delay, mean = ' num2str(avg) 'ms, range = ' num2str(theRange) 'ms']);

    offDelay = master.(fields{i}).lightDown(lightBoth) - master.(fields{i}).soundDownRaw(soundBoth);
    subplot(3,2,4);
    hist(offDelay*1000);
    theRange = range(offDelay*1000);
    avg = mean(offDelay*1000);
    xlabel('Delay (ms)');
    title([fields{i} ' light offset delay, mean = ' num2str(avg) 'ms, range = ' num2str(theRange) 'ms']);
    
    subplot(3,2,5);
    scatter(onDelay*1000,offDelay*1000);
    xlabel('Onset delay (ms)');
    ylabel('Offset delay (ms)');
    title([fields{i} ' light onset vs light offset delays']);
end

