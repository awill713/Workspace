

clear;
recFolder = uigetdir('/Users/Aaron/Documents/','Select folder with unprocessed data');

%% Load stuff

% vids = noResource(dir([recFolder '/*.tif']));
vids = dir([recFolder '/*.tif']);

infoFile = dir([recFolder '/*exptInfo.mat']);
load([recFolder '/' infoFile.name]);


%% Determine event order

f = exptInfo.framesOut;
a = f>0.2;
da = diff(a);
frameTime = find(da==1)+1;
fr = 1/mean(diff(frameTime)/exptInfo.sampleRate);

f = exptInfo.eventsOut;
a = f>0.67;
da = diff(a);
eventTimes = find(da==1)+1;
eventTimes(1) = [];
eventTimes(1) = [];
eventsOn = round(eventTimes/exptInfo.sampleRate*fr);

%% Analysis parameters
stimTypes = unique(exptInfo.order);

preToneTime = 1000; %ms
preToneFrames = round(fr*preToneTime/1000);
postToneTime = 2000; %ms
postToneFrames = ceil(fr*postToneTime/1000);
totalFrames = postToneFrames + preToneFrames + 1;
repeats = exptInfo.repeats;


%% Read video
video = zeros(512,512,0);
tstart = tic;
for v = 1:length(vids)
    file = [recFolder '/' vids(v).name]
    frames = size(imfinfo(file),1);
    temp = zeros(512,512,frames);
    
    for f = 1:frames
        temp(:,:,f) = imread(file,f);
    end
    
%     video = cat(3, video, temp);
    video{v} = temp;
    tstop = toc(tstart)
end

%% Analyze

for v = 1:length(video)
    fluor{v} = zeros(1,length(video{v}));
    for f = 1:length(fluor{v})
        fluor{v}(f) = mean(mean(video{v}(:,:,f)));
    end
    v
end

ultimate = zeros(1,0);
for v = 1:length(fluor)
    ultimate = cat(2,ultimate,fluor{v});
end

% fluor = zeros(1,length(video));
% for f = 1:length(video)
%     fluor(f) = mean(mean(video(:,:,f)));
% end

response = zeros(length(stimTypes),repeats,totalFrames);
meanResponse = zeros(length(stimTypes),totalFrames);
for type = 1:length(stimTypes)
    eventsOfType = find(exptInfo.order == type);
    for e = 1:length(eventsOfType)
        t = eventsOn(eventsOfType(e));
        response(type,e,:) = ultimate(t-preToneFrames:t+postToneFrames);
    end
    meanResponse(type,:) = squeeze(mean(response(type,:,:),2));
end

figure;hold on;
plot(meanResponse(1,:));
plot(meanResponse(2,:));
plot(meanResponse(3,:));
legend({'Sound','Light','Both'});

figure;
plot(squeeze(response(1,:,:))');
title('Sound');

figure;
plot(squeeze(response(2,:,:))');
title('Light');

figure;
plot(squeeze(response(3,:,:))');
title('Both');



% function directory = noResource(directory)
% directory(contains({directory.name},'._')) = [];
% end