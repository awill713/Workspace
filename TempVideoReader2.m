
clear
folder = '/Volumes/AARON DATA/Electrophysiology/EP010/AW165/20210106-2/Video data/';


%% Get frame-wise movement

display('Counting total frames:');
v = VideoReader(fullfile(folder,'output.mp4'));
totalFrames = 0;
while hasFrame(v)
    totalFrames = totalFrames+1;
    f = readFrame(v);
    if mod(totalFrames,5000)==0
        display([num2str(totalFrames) ' frames and counting...']);
    end
end
display([num2str(totalFrames) ' frames counted. Now processing...']);


clear v f;
v = VideoReader(fullfile(folder,'output.mp4'));
frameRate = v.frameRate;
frameMovement = zeros(1,totalFrames);
firstFrame = readFrame(v);
for f = 2:totalFrames
    secondFrame = readFrame(v);
    
    first = im2double(rgb2gray(firstFrame));
    second = im2double(rgb2gray(secondFrame));
    dd = second - first;
    frameMovement(f) = mean(abs(dd(:)));
    
    if mod(f,2000)==0
        display([num2str(f) ' frames processed (' num2str(f/totalFrames*100) '%)']);
    end
    firstFrame = secondFrame;
end
display('Processing complete');
figure;plot((1:totalFrames)/(frameRate*60),frameMovement);
xlabel('Minutes');
ylabel('Frame-wise movement');

%% Import event times

fileID = fopen(fullfile(folder,'eventTimes.txt'),'r');
eventTimes = fscanf(fileID, '%f' );
eventFrames = round(eventTimes.*frameRate);
% locomotion = [];
% for e = 1:length(eventFrames)
%     frameRange = (-15:1:45) + eventFrames(e);
%     trial = frameMovement(frameRange);
%     locomotion = [locomotion; trial];
% end
% meanLoco = mean(locomotion);
% figure;plot(-15:1:45,meanLoco);


%% Save

movementData.frameMovement = frameMovement;
movementData.eventTimes = eventTimes;
movementData.eventFrames = eventFrames;
movementData.frameRate = frameRate;
movementData.totalFrames = totalFrames;
save(fullfile(folder,'movementData.mat'),'movementData');



%% Make video of locomotion graph to accompany desired video clip

vidObj = VideoWriter(fullfile(folder,'writtenVideo.avi'));
open(vidObj);

startFrame = round((7*60)*v.frameRate);
endFrame = round(8*60*v.frameRate);
vidFrames = endFrame - startFrame + 1;
vidTime = vidFrames/v.frameRate;
figure;
xlabel('Time (s)');
ylabel('Frame-wise movement');
ymax = 0.01;
for ff = startFrame:endFrame
    ff - startFrame
    nowTime = (ff-startFrame)/v.frameRate;
    if ff-startFrame < 10*v.frameRate
        frame1 = startFrame;
        time1 = 0;
    else
        frame1 = round(ff - 10*v.frameRate)+1;
        time1 = nowTime - 10;
    end
%     time1 = ff - 30*v.FrameRate;
    plot(time1:(1/v.frameRate):nowTime, frameMovement(frame1:ff));
    if (frameMovement(ff)+0.005) > ymax
        ymax = frameMovement(ff)+0.005;
    end
    xlim([nowTime-10 nowTime]);
    xticks(0:5:vidTime);
    xlabel('Time (s)');
    ylim([0 ymax]);
    ylabel('Frame-wise motion');
    if time1 ~= nowTime
        xlim([time1 nowTime]);
    end
    
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
end
close(vidObj);


%%

for i = 1:5
    lightEOI{i} = find(ismember(stimInfo.order,12*(i-1) + (1:12)));
    soundEOI{i} = find(ismember(stimInfo.order,12*(i-1+5) + (1:12)));
    lightLocomotion{i} = [];
    soundLocomotion{i} = [];
end
for e = 1:length(lightEOI{1})
    for i = 1:5
        ev = lightEOI{i}(e);
        frameRange = (-15:1:45) + movementData.eventFrames(ev);
        trial = movementData.frameMovement(frameRange);
        lightLocomotion{i} = [lightLocomotion{i};trial];
        
        ev = soundEOI{i}(e);
        frameRange = (-15:1:45) + movementData.eventFrames(ev);
        trial = movementData.frameMovement(frameRange);
        soundLocomotion{i} = [soundLocomotion{i};trial];
    end
end
figure;
for i = 1:5
    hold on;
    plot(-15:1:45,mean(lightLocomotion{i}))
end
figure;
for i = 1:5
    hold on;
    plot(-15:1:45,mean(soundLocomotion{i}))
end


soundeoi = find(ismember(stimInfo.order,109:120));
for e = 1:length(soundeoi)
     ev = soundeoi{i}(e);
        frameRange = (-15:1:45) + movementData.eventFrames(ev);
        trial = movementData.frameMovement(frameRange);
        sL = [sL;trial];
end
meanSound = mean(sL);
figure;hold on;
plot(-15:1:45,meanSound);


% locomotion = [];
% for e = 1:length(eventFrames)
%     frameRange = (-15:1:45) + eventFrames(e);
%     trial = frameMovement(frameRange);
%     locomotion = [locomotion; trial];
% end
% meanLoco = mean(locomotion);
% figure;plot(-15:1:45,meanLoco);
