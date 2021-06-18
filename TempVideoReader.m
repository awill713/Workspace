
clear;

%Figure out how many frames the video has
display('Counting total frames:');
v = VideoReader('/Volumes/AARON USB 2/AW157-2020-12-12-2/output.mp4');
totalFrames = 0;
while hasFrame(v)
    totalFrames = totalFrames+1;
    f = readFrame(v);
    if mod(totalFrames,5000)==0
        display([num2str(totalFrames) ' frames...']);
    end
end
display([num2str(totalFrames) ' frames counted. Now processing...']);


%Calculate the frame by frame difference
clear v f;
v = VideoReader('/Volumes/AARON USB 2/AW157-2020-12-12-2/output.mp4');
frameRate = v.frameRate;
frameDiff = zeros(1,totalFrames-1);
firstFrame = readFrame(v);
for f = 1:totalFrames-1
    secondFrame = readFrame(v);
    
    first = im2double(rgb2gray(firstFrame));
    second = im2double(rgb2gray(secondFrame));
    dd = second - first;
    frameDiff(f) = mean(dd(:));
    
    if mod(f,2000)==0
        display([num2str(f) ' frames processed (' num2str(f/totalFrames*100) '%)']);
    end
    firstFrame = secondFrame;
end
display('Processing complete');
figure;plot((2:totalFrames)/frameRate,frameDiff); %plot motion over the whole stimulus, just to see 

%Show the average motion for all events, just a sanity check
%I chose this frame range (-15:1:45) because my stimuli were 1 second (30
%frames), so I looked 0.5 seconds before and 0.5 seconds after. You should
%do whatever is appropriate for you. Again, this is just a quick check
%that the video processing worked and to see if there was motion.
fileID = fopen('/Volumes/AARON USB 2/AW157-2020-12-12-1/eventTimes.txt','r');
eventTimes = fscanf(fileID, '%f' );
eventFrames = round(eventTimes.*frameRate);
locomotion = [];
for e = 1:length(eventFrames)
    frameRange = (-15:1:45) + eventFrames(e); 
    trial = frameDiff(frameRange);
    locomotion = [locomotion; trial];
end
meanLoco = mean(locomotion);
figure;plot(-15:1:45,meanLoco);
