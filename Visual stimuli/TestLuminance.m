
period = 1; %seconds
totalTime = 10; %seconds
fps = 120;

totalFrames = totalTime*fps;
periodFrames = period*fps;

video = zeros(512,512,totalFrames,3);
image = zeros(512,512,3);

for f = 1:totalFrames
    coef = 2*pi/fps/period;
    point = 0.5*sin(coef*f)+0.5;
    image(:,:,:) = point
    
    video(:,:,f,:) = point;
%     video(:,:,f,2) = 0;
%     video(:,:,f,3) = 112;
end

figure;
for i = 1:600
%     frame = video(:,:,i,:);
    imshow(squeeze(video(:,:,i,:)));
%     i
    pause(0.0083333);
end



sampleRate = 400e3;
freq = 440;

totalSamples = sampleRate*totalTime;

sound = zeros(1,totalSamples);

for i = 1:totalSamples
    fcoef = freq*2*pi/sampleRate;
    ecoef = 2*pi/sampleRate/period;
    sound(i) = cos(fcoef*i) * (0.5*sin(ecoef*i)+0.5);
end

