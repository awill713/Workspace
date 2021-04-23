
clear;
experiment = 'EP012';
mouseID = 'AW181';
date = '20210406';

pupilFrameRate = 30;

stimPath = fullfile('/Volumes/AARON DATA/Electrophysiology/',experiment,mouseID,date,'StimInfo',[date '_' mouseID '_AVmultiContrastDriftingGratingsWhiteNoise_stimInfo']);
pupilRadiiPath = fullfile('/Volumes/AARON DATA/Electrophysiology/',experiment,mouseID,date,'Pupil video','output_radii.mat');
load(stimPath);
load(pupilRadiiPath);

pupilVidEventTimesPath = fopen(fullfile('/Volumes/AARON DATA/Electrophysiology/',experiment,mouseID,date,'Pupil video','eventTimes.txt'),'r');
eventTimes = fscanf(pupilVidEventTimesPath, '%f' );
eventFrames = round(eventTimes.*30);

pupilSize = interp1(R(:,1),R(:,2),2:R(end,1)); %in case pupilMeasurement() was run with different frameInterval analysis

uniqueEvents = size(stimInfo.index,1);
repeats = stimInfo.repeats;

pupilData = zeros(uniqueEvents,repeats,pupilFrameRate*4+1);
for u = 1:uniqueEvents
    evIdx = u;
    eoi = find(stimInfo.order == evIdx);
    
    for e = 1:length(eoi)
        frameOfInterest = eventFrames(eoi(e));
        firstFrame = frameOfInterest - pupilFrameRate -1; %-1 because the pupil sizes (see variable R) start on second frame for some reason
        lastFrame = frameOfInterest + 3*pupilFrameRate -1;
        
        tempTrace = pupilSize(1,firstFrame:lastFrame);
        
        pupilData(u,e,:) = tempTrace;
    end
end

figure;hold on;
plot(-pupilFrameRate:3*pupilFrameRate,squeeze(nanmean(nanmean(pupilData(1:12,:,:)))))
plot(-pupilFrameRate:3*pupilFrameRate,squeeze(nanmean(nanmean(pupilData(61:72,:,:)))))
legend({'No stimulus','Just sound'});

figure;hold on;
plot(-pupilFrameRate:3*pupilFrameRate,squeeze(nanmean(nanmean(pupilData(49:60,:,:)))))
plot(-pupilFrameRate:3*pupilFrameRate,squeeze(nanmean(nanmean(pupilData(109:120,:,:)))))
legend({'Light','Light + sound'});

save(fullfile('/Volumes/AARON DATA/Electrophysiology/',experiment,mouseID,date,'Pupil video','pupilData.mat'),'pupilData');



%%

traces = zeros(5,2,5,121);

dp=5;
for c = 1:5
        vIdx = (c-1)*12+1 : c*12;
        avIdx = vIdx + 60;
        
        vTemp = squeeze(nanmean(nanmean(pupilData(vIdx,:,:))));
        avTemp = squeeze(nanmean(nanmean(pupilData(avIdx,:,:))));
        
        traces(c,1,dp,:) = vTemp;
        traces(c,2,dp,:) = avTemp;
end

final = zeros(5,2,121);
for c = 1:5
        vTemp = squeeze(mean(traces(c,1,:,:),3));
        avTemp = squeeze(mean(traces(c,2,:,:),3));
        
        vMean = mean(vTemp(1:30));
%         vTemp = (vTemp-vMean)./ vTemp;
        vTemp = (vTemp-vMean)./ vMean;
        
        avMean = mean(avTemp(1:30));
%         avTemp = (avTemp - avMean)./ avTemp;
        avTemp = (avTemp - avMean)./ avMean;
        
        final(c,1,:) = vTemp;
        final(c,2,:) = avTemp;
end