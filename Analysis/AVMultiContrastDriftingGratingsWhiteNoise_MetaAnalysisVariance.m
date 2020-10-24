
clear all;
% close all;

saveDir = fullfile('E:\Electrophysiology','EP004','MetaData - AVMultiContrastDriftingGratingsWhiteNoise');
if ~exist(saveDir)
    mkdir(saveDir);
end

dataPaths{1} = fullfile('EP004','AW117','20200221-1');
dataPaths{2} = fullfile('EP004','AW117','20200221-2');
dataPaths{3} = fullfile('EP004','AW118','20200221-1');
dataPaths{4} = fullfile('EP004','AW118','20200221-2');
dataPaths{5} = fullfile('EP004','AW121','20200226-1');
dataPaths{6} = fullfile('EP004','AW121','20200226-2');
dataPaths{7} = fullfile('EP004','AW124','20200303-1');
dataPaths{8} = fullfile('EP004','AW124','20200303-2');

slopes = zeros(0,5);
yint = zeros(0,5);
fano = zeros(0,5);
fanoC = cell(1,5);

for dp = 1:length(dataPaths)
    dp
    %     [randomize dp]
    %Load data
    dataFile = fullfile('E:\Electrophysiology\',dataPaths{dp},'AVMultiContrastDriftingGratingsWhiteNoise','AVMultiContrastDriftingGratingsWhiteNoiseData.mat');
    load(dataFile);
    stimPath = dir(fullfile('E:\Electrophysiology\',dataPaths{dp},'StimInfo','*AVmultiContrastDriftingGratingsWhiteNoise_stimInfo*'));
    load(fullfile(stimPath.folder,stimPath.name));
    
    orientations = stimInfo.orientations;
    contrasts = stimInfo.contrasts;
    repeats = stimInfo.repeats;
    
    for u = 1:length(unitData)
        neuronNumber = unitData(u).neuronNumber;
        if ismember(neuronNumber,responsiveUnits.soundResponsiveUnits)...
                && ismember(neuronNumber,responsiveUnits.lightResponsiveUnits)...
                && unitData(u).type~=0
            
            meanResponse = unitData(u).meanResponse;
            
            vis = meanResponse(1:length(orientations)*length(contrasts),4)';
            visAud = meanResponse(length(orientations)*length(contrasts)+1:length(orientations)*length(contrasts)*2,4)';
            stdVis = (meanResponse(1:length(orientations)*length(contrasts),5)'*sqrt(repeats));
            stdAudVis = (meanResponse(length(orientations)*length(contrasts)+1:length(orientations)*length(contrasts)*2,5)'*sqrt(repeats));
            pVis = polyfit(vis,stdVis,1);
            pVisAud = polyfit(visAud,stdAudVis,1);
            vF = nanmean(stdVis./vis);
            avF = nanmean(stdAudVis./visAud);
            
            slopes = [slopes; pVis(1) pVisAud(1) dp u neuronNumber];
            yint = [yint; pVis(2) pVisAud(2) dp u neuronNumber];
            fano = [fano; vF avF dp u neuronNumber];
            
            for c = 1:length(contrasts)
                responseVis = meanResponse(((c-1)*length(orientations))+1:((c-1)*length(orientations)+length(orientations)),4);
                responseVisAud = meanResponse(((c-1)*length(orientations)+1+length(contrasts)*length(orientations)):((c-1)*length(orientations)+length(contrasts)*length(orientations)+length(orientations)),4);
                stdVis = (meanResponse(((c-1)*length(orientations))+1:((c-1)*length(orientations)+length(orientations)),5)*sqrt(repeats));
                stdVisAud = (meanResponse(((c-1)*length(orientations)+1+length(contrasts)*length(orientations)):((c-1)*length(orientations)+length(contrasts)*length(orientations)+length(orientations)),5)*sqrt(repeats));
                
                vFc = nanmean(stdVis./responseVis);
                avFc = nanmean(stdVisAud./responseVisAud);
                
                %                     [valV indV]= max(responseVis);
                %                     [valAV indAV] = max(responseVisAud);
                %                     vFc = stdVis(indV)/valV;
                %                     avFc = stdVisAud(indAV)/valAV;
                %
                fanoC{c} = [fanoC{c}; vFc avFc dp u neuronNumber];
                
                %                     pVis = polyfit(responseVis,stdVis,1);
                %                     pVisAud = polyfit(responseVisAud,stdVisAud,1);
                %                     fanoC{c} = [fanoC{c}; pVis(1) pVisAud(1) dp u neuronNumber];
                
            end
        end
    end
end

anovaMat = [];
for c = 1:length(contrasts)
    contFano(1,c) = mean(fanoC{c}(:,1));
    contFano(2,c) = mean(fanoC{c}(:,2));
    stdFano(1,c) = std(fanoC{c}(:,1))/sqrt(size(fanoC{c}(:,1),1));
    stdFano(2,c) = std(fanoC{c}(:,2))/sqrt(size(fanoC{c}(:,2),1));
    
    deltFano(1,c) = mean(fanoC{c}(:,2) - fanoC{c}(:,1));
    stdDeltF(1,c) = std(fanoC{c}(:,2) - fanoC{c}(:,1))/sqrt(size(fanoC{c}(:,2)-fanoC{c}(:,1),1));
    
    anovaMat = [anovaMat; fanoC{c}(:,1) fanoC{c}(:,2)];
end
[p,table,stats] = anova2(anovaMat,size(fanoC{1},1),'off');
figure;hold on;
hA = area(contrasts,[contFano(1,:)-stdFano(1,:); 2*stdFano(1,:)]')
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
hA = area(contrasts,[contFano(2,:)-stdFano(2,:); 2*stdFano(2,:)]')
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 1];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(contrasts,contFano(1,:),'Color',[0 0 0]);
plot(contrasts,contFano(2,:),'Color',[0 0 1]);
xlabel('Contrasts');
ylabel('Fano factor');
% ylim([4 6]);

figure;hold on;
hA = area(contrasts,[deltFano-stdDeltF; 2*stdDeltF]')
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.4;
hA(2).EdgeColor = [1 1 1];
plot(contrasts,deltFano);
xlabel('Contrasts');
ylabel('\Delta Fano factor');
% ylim([-.1, -.04])

f8 = figure;
scatter(fano(:,1),fano(:,2))
dd = mean(fano(:,2)) - mean(fano(:,1));
perc = dd / mean(fano(:,1))*100;
[hh pp] = ttest2(fano(:,1),fano(:,2));
title(['\Delta FF = ' num2str(dd) '(' num2str(perc) ' %)']);

f1 = figure;hold on;
scatter(slopes(:,1),yint(:,1),[],[0 0 0],'filled');
scatter(slopes(:,2),yint(:,2),[],[0 0 1],'filled');
[hslopes pslopes] = ttest(slopes(:,1),slopes(:,2));
[hint pint] = ttest(yint(:,1),yint(:,2));
xlabel('Slope');
ylabel('Y-intercept');
title(['Slopes p = ' num2str(pslopes) ', y-int p = ' num2str(pint)]);
suptitle({'Mean response to variance relationship'});
% saveas(f1,fullfile(saveDir,'Response variance across sound UNresponsive units'));
% save(fullfile(saveDir,'Response variance slopes and intercepts.mat'),'slopes','yint');


% figure; hold on;
% errorbar([mean1 mean2],[sem1 sem2],'Color',[0 0 0],'LineWidth',3);
% for n = 1:size(slopes,1)
%     line([1 2],[slopes(n,1) slopes(n,2)],'Color',[0 0 0 0.1]);
% end

