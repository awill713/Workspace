% plot red versus green channel for each cell
% upload a movie not an average image



clear all
close all

mouse='MT032';
date=20190208;
folder_Im2ch='2';
folder_proc='3';

% load(['E:\2Pdata\data\',mouse,'_',num2str(date),'_OptoStim_0',folder_proc,'.mat']);
Stimtype='_FRA_OptoStim_0';
load(['E:\2Pdata\data\',mouse,'_',num2str(date),Stimtype,folder_proc,'.mat']);

%load proc file
tiffolder = '_tifStacks';
procFile=['F_',mouse,'_',num2str(date),mouse,tiffolder,'_plane1_proc.mat'];
load(['E:\2Pdata\',mouse,'\',num2str(date),mouse,tiffolder,'\',folder_proc,'\',procFile]); 
% load non registered image with 2 channels
cd(['E:\2Pdata\',mouse,'\',num2str(date),mouse,'_tifStacks','\',folder_Im2ch])
[file_Im2Ch,path_Im2Ch] = uigetfile('*RED*.tif','Select the red file to open');
Im2chFile{1}=[file_Im2Ch(1:end-7),'GREEN.tif'];
Im2chFile{2}=file_Im2Ch;

info_gr = imfinfo(Im2chFile{1});
info_red = imfinfo(Im2chFile{2});
num_images = numel(info_gr);
for k = 1:num_images
    Imgreen(:,:,k) = imread(Im2chFile{1}, k, 'Info', info_gr);
    Imred(:,:,k) = imread(Im2chFile{2}, k, 'Info', info_red);
    % ... Do something with image A ...
end

Imgreen=double(Imgreen);
Imred=double(Imred);

avg_Imgreen = importdata([path_Im2Ch,'AVG_',Im2chFile{1}]); 
avg_Imred = importdata([path_Im2Ch,'AVG_',Im2chFile{2}]); 

avg_Imgreen=double(avg_Imgreen);
avg_Imred=double(avg_Imred);

regIm=dat.mimg(:,:,2);

% normalise avg_Imgreen and avg_Imred to regIm range of grayscale values
%conversion to dat_img scale : 16-bit image is 0-65535; 15-bit image is
%0-32767 ; it seems from measures im imageJ comparing the scales that by
%just subtracting 32767 we'll get the same scale used in dat.mimg

% Imgreen=(Imgreen-min(min(Imgreen))).*(max(max(regIm))-min(min(regIm)))./(max(max(Imgreen))-min(min(Imgreen)))+min(min(regIm));
% Imred=(Imred-min(min(Imred))).*(max(max(regIm))-min(min(regIm)))./(max(max(Imred))-min(min(Imred)))+min(min(regIm));
% 
% avg_Imgreen=(avg_Imgreen-min(min(avg_Imgreen))).*(max(max(regIm))-min(min(regIm)))./(max(max(avg_Imgreen))-min(min(avg_Imgreen)))+min(min(regIm));
% avg_Imred=(avg_Imred-min(min(avg_Imred))).*(max(max(regIm))-min(min(regIm)))./(max(max(avg_Imred))-min(min(avg_Imred)))+min(min(regIm));

Imgreen = Imgreen - 32767;
Imred = Imred - 32767;
avg_Imgreen = avg_Imgreen - 32767;
avg_Imred = avg_Imred - 32767;

f1=figure;
colormap gray
imb = brighten(double(regIm),1);
imb = brighten(double(imb),0.25);
imagesc(imb)

f2=figure;
colormap gray
im2b = brighten(double(mean(avg_Imgreen,3)),1);
im2b = brighten(double(im2b),0.25);
imagesc(im2b)

f3=figure;
colormap gray
im3b = brighten(double(mean(avg_Imred,3)),1);
im3b = brighten(double(im3b),0.25);
imagesc(im3b)

%% AUTOMATIC VERSION FOR IMAGE ALIGNMENT
% if it doesn't work, do manual alignment

method = 'monomodal';
[optimizer,metric] = imregconfig(method);
% optimizer.MaximumStepLength = 6.25e-10;
% optimizer.MinimumStepLength = 1e-15;
% optimizer.GradientMagnitudeTolerance = 1e-10;
% optimizer.MaximumIterations = 1000;
tform = imregtform(avg_Imgreen,regIm,'rigid',optimizer,metric);
Roriginal = imref2d(size(regIm));
recovered_green = imwarp(Imgreen,tform,'OutputView',Roriginal);
recovered_red = imwarp(Imred,tform,'OutputView',Roriginal);

f4=figure;
imshowpair(mean(recovered_green,3),regIm)
f5=figure;
imshowpair(mean(recovered_red,3),regIm)

%% measure green and red frame per frame


micronsPerPixel=exptInfo.recInfo.micronsPerPixel(1);

img = regIm; %recovered_green;% regIm;

     % create a mask for area outside all cells
     BW_bg=zeros(size(regIm,1),size(regIm,2));
 for kk = 1:size(spatialInfo.ROIs,2)     
    BW = poly2mask(spatialInfo.ROIs{kk}(:,2),spatialInfo.ROIs{kk}(:,1),size(regIm,1),size(regIm,2));
    BW_bg = BW_bg+BW;
 end
    BW_bg = double(BW_bg);
    BW_bg = ones(size(regIm,1),size(regIm,2))-BW_bg;
    BW_bg(BW_bg==0) = NaN; 
 
    filtROI_green_bg = mean(recovered_green,3).*BW_bg;
    filtROI_red_bg = mean(recovered_red,3).*BW_bg;
    mean_g=mean(filtROI_green_bg(~isnan(filtROI_green_bg)));
    std_g=std(filtROI_green_bg(~isnan(filtROI_green_bg)));
    tic
 for kk = 1:size(spatialInfo.ROIs,2)
     % create a mask and convert it to NaN's and 1's
    BW = poly2mask(spatialInfo.ROIs{kk}(:,2),spatialInfo.ROIs{kk}(:,1),size(regIm,1),size(regIm,2));
    BW = double(BW);
    BW(BW==0) = NaN;  
    [a b]=find(BW==1);
     % filter
    filtROI_green = recovered_green(min(a):max(a),min(b):max(b),:).*BW(min(a):max(a),min(b):max(b),:);
    filtROI_red = recovered_red(min(a):max(a),min(b):max(b),:).*BW(min(a):max(a),min(b):max(b),:);
     % cutting matrix to ROI around given cell - to save computing time
    mean_r=mean(filtROI_red_bg(~isnan(filtROI_red_bg)));
    std_r=std(filtROI_red_bg(~isnan(filtROI_red_bg)));
    for ff=1:num_images
        A_g=filtROI_green(:,:,ff);
        A_r=filtROI_red(:,:,ff);
    AvGrayValue_raw_green(kk,ff) = mean(A_g(~isnan(A_g)));
    AvGrayValue_raw_red(kk,ff) = mean(A_r(~isnan(A_r)));         
    AvGrayValue_green(kk,ff) = (mean(A_g(~isnan(A_g)))-mean_g)/std_g;
    AvGrayValue_red(kk,ff) = (mean(A_r(~isnan(A_r)))-mean_r)/std_r;    
    end
    % cell surface in um2
    AvSurf_um2_green(kk)=length(BW(~isnan(BW)))*micronsPerPixel^2;
    clear BW filtROI_green filtROI_red A_g A_r a b
 end
toc

%% plot cell outlines with their cell number

% figure with cell outlines and cell number (=ID)

figure
subplot(121)
imagesc(brighten(brighten(dat.mimg(:,:,2),1),0.5))
colormap gray
hold on
for kk = 1:size(spatialInfo.ROIs,2)
    h = patch(spatialInfo.ROIs{kk}(:,2),spatialInfo.ROIs{kk}(:,1),'y','facealpha',0,'edgealpha',1,'LineWidth',1,'edgecolor',[0.7 0.7 0.7]);
    text(mean(spatialInfo.ROIs{kk}(:,2)),mean(spatialInfo.ROIs{kk}(:,1)),num2str(kk),'Color','m')
end
axis tight equal
set(gca,'TickDir','out','FontSize',14)
xlabel('Dorsal-ventral (pixels)')
ylabel('Rostral-caudal (pixels)')
title({[num2str(date),' - ', mouse],['Green channel'] })

 avg_recovered_red=mean(recovered_red,3);

subplot(122)
imagesc(brighten(brighten(avg_recovered_red,1),0.5))
colormap gray
hold on
for kk = 1:size(spatialInfo.ROIs,2)
    h = patch(spatialInfo.ROIs{kk}(:,2),spatialInfo.ROIs{kk}(:,1),'y','facealpha',0,'edgealpha',1,'LineWidth',1,'edgecolor',[0.7 0.7 0.7]);
    text(mean(spatialInfo.ROIs{kk}(:,2)),mean(spatialInfo.ROIs{kk}(:,1)),num2str(kk),'Color','m')
end
axis tight equal
set(gca,'TickDir','out','FontSize',14)
xlabel('Dorsal-ventral (pixels)')
ylabel('Rostral-caudal (pixels)')
title({[num2str(date),' - ', mouse],['Red channel'] })

%% plot for a given cell the green signal versus red signal
 
 % potential pitfalls given there is no registration used here : if mouse
 % moves, intensity will decrease in both independently of neuronal
 % activity change
 
 test_cell =40;
 
 figure
 subplot(121)
 plot(AvGrayValue_green(test_cell,:),AvGrayValue_red(test_cell,:),'+')
 title(num2str(test_cell))
 
 subplot(122)
  plot(AvGrayValue_green(test_cell,:),'g')
  hold on,
    plot(AvGrayValue_red(test_cell,:),'r')
 title(num2str(test_cell))
 
 %% Find baseline graysscale value in red channel
 
 % strategy is to calculate average red value over timepoints where green is low 
 % what is "low green ?" 
 
 % 1st strategy : look at red values when green is lower than a threshold,
 % say 10 sigmas
 thresh_green=10;
 
 for i=1:size(AvGrayValue_red,1)
 AvGrayValue_grbaseline_red(i) = mean(AvGrayValue_red(i,find(AvGrayValue_green(i,:)<thresh_green)),2);
%  if isnan(AvGrayValue_grbaseline_red(i)) %if there are no green points less than 10 sigmas, take min and calculate red baseline over min green and min+10sigmas
%      AvGrayValue_grbaseline_red(i) = mean(AvGrayValue_red(i,find(AvGrayValue_green(i,:)<min(AvGrayValue_green(i,:))+thresh_green)),2);
%  end
 end
 
 %% 
 
 figure
plot(mean(AvGrayValue_green,2),AvGrayValue_grbaseline_red,'+')
hold on,
hline=refline(1,0);
hline.Color=[0.5 0.5 0.5];
annotation('textbox', [0.2, 0.8, 0.1, 0.1], 'String', {['\sigma = 1 (0.37); n_{red} = ',num2str(length(find(AvGrayValue_grbaseline_red>=1)))]; ...
['\sigma = 2 (0.14); n_{red} = ',num2str(length(find(AvGrayValue_grbaseline_red>=2)))]; ...
['\sigma = 3 (0.05); n_{red} = ',num2str(length(find(AvGrayValue_grbaseline_red>=3)))]})


%% plot "red" cells for different thresholds

%avg_recovered_red=mean(recovered_red,3);

% figure with different thresholds to decide on which cells are red

figure
subplot(231)
imagesc(brighten(brighten(dat.mimg(:,:,2),1),0.5))
colormap gray
hold on
 for kk = 1:size(spatialInfo.ROIs,2)
     if AvGrayValue_grbaseline_red(kk)>=1.5;
         h = patch(spatialInfo.ROIs{kk}(:,2),spatialInfo.ROIs{kk}(:,1),'y','facealpha',0,'edgealpha',1,'LineWidth',2,'edgecolor',[1 0 0]);
     else
         h = patch(spatialInfo.ROIs{kk}(:,2),spatialInfo.ROIs{kk}(:,1),'y','facealpha',0,'edgealpha',1,'LineWidth',1,'edgecolor',[0.7 0.7 0.7]);
     end
 end
axis tight equal
set(gca,'TickDir','out','FontSize',14)
xlabel('Dorsal-ventral (pixels)')
ylabel('Rostral-caudal (pixels)')
title({[num2str(date),' - ', mouse],['red cells : \sigma = 1.5; n = ',num2str(length(find(AvGrayValue_grbaseline_red>=1.5))), '/', num2str(length(AvGrayValue_grbaseline_red))] })


subplot(232)
imagesc(brighten(brighten(dat.mimg(:,:,2),1),0.5))
colormap gray
hold on
 for kk = 1:size(spatialInfo.ROIs,2)
     if AvGrayValue_grbaseline_red(kk)>=2;
         h = patch(spatialInfo.ROIs{kk}(:,2),spatialInfo.ROIs{kk}(:,1),'y','facealpha',0,'edgealpha',1,'LineWidth',2,'edgecolor',[1 0 0]);
     else
         h = patch(spatialInfo.ROIs{kk}(:,2),spatialInfo.ROIs{kk}(:,1),'y','facealpha',0,'edgealpha',1,'LineWidth',1,'edgecolor',[0.7 0.7 0.7]);
     end
 end
axis tight equal
set(gca,'TickDir','out','FontSize',14)
xlabel('Dorsal-ventral (pixels)')
ylabel('Rostral-caudal (pixels)')
title({[num2str(date),' - ', mouse],['red cells : \sigma = 2; n = ',num2str(length(find(AvGrayValue_grbaseline_red>=2))), '/', num2str(length(AvGrayValue_grbaseline_red))] })

subplot(233)
imagesc(brighten(brighten(dat.mimg(:,:,2),1),0.5))
colormap gray
hold on
 for kk = 1:size(spatialInfo.ROIs,2)
     if AvGrayValue_grbaseline_red(kk)>=3;
         h = patch(spatialInfo.ROIs{kk}(:,2),spatialInfo.ROIs{kk}(:,1),'y','facealpha',0,'edgealpha',1,'LineWidth',2,'edgecolor',[1 0 0]);
     else
         h = patch(spatialInfo.ROIs{kk}(:,2),spatialInfo.ROIs{kk}(:,1),'y','facealpha',0,'edgealpha',1,'LineWidth',1,'edgecolor',[0.7 0.7 0.7]);
     end
 end
axis tight equal
set(gca,'TickDir','out','FontSize',14)
xlabel('Dorsal-ventral (pixels)')
ylabel('Rostral-caudal (pixels)')
title({[num2str(date),' - ', mouse],['red cells : \sigma = 3; n = ',num2str(length(find(AvGrayValue_grbaseline_red>=3))), '/', num2str(length(AvGrayValue_grbaseline_red))] })

% red channel

subplot(234)
imagesc(brighten(brighten(avg_recovered_red,1),0.5))
colormap gray
hold on
 for kk = 1:size(spatialInfo.ROIs,2)
     if AvGrayValue_grbaseline_red(kk)>=1.5;
         h = patch(spatialInfo.ROIs{kk}(:,2),spatialInfo.ROIs{kk}(:,1),'y','facealpha',0,'edgealpha',1,'LineWidth',2,'edgecolor',[1 0 0]);
     else
         h = patch(spatialInfo.ROIs{kk}(:,2),spatialInfo.ROIs{kk}(:,1),'y','facealpha',0,'edgealpha',1,'LineWidth',1,'edgecolor',[0.7 0.7 0.7]);
     end
 end
axis tight equal
set(gca,'TickDir','out','FontSize',14)
xlabel('Dorsal-ventral (pixels)')
ylabel('Rostral-caudal (pixels)')


subplot(235)
imagesc(brighten(brighten(avg_recovered_red,1),0.5))
colormap gray
hold on
 for kk = 1:size(spatialInfo.ROIs,2)
     if AvGrayValue_grbaseline_red(kk)>=2;
         h = patch(spatialInfo.ROIs{kk}(:,2),spatialInfo.ROIs{kk}(:,1),'y','facealpha',0,'edgealpha',1,'LineWidth',2,'edgecolor',[1 0 0]);
     else
         h = patch(spatialInfo.ROIs{kk}(:,2),spatialInfo.ROIs{kk}(:,1),'y','facealpha',0,'edgealpha',1,'LineWidth',1,'edgecolor',[0.7 0.7 0.7]);
     end
 end
axis tight equal
set(gca,'TickDir','out','FontSize',14)
xlabel('Dorsal-ventral (pixels)')
ylabel('Rostral-caudal (pixels)')

subplot(236)
imagesc(brighten(brighten(avg_recovered_red,1),0.5))
colormap gray
hold on
 for kk = 1:size(spatialInfo.ROIs,2)
     if AvGrayValue_grbaseline_red(kk)>=3;
         h = patch(spatialInfo.ROIs{kk}(:,2),spatialInfo.ROIs{kk}(:,1),'y','facealpha',0,'edgealpha',1,'LineWidth',2,'edgecolor',[1 0 0]);
     else
         h = patch(spatialInfo.ROIs{kk}(:,2),spatialInfo.ROIs{kk}(:,1),'y','facealpha',0,'edgealpha',1,'LineWidth',1,'edgecolor',[0.7 0.7 0.7]);
     end
 end
axis tight equal
set(gca,'TickDir','out','FontSize',14)
xlabel('Dorsal-ventral (pixels)')
ylabel('Rostral-caudal (pixels)')




%% DECIDE ON THRESHOLD AND SAVE VECTORS

Threshold = 2;

cd(path_Im2Ch)


RedCells.Threshold=Threshold;
RedCells.AvGrayValue_red=AvGrayValue_grbaseline_red;
isCell_red=zeros(size(AvGrayValue_grbaseline_red));
isCell_red(AvGrayValue_grbaseline_red>=Threshold)=1;
RedCells.isCell_red=isCell_red;
RedCells.SourceData=file_Im2Ch;
RedCells.AnalysisFile='IdentificationRedCells_Movie_MT.m';
RedCells.Thresh_green=thresh_green;

save workspace
save('RedCells','RedCells')