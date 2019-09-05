
%Correlate intensities of red and green stains to see if it's bleedthrough

redImage = imread('/Volumes/AARON FILES/Retrograde/Retro 4 (antero)/Confocal/New/AW083_s2r2s6_IC_10x_green/AW083_s2r2s6_IC_10x_stackMaxRed.tif');
greenImage = imread('/Volumes/AARON FILES/Retrograde/Retro 4 (antero)/Confocal/New/AW083_s2r2s6_IC_10x_green/AW083_s2r2s6_IC_10x_stackMaxGreen.tif');
figure;scatter(reshape(redImage,[1024*1024 1]),reshape(greenImage,[1024*1024 1]))


redControl = rgb2gray(imread('/Volumes/AARON FILES/IHC/IHC 5/PV-488 9-6-18/Confocal/AW046_AC_s1r2s3/AW046_AC_s1r2s3/AW046_AC_s1r2s3_z05c2.tif'));
greenControl = rgb2gray(imread('/Volumes/AARON FILES/IHC/IHC 5/PV-488 9-6-18/Confocal/AW046_AC_s1r2s3/AW046_AC_s1r2s3/AW046_AC_s1r2s3_z05c3.tif'));
figure;scatter(reshape(redControl,[512*512 1]),reshape(greenControl,[512*512 1]))

redNew = rgb2gray(imread('/Volumes/AARON FILES/Retrograde/Retro 4 (antero)/Confocal/Newest/AW083_s2r2s6_IC_20x/AW083_s2r2s6_IC_20x/AW083_s2r2s6_IC_20x_redMax.tif'));
greenNew = rgb2gray(imread('/Volumes/AARON FILES/Retrograde/Retro 4 (antero)/Confocal/Newest/AW083_s2r2s6_IC_20x/AW083_s2r2s6_IC_20x/AW083_s2r2s6_IC_20x_greenMax.tif'));
figure;scatter(reshape(redNew,[512*512 1]),reshape(greenNew,[512*512 1]))

threshold = 113;
threshRed = redNew;
threshRed(threshRed<threshold) = 0;
figure;imshow(redNew)
figure;imshow(threshRed)

threshold = 80;
threshGreen = greenNew;
threshGreen(threshGreen<threshold) = 0;
figure;imshow(greenNew)
figure;imshow(threshGreen)
