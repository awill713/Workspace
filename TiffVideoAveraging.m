

file_path = '/Volumes/AARON FILES/Two photon/2PM 001/AW017/2P 20180519/20180519_AW017_pureTones-001_01.tif';
vid = zeros(512,512,2000);
for frame = 1:2000
    vid(:,:,frame) = imread(file_path,frame);
end

disp('Done');