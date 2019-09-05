clear
% close all
mouse = {'K048','K056','K070','K073','K074'};
% mouse = {'K048','K056'};

for mm = 1:length(mouse)
    disp(mm)
    x = xlsread('D:\2Pdata\mouseCellTracking.xlsx',mouse{mm});
    r = find(~isnan(x(:,1))); % dates to analyse
    c = x(r,1); % condition
    dates = x(r,3); % dates
    dataLoc = 'D:\2Pdata\data\';
    
    for ii = 1:length(r)
        fn = dir([dataLoc mouse{mm} '_' num2str(dates(ii)) '_2P_FRA_*.mat']);
        
        if ~isempty(fn)
            load([fn.folder filesep fn.name],'netanal_tone')
            if ~exist('netanal_tone','var')
                netanal_tone = networkAssignmentCode_fcn([fn.folder filesep fn.name]); 
                save([fn.folder filesep fn.name],'netanal_tone','-append')
            end
            clear netanal_tone
        end
    end
end
