clear

mouse = 'K074';


x = xlsread('D:\2Pdata\mouseCellTracking.xlsx',mouse);
r = find(x(:,2)==1); % dates to analyse
c = x(r,1); % condition
dates = x(r,3); % dates
uC = unique(c); % unique conditions
fLabels = {'pre','first','second'}; % file labels

for ii = 1:length(uC)
    
    fnCore = [cd filesep mouse filesep 'coreness.' mouse '.toneEvokedOffset.' fLabels{ii} '.mat'];
    load(fnCore)
    
    % Check dates match up
    dsn = zeros(size(ds)); for jj = 1:length(ds); dsn(jj) = str2double(ds{jj}); end
    dateRows = find(ismember(dates,dsn));
    if length(dateRows)~=length(dsn) % break if can't get same sessions
        disp('sessions do not align!!')
        return
    end
    
    Xm = squeeze(mean(squeeze(mean(squeeze(mean(X,4)),2)),2)); % average across repeats, and paramaeter values alpha and beta
%     Xs = sort(Xm); % sort coreness values
    
    if ii==1
        topcells_prc = prctile(Xm,80); % Look at top 20% of cells
        topcells = Xm>repmat(topcells_prc,size(Xm,1),1);
        
    else
       
    end
    
    % Load session data
    raster = [];
    cellIndex = zeros(size(cs)); %raster = zeros(size(proc.raster,1),length(proc.preEv+1:proc.preEv+floor(exptInfo.fr)),length(cs),size(cs,2));
    for jj = 1:length(ds)
        fnSession = dir(['D:\2Pdata\data\' mouse '_' ds{jj} '_2P_FRA_*.mat']);
        load([fnSession.folder filesep fnSession.name],'proc','exptInfo','stimInfo','proc','spatialInfo')
        cellIndex(:,jj) = spatialInfo.masterCellIndex(cs(:,jj),2);
        ci = find(topcells(:,jj)==1);
        cellsm(:,jj) = spatialInfo.masterCellIndex(ci,2); % master index for tracking same cells
        raster(:,:,:,jj) = proc.raster(:,proc.preEv+1:proc.preEv+floor(exptInfo.fr),cs(:,jj));
        
        order = stimInfo.index(stimInfo.order,:);
        if length(order)<size(raster,1)
            order = repmat(order,ceil(size(raster,1)/length(order)),1);
            order = order(1:size(raster,1),:);
        end
        [~,it] = sort(order(:,1));
        [~,ia] = sort(order(:,2));
        uT = unique(order(:,1));
        uA = unique(order(:,2));
        
        mca = zeros(length(uA),length(uT),size(topcells,2)); % preallocate
        mcap = zeros(length(uA),length(uT),size(topcells,2)); % preallocate
        for gg = 1:size(topcells,2)
            for aa = 1:length(uA)
                for tt = 1:length(uT)
                    rows = find(order(:,1)==uT(tt) & order(:,2)==uA(aa));
                    mca(aa,tt,gg) = mean(mean(nanmean(raster(rows,:,topcells(:,gg),jj),3)));
                    mcap(aa,tt,gg) = mean(mean(nanmean(raster(rows,:,~topcells(:,gg),jj),3)));
                end
            end            
        end
        
        m_resp{ii,jj} = squeeze(mean(mca));
        
    end
        
end
      
index = 1;
for ii  = 1:3
    for jj = 1:3
        subplot(3,3,index)
        plot(mean(m_resp{ii,jj},2))
        index = index+1;
    end
end


