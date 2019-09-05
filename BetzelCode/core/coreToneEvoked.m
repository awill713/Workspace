clear
% close all
mouse = {'K056','K070','K073','K074'};
% mouse = {'K048','K056'};

for mm = 1:length(mouse)
    disp(mm)
    x = xlsread('D:\2Pdata\mouseCellTracking.xlsx',mouse{mm});
    r = find(x(:,2)==1); % dates to analyse
    c = x(r,1); % condition
    dates = x(r,3); % dates
    folderNos = x(r,5);
    dataLoc = 'D:\2Pdata\data\';
    cellN = cell(1,length(r));
    for ii = 1:length(r)
        disp(ii)
        fnt = dir([dataLoc mouse{mm} '_' num2str(dates(ii)) '_2P_FRA_' sprintf('%02d',folderNos(ii)) '.mat']);
        load([fnt.folder filesep fnt.name],'spatialInfo')
        cellN{ii} = spatialInfo.masterCellIndex(:,2);
    end
    
    
    %% Which cells to look at?
    
    cells=[];
    for ii = 1:length(cellN)
        cells = [cells;cellN{ii}];
    end
    uCells = unique(cells);
    cm = zeros(length(uCells),length(r));
    for jj = 1:length(cellN)
        cm(:,jj) = ismember(uCells,cellN{jj});
    end
    
    trackedCells = uCells(sum(cm,2)==length(cellN));
    
    
       for ii = 1:length(r)
           
        fnt = dir([dataLoc mouse{mm} '_' num2str(dates(ii)) '_2P_FRA_' sprintf('%02d',folderNos(ii)) '.mat']);
        load([fnt.folder filesep fnt.name],'netanal_tone','spatialInfo')
        
        
    
    
    
    
    
    
    