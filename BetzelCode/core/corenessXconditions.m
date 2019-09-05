% Look at how core cells are across the conditions, do they become more or
% less core? Make a matrix for each cell with core values

clear

mouse = 'K070';


x = xlsread('D:\2Pdata\mouseCellTracking.xlsx',mouse);
r = find(x(:,2)==1); % dates to analyse
c = x(r,1); % condition
dates = x(r,3); % dates
uC = unique(c); % unique conditions
fLabels = {'pre','first','second'}; % file labels
fLabels = fLabels(uC+1);
cellIndex = cell(1,length(uC));
Xm = cellIndex;

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
    
    X = squeeze(mean(X,4)); % average across repeats(10)
    
    alphas = size(X,2);
    betas = size(X,3);
    gammas = size(X,4);
    
    for aa = 1:alphas
        for bb = 1:betas
            Xm{ii}(:,:,aa,bb) = squeeze(X(:,aa,bb,:)); % average across repeats, and paramaeter values alpha and beta
        end
    end
    
    %     topcells_prc = prctile(Xm,80); % Look at top 20% of cells
    %     topcells = Xm>repmat(topcells_prc,size(Xm,1),1);
    
    % Load session data
    
    for jj = 1:length(ds)
        fnSession = dir(['D:\2Pdata\data\' mouse '_' ds{jj} '_2P_FRA_*.mat']);
        load([fnSession.folder filesep fnSession.name],'proc','exptInfo','stimInfo','proc','spatialInfo')
        cellIndex{ii}(:,jj) = spatialInfo.masterCellIndex(cs(:,jj),2);
    end
    
end


%% find cells in all conditions
cells=[];
for ii = 1:length(cellIndex)
    cells = [cells;cellIndex{ii}(:,1)];
end
uCells = unique(cells);

c2p = zeros(length(uCells),length(cellIndex));
for jj = 1:length(cellIndex)
    c2p(:,jj) = ismember(uCells,cellIndex{jj}(:,1));
end

trackedCells = uCells(sum(c2p,2)==size(c2p,2));

% now extract data for those cells
for gg = 1:gammas
    for aa = 1:alphas
        for bb = 1:betas
            index = 1; coreness=[];
            for ii = 1:length(Xm)
                coreness(:,:,index) = Xm{ii}(ismember(cellIndex{ii}(:,1),trackedCells),:,aa,bb);
                index = index+1;
            end
            
            
            plot(squeeze(coreness(:,gg,:))','k')
            colormap gray
            title(['alpha: ' num2str(aa) ', beta: ' num2str(bb) ', gamma: ' num2str(gg)])
            drawnow
            pause(0.05)
        end
    end
end








