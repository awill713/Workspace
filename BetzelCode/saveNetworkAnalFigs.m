clear
% close all
mouse = {'K048','K056','K070','K073','K074'};
% mouse = {'K048','K056'};

for mm = 3:length(mouse)
    disp(mm)
    x = xlsread('D:\2Pdata\mouseCellTracking.xlsx',mouse{mm});
    r = find(~isnan(x(:,1))); % dates to analyse
    c = x(r,1); % condition
    dates = x(r,3); % dates
    folderNos = x(r,5);
    dataLoc = 'D:\2Pdata\data\';
    
    for ii = 1:length(r)
        disp(ii)
        fns = dir([dataLoc mouse{mm} '_' num2str(dates(ii)) '_2P_FRA_*.mat']);
        
        if ~isempty(fns)
%              return
            fnt = dir([dataLoc mouse{mm} '_' num2str(dates(ii)) '_2P_FRA_' sprintf('%02d',folderNos(ii)) '.mat']);
            load([fnt.folder filesep fnt.name])
            load([fns.folder filesep fns.name],'netanal_tone')
            
            for jj = 1:size(netanal_tone.conCom,2)
                hold on
                subplot(2,size(netanal_tone.conCom,2),jj)
                uMod = unique(netanal_tone.conCom(:,jj));
                nMod = zeros(1,length(uMod));
                cm = colormap('hsv');
                if length(cm)<length(uMod)
                    cm = [cm;colormap('summer')];
                end
                cm = cm(1:floor(length(cm)/length(uMod)):floor(length(cm)/length(uMod))*length(uMod),:);
                
                index = fcn_order_partition(netanal_tone.rho,netanal_tone.conCom(:,jj));
                imagesc(netanal_tone.rho(index,index))
                colormap gray
                title(['# of modules: ' num2str(length(uMod))]);
                hold on
                
                for kk = 1:length(uMod)
                    nMod(kk) = sum(netanal_tone.conCom(:,jj)==uMod(kk));
                    cs = cumsum(nMod(1:kk));
                    if kk==1
                        outline = [0.5,cs(end)+0.5];
                    else
                        outline = [cs(end-1)+0.5,cs(end)+0.5];
                    end
                    patch([outline(1) outline(2) outline(2) outline(1)],[outline(2) outline(2) outline(1) outline(1)],cm(kk,:),'facealpha',0,...
                        'edgeColor',cm(kk,:),'LineWidth',1)
                end
                axis equal tight 
                subplot(2,size(netanal_tone.conCom,2),jj+size(netanal_tone.conCom,2))
               
                imagesc(brighten(spatialInfo.im,1))
                colormap gray
                axis equal tight off
                hold on
                
                for kk = 1:size(netanal_tone.conCom,1)
                    patch(spatialInfo.ROIs{kk}(:,2),spatialInfo.ROIs{kk}(:,1),cm(netanal_tone.conCom(kk,jj),:),'edgealpha',0,'facealpha',0.5)
                end
               
            end
             print(gcf, '-dpdf', '-r400','-bestfit', ['D:\2Pdata\networkFigs\' fnt.name(1:end-4) '_toneEvoked.pdf'])
                clf
        end
    end
end