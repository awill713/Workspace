% spatial location of cells within vs across modules

clear

load('D:\2Pdata\communityData.mat')
mouse = {'K056','K070','K073','K074'};

for mm = 2:length(mouse)

    disp(mm)
    x = xlsread('D:\2Pdata\mouseCellTracking.xlsx',mouse{mm});
    recs = x(:,3);
    r = ismember(recs,str2double(data(mm).days)); % dates to analyse
    c = x(r,1); % condition
    dates = x(r,3); % dates
    folderNos = x(r,5);
    dataLoc = 'D:\2Pdata\data\';
    
    if strcmp(mouse{mm},'K070')
        folderNos(ismember(dates,20170731))=[];
        dates(ismember(dates,20170731))=[];
        folderNos(ismember(dates,20170812))=[];
        dates(ismember(dates,20170812))=[];
    end
    
       for ii = 1:length(dates) % for each recording...
        
        disp(ii)
        fnt = dir([dataLoc mouse{mm} '_' num2str(dates(ii)) '_2P_FRA_' sprintf('%02d',folderNos(ii)) '.mat']);
        load([fnt.folder filesep fnt.name])
        
        distMat = zeros(size(proc.raster,3));
        for jj = 1:size(proc.raster,3)
            for kk = 1:size(proc.raster,3)
                x = spatialInfo.centroid(jj,1)-spatialInfo.centroid(kk,1);
                y = spatialInfo.centroid(jj,2)-spatialInfo.centroid(kk,2);
                distMat(jj,kk) = sqrt(x.^2 + y.^2);
            end
        end
        micronsPerPixel = 1.26890590123792;
        distMat = distMat/micronsPerPixel;
        comInd = find(strcmp(data(mm).days,num2str(dates(ii))));
%         drivenCellModules = data(mm).communities{comInd}(sig==1,:); % which modules are the driven neurons in?
        
        % How many of ANY cells are in each module? And how many driven
        % in each module?
        uModules = cell(1,size(data(mm).communities{comInd},2));
        din = zeros(1,length(uModules));
        dacross = din;
        tt_p = din;
        for jj = 1:size(data(mm).communities{comInd},2) % for each gamma value...
            M = false(size(distMat));
            uModules{jj} = unique(data(mm).communities{comInd}(:,jj)); % what are the unique modules?
            for kk = 1:length(uModules{jj}) % for each module...
                cim = find(data(mm).communities{comInd}(:,jj)==kk); % cells in the module
                for zz = 1:length(cim)
                    M(cim(zz),cim) = true;
                end
            end
            M(logical(eye(size(M,1)))) = false;
            M = triu(M);
            din(jj) = mean(distMat(M));
            distsIn = distMat(M);
            histogram(distMat(M),0:5:max(distMat(:)))
            M = triu(~M);
            M(logical(eye(size(M,1)))) = false;
            dacross(jj) = mean(distMat(M));
            distsAcross = distMat(M);
            hold on
            histogram(distMat(M),0:5:max(distMat(:)))
            [~,tt_p(jj)] = ttest2(distsIn,distsAcross);
            title(['ttest p = ' num2str(tt_p(jj))])
            figure
            cs = get(groot,'DefaultAxesColorOrder');
            b = bar(1,mean(distsIn));
            b.FaceColor = cs(5,:);
            b.EdgeAlpha = 0;
            hold on
            b2 = bar(2,mean(distsAcross));
            b2.FaceColor = cs(7,:);
            b2.EdgeAlpha = 0;
            set(gca,'XTick',1:2,'XTickLabel',{'Within cluster','Across cluster'},'FontSize',14)
            ylabel('Mean distance (\mum)')
            errorbar([mean(distsIn) mean(distsAcross)],[mean(distsIn)/sqrt(length(distsIn-1)) mean(distsAcross)/sqrt(length(distsAcross-1))],'LineStyle','none','LineWidth',1,'Color','k')
            pause()
            clf
        end
%         
%         for jj = 1:length(din)
%             plot([1 2],[din(jj) dacross(jj)],'Color',[0.7 0.7 0.7])
%             hold on
%         end
%         b = bar([mean(din),nanmean(dacross)]);
%         b.FaceAlpha = 0.7;
%         b.EdgeColor = 'none';
%         set(gca,'FontSize',14,'XTickLabel',{'within modules','across modules'})
%         ylabel('Mean distance between cells (\mum)')
%         xlim([0.5 2.5])
        
        for jj = 1:size(data(mm).communities{comInd},2) % for each gamma value...
            uModules{jj} = unique(data(mm).communities{comInd}(:,jj)); % what are the unique modules?
            din = []; dout = [];
            for kk = 1:length(uModules{jj}) % for each module...
                cim = find(data(mm).communities{comInd}(:,jj)==kk); % cells in the module
                cnim = find(data(mm).communities{comInd}(:,jj)~=kk); % cells NOT in the module
                
                for zz = 1:length(cim)
                    d = distMat(cim(zz),:); % distance to all other cells
                    dd = d(cim((zz+1):end));% dd(dd==0)=[]; dd(
                    din = [din, dd];
                    dout = [dout, d(cnim)];
                end
                
%                 n{jj}(kk) = sum(data(mm).communities{comInd}(:,jj)==uModules{jj}(kk));
% %                 n_driven{jj}(kk) = sum(data(mm).communities{comInd}(sig==1,jj)==uModules{jj}(kk));
            end
        end
        
        
%         sigCells = find(sig==1);

%         for jj = 1:length(uModules)
%             
%             
% 
% 
%         end
% 
%         axis tight
%         set(gca,'TickDir','out','FontSize',14)
%         xlabel('Dorsal-ventral (\mum)')
%         ylabel('Rostral-caudal (\mum)')
% 
%         
%         % Now perform a shuffle and randomly select the number of tuned
%         % cells from the modules you have and see if any are
%         % significantly different than chance (2.5/97.5%)
%         nShuf = 2000;
%         n_driv_shuf = cell(1,length(uModules));
%         n_shuf_sig = n_driv_shuf;
%         
%         for jj = 1:length(uModules)
%             n_shuf = zeros(nShuf,length(drivenCellModules));
%             for ss = 1:nShuf
%                 s = randperm(sum(n{jj}),length(drivenCellModules));
%                 for kk = 1:length(s)
%                     n_shuf(ss,kk) = find(s(kk) <= cumsum(n{jj}),1,'first');
%                 end
%             end
%             
%             for kk = 1:length(uModules{jj})
%                 n_driv_shuf{jj}(:,kk) = sum(n_shuf == uModules{jj}(kk),2);
%                 n_shuf_sig{jj}(kk) = n_driven{jj}(kk) < prctile(n_driv_shuf{jj}(:,kk),2.5) | n_driven{jj}(kk) > prctile(n_driv_shuf{jj}(:,kk),97.5);
%             end
%         end
%         
%         %         newNet.n_shuf_driv = n_shuf_sig;
%         
%         %         save([fnt.folder filesep fnt.name],'newNet','-append')
%         dc  = find(sig==1);
%         for jj = 1:length(n_shuf_sig)
%             s = find(n_shuf_sig{jj}==1);
%             for ss = 1:length(s)
%                 cells = find(drivenCellModules(:,jj) == s(ss));
%                 disp(cells)
%                 for cc = 1:length(cells)
%                     plot(mean(mresp(dc(cells(cc)),:),1))
%                     hold on
%                 end
%                 pause()
%                 clf
%             end
%         end
        
       end
    
end