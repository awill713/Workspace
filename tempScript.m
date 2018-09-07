

% % TransferEntropy(logical(sp{1,1}(1,:)),logical(sp{1,1}(2,:)),1)
% 
% xcorrMat = zeros(size(Fcell{1,1},1));
% zScoreMat = zeros(size(Fcell{1,1},1));
% teMat = zeros(size(Fcell{1,1},1));
% 
% for i = 1:size(Fcell{1,1},1)
%     for j = i:size(Fcell{1,1},1)
%         if i~=j
%             xcorrMat(i,j) = max(xcorr(Fcell{1,1}(i,:),Fcell{1,1}(j,:),'coeff'));
%             zScoreMat(i,j) = max(xcorr(zscore(Fcell{1,1}(i,:)),zscore(Fcell{1,1}(j,:)),'coeff'));
%             teMat(i,j) = TransferEntropy(logical(sp{1,1}(i,:)),logical(sp{1,1}(j,:)),1);
%             
%             xcorrMat(j,i) = max(xcorr(Fcell{1,1}(j,:),Fcell{1,1}(i,:),'coeff'));
%             zScoreMat(j,i) = max(xcorr(zscore(Fcell{1,1}(j,:)),zscore(Fcell{1,1}(i,:)),'coeff'));
%             teMat(j,i) = TransferEntropy(logical(sp{1,1}(j,:)),logical(sp{1,1}(i,:)),1);
%         end
%     end
%     i
% end
% 
% figure;
% imagesc(xcorrMat);
% 
% figure;
% imagesc(zScoreMat);
% 
% count = 1;
% theMetrics = zeros(2,555*554/2);
% raw2t = zeros(2,555*554/2);
% z2t = zeros(2,555*554/2);
% for i = 1:555
%     for j = i+1 :555
%         theMetrics(1,count) = xcorrMat(i,j);
%         theMetrics(2,count) = zScoreMat(i,j);
%         theMetrics(3,count) = teMat(i,j);
%         count = count+1;
%     end
% end
% 
% figure;
% plot(theMetrics(1,:),theMetrics(2,:),'.');
% xlabel('Raw');
% ylabel('Z score');
% 
% figure;
% imagesc(teMat);
% 
% figure;
% plot(theMetrics(1,:),theMetrics(3,:),'.');
% xlabel('Raw');
% ylabel('Transfer entropy');
% 
% figure;
% plot(theMetrics(2,:),theMetrics(3,:),'.');
% xlabel('Z score');
% ylabel('Transfer entropy');


realTE = zeros(1,99);
normTE = zeros(1,99);
trials = 100;

for theTau = 1:99
    tempReal = 0;
    tempNorm = 0;
    
    for repeat = 1:trials
        
        x = zeros(1,100);
        y = zeros(1,100);
        for i = 1:25
            if randi(10)<8
                x(i) = 1;
            end
            if randi(10)>9
                y(i) = 1;
            end
        end
        for i = 26:50
            if randi(10)>9
                x(i) = 1;
            end
            if randi(10)>9
                y(i) = 1;
            end
        end
        for i = 51:75
            if randi(10)>9
                x(i) = 1;
            end
            if randi(10)<8
                y(i) = 1;
            end
        end
        for i = 76:100
            if randi(10)>9
                x(i) = 1;
            end
            if randi(10)>9
                y(i) = 1;
            end
        end
        
        tempReal = TransferEntropy(x,y,theTau) + tempReal;
        
        randX = randi([0 1],1,100);
        randY = randi([0 1],1,100);
        
        tempNorm = TransferEntropy(randX,randY,theTau) + tempNorm;
    end
    tempReal = tempReal / trials;
    tempNorm = tempNorm / trials;
    
    realTE(theTau) = tempReal;
    normTE(theTau) = tempNorm;
end

scaledTE = realTE./normTE;

figure;
plot(realTE);
figure;
plot(normTE);
figure;
plot(scaledTE);
% 
% teValues = zeros(1,99);
% trials = 100;
% for tau = 1:99
%     tempTE = 0;
%     for repeat = 1:trials
%         x = zeros(1,100);
%         y = zeros(1,100);
% %         x = randi([0 1],1,100);
% %         y = randi([0 1],1,100);
% %         y = zeros(1,100);
% %         y(51:100) = x(1:50);
% %         y(51:75) = 1;
%         for i = 1:25
%             if randi(10)<8
%                 x(i) = 1;
%             end
%             if randi(10)>9
%                 y(i) = 1;
%             end
%         end
%         for i = 26:50
%             if randi(10)>9
%                 x(i) = 1;
%             end
%             if randi(10)>9
%                 y(i) = 1;
%             end
%         end
%         for i = 51:75
%             if randi(10)>9
%                 x(i) = 1;
%             end
%             if randi(10)<8
%                 y(i) = 1;
%             end
%         end
%         for i = 76:100
%             if randi(10)>9
%                 x(i) = 1;
%             end
%             if randi(10)>9
%                 y(i) = 1;
%             end
%         end
%         
%         % x = [0 0 1 1 0 0 0 1 0 0];
%         % y = [0 0 0 1 1 0 0 0 1 0];
%         
%         tMinusTau = length(x)-tau;
%         
%         jCounts = zeros(2);
%         probs = zeros(2);
%         for t = 1:tMinusTau
%             t0 = y(t)+1;
%             t1 = y(t+tau)+1;
%             jCounts(t0,t1) = jCounts(t0,t1) + 1;
%         end
%         for c = 1:size(jCounts,2)
%             for r = 1:size(jCounts,1)
%                 probs(r,c) = jCounts(r,c) / sum(jCounts(r,:));
%             end
%         end
%         ijCounts = zeros(2,2,2);
%         ijProbs = zeros(2,2,2);
%         for t = 1:tMinusTau
%             i0 = x(t)+1;
%             j0 = y(t)+1;
%             j1 = y(t+tau)+1;
%             ijCounts(i0,j0,j1) = ijCounts(i0,j0,j1) + 1;
%         end
%         for z = 1:size(ijCounts,3)
%             for c = 1:size(ijCounts,2)
%                 for r = 1:size(ijCounts,1)
%                     ijProbs(r,c,z) = ijCounts(r,c,z) / sum(ijCounts(r,c,:));
%                 end
%             end
%         end
%         
%         TE = 0;
%         
%         for gI = 1:size(ijProbs,1)
%             for gJ = 1:size(ijProbs,2)
%                 for r = 1:size(ijProbs,3)
%                     gI;
%                     gJ;
%                     r;
%                     overallProb = ijCounts(gI,gJ,r)/sum(sum(sum(ijCounts)));
%                     if overallProb == 0
%                         dTE = 0;
%                     else
%                         dTE = overallProb * log2(ijProbs(gI,gJ,r)/probs(gJ,r));
%                     end
%                     TE = TE + dTE;
%                 end
%             end
%         end
%         tempTE = tempTE + TE;
%     end
%     tempTE = tempTE / trials;
%     teValues(tau) = tempTE;
% end
% figure;
% plot(teValues);
