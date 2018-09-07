

teMat = zeros(size(sp{1,2},1));
trials = 1;
trainLength = length(sp{1,2});
neurons = size(sp{1,2},1)
offset = 1;
bestTau = zeros(size(sp{1,2},1));

for i = 1:neurons
    for j = 1:neurons
        if i~=j
            teTau = zeros(1,offset);
            
            teMax = 0;
            maxIndex = 0;
            
            for t = 1:offset
                teTemp = 0;
                randTemp = 0;
                
                for repeat = 1:trials
                    xSpikes = logical(sp{1,2}(i,:));
                    ySpikes = logical(sp{1,2}(j,:));
                    dTE = TransferEntropy(xSpikes,ySpikes,t);
                    teTemp =  teTemp + dTE;
                    
                    xRand = randi([0 1],1,trainLength);
                    yRand = randi([0 1],1,trainLength);
                    dRand = TransferEntropy(xRand,yRand,t);
                    randTemp = randTemp + dRand;
                end
                
                teTemp = teTemp/trials;
                randTemp = randTemp/trials;
                
                teTau(t) = teTemp / randTemp;
            end
            
            [teMax maxIndex] = max(teTau);
            teMat(i,j) = teMax;
            bestTau(i,j) = maxIndex;
        end
    end
    i
end

figure;
imagesc(teMat);
colorbar;

% teValues = zeros(1,555*554/2);
% count = 1;
% for i = 1:555
%     for j = i+1:555
%         teValues(count) = teMat(i,j);
%         count = count+1;
%     end
% end
% figure;
% histogram(teValues);
