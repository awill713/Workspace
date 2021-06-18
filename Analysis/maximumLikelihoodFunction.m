function [outputDir] = maximumLikelihoodFunction(distData,probe)
%MAXIMUMLIKELIHOODFUNCTION Summary of this function goes here
%   Detailed explanation goes here

likelihoods = zeros(1,size(distData,2));
for d = 1:length(likelihoods)
%     prob = 1;
%     for nn = 1:size(distData,1)
%         post = pdf(distData(nn,d),probe(nn));
%         if post>0
%             prob = pdf(distData(nn,d),probe(nn)) * prob;
%         end
%         [nn prob];
%     end
%     likelihoods(d) = prob;
    
    prob = 0;
    for nn = 1:size(distData,1)
%         post = pdf(distData(nn,d),probe(nn));
        post = pdf('Poisson',probe(nn),distData(nn,d));
        if post>0
%             prob = log(pdf(distData(nn,d),probe(nn))) + prob;
            prob = log(post) + prob;
        end
        [nn prob];
    end
    likelihoods(d) = prob;
end
[outputProb, outputDir] = max(likelihoods);
% likelihoods
% outputProb

end

