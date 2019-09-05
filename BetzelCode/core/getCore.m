clear clear all
close all
clc
%% load consensus partitions
load('cicon');
[N,T] = size(cicon);
%% generate some null partitions -- same number/size comdnities
d = agreement(cicon);
nrand = 100;
dr = zeros(N,N,nrand);
for irand = 1:nrand
    [~,r] = sort(rand(size(cicon)));
    r = bsxfun(@plus,r,(0:(T - 1))*N);
    ciconr = cicon(r);
    dr(:,:,irand) = agreement(ciconr);
end
z = (d - mean(dr,3))./std(dr,[],3);
%% make a figure
% show that some nodes/cells are co-assigned to the same module more often
% than expected.
f = figure;
hist(z(triu(ones(N),1) > 0),11);
xlabel('z-score co-assignment probability');
ylabel('count');
%% find core/periphery:
% method adapted from "Task-based core-periphery organization of human 
% brain dynamics".

% the goal is to assign nodes to a core and a periphery. usually, these
% assignments are binary -- a node is or is not part of core/periphery. in
% the bassett paper they parameterized this process so that 1) the size of
% the core is variable, and 2) the distinction between core and periphery
% can be sharp/soft/anywhere in between. there are two parameters that
% control these features.

% we don't know the values of those parameters a priori, so we test a
% range. we still need a good way of choosing the parameters, but that's a
% bridge i'll cross next week.

m = 5;                 % number of steps between parameter values
niter = 2;              % the algorithm is non-deterministic so you need to run a few times, in practice
X = cell(m,m,niter);    % save the node "coreness" scores
R = zeros(m,m,niter);   % save the quality of the core (note: sign needs to be flipped)
for i = 1:m
    for j = 1:m
        for iter = 1:niter
            [X{i,j,iter},R(i,j,iter)] = ...
                corepercont_productnorm(d,i/m,j/m);
            imagesc(nanmean(R,3)); drawnow;
        end
    end
end
Rmean = nanmean(R,3);       % mean quality across iterations
mm = min(Rmean(:));         % find them minimum value (top quality)
[u,v] = find(Rmean == mm);  % find parameters where that occurred
[~,idx] = min(R(u,v,:));    % of those, find the smallest
x = X{u,v,idx};             % these are node "coreness" scores. larger value are more core-like, smaller values more periphery-like. you              

