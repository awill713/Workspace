function [Z,med,medsubsample] = fcn_choose_gamma(ci,subsample,nsubsample)
[~,nreps,ngam] = size(ci);
mask = triu(ones(nreps),1) > 0;         % upper triangle mask
masksub = triu(ones(subsample),1) > 0;  % a different upper triangle mask
med = zeros(1,ngam);                    % vector to store median similarity over ALL cluster partitions
medsubsample = zeros(nsubsample,ngam);  % matrix to store median similarity from subsamples
Z = zeros(nreps,nreps,ngam);
for igam = 1:ngam                       % loop over gamma
    fprintf('gamma %i/%i',igam,ngam);
    tic;
    z = zeros(nreps);                   % we'll temporarily store similarity scores here.
    for irep = 1:(nreps - 1)            % loop over pairs of partitions
        for jrep = (irep + 1):nreps
            z(irep,jrep) = fcn_zrand(ci(:,irep,igam),ci(:,jrep,igam));  % compute similarity
            z(jrep,irep) = z(irep,jrep);                                % z-score rand index is symmetric, so z(i,j) is the same as z(j,i)
        end
    end
    med(igam) = median(z(mask));        % get median similarity
    for isub = 1:nsubsample             % randomly subsample
        r = randperm(nreps,subsample);  % select pairs at random
        zsub = z(r,r);                  % extract subsampled matrix
        medsubsample(isub,igam) = median(zsub(masksub));    % compute median z-score rand index
    end
    Z(:,:,igam) = z;
    fprintf(' ... %.2f s\n',toc);
end

pks = zeros(size(medsubsample));
for isub = 1:nsubsample                 % loop over subsamples and find peaks
    idx = findpeaks(medsubsample(isub,:));
    pks(isub,idx) = 1;
end
mu = mean(pks,1) >= 0.95;               % keep peaks that are consistent in 95% of subsamples
cipks = ci(:,:,mu);
cicon = zeros(N,size(cipks,3));         % do consensus clustering to get "representative" clusters (lancichinetti & fortunato 2012, sci. rep.)
for j = 1:size(cipks,3)
    cicon(:,j) = fcn_consensus_communities(cipks(:,:,j),nreps,true); %  these represent the find clusters
end

function ciu = fcn_consensus_communities(ci,niter,vis)
if ~exist('vis','var')
    vis = false;
end
N = size(ci,1);
mask = triu(ones(N),1) > 0;
d = agreement(ci);
goFlag = length(unique(d));
if goFlag <= 2
    CiCon = ci;
end
while goFlag > 2
    
    mu = mean(d(mask));
    b = d - mu;
    b(1:(N + 1):end) = 0;
    CiCon = zeros(N,niter);
    for iRep = 1:niter
        CiCon(:,iRep) = genlouvain(b);
        if vis
            imagesc(CiCon); drawnow;
        end
    end
    d = agreement(CiCon);
    goFlag = length(unique(d));
    
end
ciu = CiCon(:,1);

function D = agreement(ci,buffsz)
%AGREEMENT      agreement matrix from clusters
%
%   D = AGREEMENT(CI) takes as input a set of vertex partitions CI of
%   dimensions [vertex x partition]. Each column in CI contains the
%   assignments of each vertex to a class/community/module. This function
%   aggregates the partitions in CI into a square [vertex x vertex]
%   agreement matrix D, whose elements indicate the number of times any two
%   vertices were assigned to the same class.
%
%   In the case that the number of nodes and partitions in CI is large
%   (greater than ~1000 nodes or greater than ~1000 partitions), the script
%   can be made faster by computing D in pieces. The optional input BUFFSZ
%   determines the size of each piece. Trial and error has found that
%   BUFFSZ ~ 150 works well.
%
%   Inputs,     CI,     set of (possibly) degenerate partitions
%               BUFFSZ, optional second argument to set buffer size
%
%   Outputs:    D,      agreement matrix
%
%   Richard Betzel, Indiana University, 2012
%

%modification history
%09.24.2012 - added loop for big N that makes the function slower but also
% prevents it from maxing out memory.

n = size(ci,2);

if nargin < 2
    buffsz = 1000;
end

if n <= buffsz
    
    ind = dummyvar(ci);
    D = ind*ind';
    
else
    
    a = 1:buffsz:n;
    b = buffsz:buffsz:n;
    
    if length(a) ~= length(b)
        b = [b, n];
    end
    
    x = [a' b'];
    nbuff = size(x,1);
    
    D = zeros(size(ci,1));
    for i = 1:nbuff
       y = ci(:,x(i,1):x(i,2));
       ind = dummyvar(y);
       D = D + ind*ind';
    end
    
end

D = D.*~eye(length(D));


function xmax=findpeaks(data,threshold)
% Helper function to find peaks in a given continuous valued time series x
% Usage: xmax=findpeaks(data,threshold)
% Input:
%      data     (data in time x channels/trials form or a single vector)
%      threshold (if specified returns locations of peaks at which data exceeds threshold) - optional
% Output:
%      xmax     (locations of local maxima of data in a structure array of dimensions channels/trials)
if nargin < 1; error('Need data'); end;
data=change_row_to_column(data);
C=size(data,2);
pp1=[data(1,:);data(1:end-1,:)];
pp2=[data(2:end,:);data(end,:)];
xmax(1:C)=struct('loc',[]);
% for ch=1:C,
%    if nargin ==1
%      xmax(ch).loc=[xmax(ch).loc; find(data(:,ch)-pp1(:,ch)>0 & data(:,ch)-pp2(:,ch)>0)];
%    else
%      xmax(ch).loc=[xmax(ch).loc; find(data(:,ch)-pp1(:,ch)>0 & data(:,ch)-pp2(:,ch)>0 & data(:,ch)>threshold)];
%    end
% end

for ch=1:C,
   if nargin ==1
     xmax(ch).loc=[xmax(ch).loc; find(data(:,ch)-pp1(:,ch)>=0 & data(:,ch)-pp2(:,ch)>=0)];
   else
     xmax(ch).loc=[xmax(ch).loc; find(data(:,ch)-pp1(:,ch)>=0 & data(:,ch)-pp2(:,ch)>=0 & data(:,ch)>=threshold)];
   end
end

function zr = fcn_zrand(x,y)
% FCN_ZRAND         z-score rand index
%
%   Reference: Traud et al (2010). Comparing community structure to
%   characteristics in online collegiate social networks. arxiv:0809.0690v3
%



% clear all
% close all
% clc

% two partitions
% x = [1 1 1 2 2 2 3 3 3 4 4 5]';
% y = [1 1 2 2 3 3 4 4 5 5 6 6]';

% build agreement matrices
indx = dummyvar(x);
indy = dummyvar(y);

% get size of partitions
n = length(x);

% build consistency matrix
nxy = indx'*indy;
vxy = nxy(nxy >= 2);
w11 = 0;
for i = 1:length(vxy)
    w11 = w11 + vxy(i)*(vxy(i) - 1)/2;
end

% row and column sum of consistency matrix
nx = sum(nxy,2);
ny = sum(nxy);
M1 = sum(nx.*(nx - 1)/2);
M2 = sum(ny.*(ny - 1)/2);
M = n*(n - 1)/2;

c1 = n*(n.^2 - 3*n - 2) - 8*(n + 1)*M1 + 4*sum(nx.^3);
c2 = n*(n.^2 - 3*n - 2) - 8*(n + 1)*M2 + 4*sum(ny.^3);

a = M/16;
b = (4*M1 - 2*M).^2;
c = (4*M2 - 2*M).^2;
d = c1*c2/(16*n*(n - 1)*(n - 2));
e = (b - 4*c1 - 4*M)*(c - 4*c2 - 4*M)./(64*n*(n - 1)*(n - 2)*(n - 3));

s = a - (b*c)/(256*(M^2)) + d + e;
zr = (sqrt(s)^-1)*(w11 - (M1*M2/M));