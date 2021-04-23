

x = 0:0.01:1;
y = 0.2*exp(x*10-5)./(exp(x*10-5)+1);
yy = 0.2*exp(x*10-5)./(exp(x*10-5)+1);
% figure;
% plot(x,y);
% title('Probability of firing');
% 


% x1 = 0.2;
% x2 = 0.3;
% 
% 
% y1 = y(find(x==x1));
% y2 = y(find(x==x2));
% std1 = yy(find(x==x1));
% std2 = yy(find(x==x2));
% 
% dprime = (y2 - y1) / sqrt( std1^2 + std2^2);
% 
% indx = 0:0.01:0.9;
% primes = zeros(1,length(indx));
% for i = 1:length(primes)
%     x1 = indx(i);
%     x2 = x1+0.2;
%     
%     y1 = 0.2*exp(x1*10-5)./(exp(x1*10-5)+1);
%     y2 = 0.2*exp(x2*10-5)./(exp(x2*10-5)+1);
%     yy1 = 0.1*exp(x1*10-5)./(exp(x1*10-5)+1);
%     yy2 = 0.1*exp(x2*10-5)./(exp(x2*10-5)+1);
%     
%     dprime = (y2 - y1) / sqrt( 0.5*(yy1^2 + yy2^2));
%     primes(i) = dprime;
% end
% 
% figure;plot(indx,primes);
% 
