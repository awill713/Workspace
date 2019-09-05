

dim = 1000;
square = zeros(dim);

tm = 0;
sm = 3;
xMu = 300;
yMu = 500;
xStd = 100;
yStd = 100;

for i  = 1:size(square,1)
    for j = 1:size(square,2)
%         x = cos(2*pi*j*tm/dim);
%         y = cos(2*pi*i*sm/dim);
%         x=0;
%         square(i,j) = x+y;
        z = cos(2*pi*(tm*j + sm*i)/dim);
        x = 1/(xStd*sqrt(2*pi))*exp(-0.5*((j-xMu)/xStd)^2);
        y = 1/(yStd*sqrt(2*pi))*exp(-0.5*((i-yMu)/yStd)^2);
        
        square(i,j) = z*x*y;
    end
end
figure;
imagesc(square)


time = 1000;
for t = 1:1000
    curve(t) = exp((t-500)/100)/(exp((t-500)/100)+1);
end
figure;plot(curve)





