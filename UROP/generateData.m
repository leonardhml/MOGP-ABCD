function [ X,Y, xpred, ypred, output ] = generateData()

%%%
% Sample training and test inputs
%%%
x1 = rand(100,1);
x1 = sort(x1);
x1 = [x1(1:30); x1(60:100)];
x2 = rand(100,1);                      
f1 = @(x) 0.5 * x.^3;
f2 = @(x) 0.25 * x.^3;
%f1 = @(x) sin(6*x1);
%f2 = @(x) sin(4*x2);
y1 = f1(x1) + normrnd(0, 0.01, length(x1),1);
y2 = f2(x2) + normrnd(0, 0.01, length(x2),1);

X.x1 = x1;
X.x2 = x2;
Y.y1 = y1;
Y.y2 = y2;
% 
% a = preprocessTrafficData('data/nus_traffic47/e103021891.csv', 1, 1, 12);
% X.x1 = a(1:2:end, 1);
% X.x2 = a(2:2:end, 1);
% Y.y1 = a(1:2:end, 2);
% Y.y2 = a(2:2:end, 3);

xpred = linspace(0, 1, 20)';                      
ypred = f1(xpred);
output = 1;

end

