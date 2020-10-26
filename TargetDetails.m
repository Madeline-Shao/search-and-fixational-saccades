clearvars;
load('TargetDetails.mat');

%initializes variables
heightPixels = zeros(length(Y1), 1);
widthPixels = zeros(length(X1), 1);
areaPixels = zeros(length(X1), 1);
areaDegrees = zeros(length(X1), 1);
heightDegrees = zeros(length(X1), 1);
widthDegrees = zeros(length(X1), 1);

%calculates values
for i = 1:length(X1)
    heightPixels(i,1) = Y2(i,1) - Y1(i,1);
    widthPixels(i,1) = X2(i,1) - X1(i,1);
    areaPixels(i,1) = heightPixels(i,1) * widthPixels(i,1);
    areaDegrees(i,1) = areaPixels(i,1) * 0.0382;
    heightDegrees(i,1) = heightPixels(i,1) * 0.0382;
    widthDegrees(i,1) = widthPixels(i,1) * 0.0382;
    if (heightDegrees(i,1) < 1 || widthDegrees(i,1) < 1) && i < 446 
        disp(i)
    end
end
figure(1)
histogram(areaDegrees);
figure(2)
histogram(heightDegrees);
figure(3)
histogram(widthDegrees);