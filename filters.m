function [Lx_filt, Ly_filt, velocityL, accelerationL, fig] = filters(Lx,Ly,f,Time)
%median
if f == 1
    Lx_filt = medfilt1(Lx);
    Ly_filt = medfilt1(Ly);
%gaussian
elseif f == 2
    h = fspecial('gaussian', [11 1], 11/16);
    Lx_filt = conv2(Lx, h, 'same');
    Ly_filt = conv2(Ly, h, 'same');
%smoothing
else
    Lx_filt=smooth(Lx);
    Ly_filt=smooth(Ly);
end

%sgolay filter and plots data
order = 2;
frame = [21];
for q = 1:length(frame)
    framelen = frame(q);
    Lx_filt = sgolayfilt(Lx_filt,order,framelen);
    Ly_filt = sgolayfilt(Ly_filt,order,framelen);
    
    fig = framelen*10+f
    figure(fig);
    h1 = subplot(3,1,1);
    plot(Time*1000,Lx_filt,'b',Time*1000,Ly_filt,'r')
    subplot(3,1,2)
    v = vecvel([Lx_filt Ly_filt], 500); % vecvel is to calculate velocity by using a running window
    velocityL = sqrt(v(:,1).^2 + v(:,2).^2); % vertical and horizontal component
    plot(Time*1000,velocityL)
    hold on;
    plot(5.*ones(1, round(Time(end-1)*1000, 0)));
    plot(7.*ones(1, round(Time(end-1)*1000, 0)));
    subplot(3,1,3)
    distanceL = diff(velocityL);
    accelerationL = distanceL./diff(Time);
    plot(Time(1:end-1)*1000,accelerationL)
end