%%
% this script shows how to use the microsaccade detection method
% published in Otero-Millan et al. Journal of Vision 2014

% Set up variables --------------------------------------------------------
clear all;
close all;

p=1;
results = [];

foldername = uipickfiles;
if ~iscell(foldername)
    if foldername == 0
        fprintf('User cancelled folder selection. Silently exiting...\n');
        return;
    end
end

for folders=1:length(foldername)
    currentFolder = foldername{folders};
    if isdir(currentFolder)
        listing = dir(currentFolder);
        indices = find(~[listing.isdir]); % get indices of all files (but not folders)
    else
        indices = 1;
        [currentFolder,listing.name,ext] = fileparts(currentFolder);
        listing.name = [listing.name ext];
    end
    
    %     try
    % loop over all files in the selected folder
    for i=1:length(indices)
        %finds file name
        [foldernameparse,matches] = split(foldername(1,folders),'/');
        filename = foldernameparse(length(foldernameparse),1);
        [folder,matches] = split(foldername,filename{1});
        folder=folder{1};
        data = importdata(foldername{1});
%         fid = fopen('Lx_filt_AMA.txt', 'w');
%         fclose(fid);
%         fid = fopen('Ly_filt_AMA.txt', 'w');
%         fclose(fid);
%         fid = fopen('dataLabels_AMA.txt', 'w');
%         fclose(fid);
        for trial = 1:max(data(:,4))
            indx = find(data(:,4)==trial);
            currTrial = data(min(indx):max(indx),:);
            
            pixelSizeMm = .38;
            displayDistanceMm = 570;
            pix2deg = atand(pixelSizeMm/displayDistanceMm);
            resolutionWidth = 1024;
            resolutionHeight = 768;
            cx = resolutionWidth/2;
            cy = resolutionHeight/2;
            
            
            %FIX THIS
            Lx = currTrial(:,2);
            Ly = currTrial(:,3);
            Rx = Lx;
            Ry= Ly;
            Lx = pix2deg*(Lx - cx);
            Rx = pix2deg*(Rx - cx);
            Ly = pix2deg*(-Ly + cy);
            Ry = pix2deg*(-Ry + cy);
            Time = (currTrial(:,1) - currTrial(1,1))/1000; % from milliseconds to seconds
            samplerate = round(1/diff(Time(1:2)));
            
            %removing blink artifacts
            Lv = [0; diff(sqrt(Lx.^2 + Ly.^2))./diff(Time)];
            Rv = [0; diff(sqrt(Rx.^2 + Ry.^2))./diff(Time)];
            positionThreshold = 20;  % 30 deg
            
            velocityThreshold = 1000; % deg/sec % do not make it smaller than 1000
            toremoveL = abs(Lx)>positionThreshold | abs(Ly)>positionThreshold | Lv>velocityThreshold;
            toremoveR = abs(Rx)>positionThreshold | abs(Ry)>positionThreshold | Rv>velocityThreshold;
            
            dangerZone = 300; % ms, Before and After blink
            % this is vector with 0 no blink, 1 blink
            toremoveL = conv(double(toremoveL),ones(round(dangerZone*samplerate/1000),1),'same')>0;
            toremoveR = conv(double(toremoveR),ones(round(dangerZone*samplerate/1000),1),'same')>0;
            % variable here called 'blinks' where a blink is 1 and no blink
            % is 0 for each time point. 
            
            % =========================================================================
            % EDIT THIS ===============================================================
            % Make sure to fill the variables samples with the eye movement data in the
            % proper colums. In blinks mark all the samples that are not good data. In
            % most VOG systems it is best to remove 100 ms before and after the blink
            % to remove all artifacts.
            % =========================================================================
            session = 'test';
            % TODO: fix this so can accomodate binocular and monocular
            samples = [];
            samples(:,1) = Time;
            samples(:,2) = Lx;
            samples(:,3) = Ly;
            samples(:,4) = Rx;
            samples(:,5) = Ry;
            
            blinks = [];
            for j = 1:size(toremoveL)
                if toremoveL(j) == 1 || toremoveR(j) == 1
                    blinks(j) = 1;
                else
                    blinks(j) = 0;
                end
            end
            
            % XY values
            figure(1);
            plot(Time*1000, Lx, 'b',Time*1000, Ly, 'r');
            xticks(0:200:(ceil(max(Time*1000)/200)*200))
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.02, .5, 0.96/2]);
            
            % test for different filters, Avi added
%             figure(50)
%             Lx=smooth(Lx);
%             Ly=smooth(Ly);
%             h1 = subplot(3,1,1);
%             plot(Time,Lx,'b',Time,Ly,'r')
%             subplot(3,1,2)
%             v = vecvel([Lx Ly], 500) % vecvel is to calculate velocity by using a running window
%             vel = sqrt(v(:,1).^2 + v(:,2).^2) % vertical and horizontal component
%             plot(Time*1000,vel)
%             subplot(3,1,3)
%             distanceL = diff(vel)
%             accelerationL = distanceL./diff(Time); 
%             plot(Time(1:end-1)*1000,accelerationL)
             
             %removes beginning and end of data
             Time = Time(12:end-12);
             Lx = Lx(12:end-12);
             Ly = Ly(12:end-12);
             
            %1 = median, 2 = gaussian, 3 = smoothing
            for f = 1:3
                %filters data
                [Lx_filt, Ly_filt, velocityL, accelerationL, fig] = filters(Lx, Ly, f, Time);
%                 order = 2;
%                 frame = 21:20:61;
%                 for q = 1:length(frame)
%                     framelen = frame(q);
%                     Lx_filt = sgolayfilt(Lx_filt,order,framelen);
%                     Ly_filt = sgolayfilt(Ly_filt,order,framelen);
% 
%                     
%                     figure(framelen);
%                     h1 = subplot(3,1,1);
%                     plot(Time,Lx_filt,'b',Time,Ly_filt,'r')
%                     subplot(3,1,2)
%                     v = vecvel([Lx_filt Ly_filt], 500); % vecvel is to calculate velocity by using a running window
%                     vel = sqrt(v(:,1).^2 + v(:,2).^2); % vertical and horizontal component
%                     plot(Time*1000,vel)
%                     subplot(3,1,3)
%                     distanceL = diff(vel);
%                     accelerationL = distanceL./diff(Time);
%                     plot(Time(1:end-1)*1000,accelerationL)
                    
                    
                %end
            %end
            
            %Filter from Hafed
%             
%             b = zeros(15,1);
%             b(1) =  -4.3353241e-04 * 2*pi;
%             b(2) =  -4.8506188e-04 * 2*pi;
%             b(3) =  -2.0984645e-05 * 2*pi;
%             b(4) =   1.3669190e-03 * 2*pi;
%             b(5) =   3.0795928e-03 * 2*pi;
%             b(6) =   3.8369002e-03 * 2*pi;
%             b(7) =   2.6923104e-03 * 2*pi;
%             b(8) =   0.0;
%             b(9:15) = -b(7:-1:1);
             
%             order = 2;
%             frame = 21:20:101;
%             frame = [frame, 167];
%             for f = 1:length(frame)
%                 framelen = frame(f);
%                 Lx_filt = sgolayfilt(Lx,order,framelen);
%                 figure(framelen);
% %                 plot(Lx,':')
% %                 hold on
% %                 plot(Lx_filt,'.-')
% %                 legend('signal','sgolay')
%                 
%                 Ly_filt = sgolayfilt(Ly,order,framelen);
%                 figure(framelen);
%                 plot(Ly,':')
%                 hold on
%                 plot(Ly_filt,'.-')
%                 legend('signal','sgolay')
                

%                 h1 = subplot(3,1,1);
%                 plot(Time*1000, Lx_filt, 'b',Time*1000, Ly_filt, 'r');
% %                 hold on;
% %                 plot(Time*1000,Lx)
% %                 plot(Time*1000,Ly)
%                 
%                 hold off
%                 subplot(3,1,2)
%                 distanceL = hypot(diff(Lx_filt), diff(Ly_filt));
%                 velocityL = distanceL./diff(Time);
%                 plot(Time(1:end-1)*1000,velocityL)
%                 hold on
%                 plot(7.*ones(1, Time(end-1)*1000));
%                 
%                 title('velocity')
%                 set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.02, .5, 0.96/2]);
%                 subplot(3,1,3)
%                 hold off
%                 accelerationL = diff(velocityL)./diff(Time(1:end-1));
%                 plot(Time(1:end-2)*1000,accelerationL);
%                 hold on
%                 plot(threshold.*ones(1, Time(end-2)*1000));
%                 plot(-threshold.*ones(1, Time(end-2)*1000));
%             end
%             pause;
            
%             b = fspecial('gaussian',[11 1],11/6);
%             Lx_filt = conv2(Lx,b,'same');
%             Ly_filt = conv2(Ly,b,'same');            

            dataLabels = zeros(1, length(Time));
            
            threshold = [350]; %acceleration threshold
            for n = 1:length(threshold)    
%                 figure(3)
%                 h1 = subplot(3,1,1);
%                 plot(Time*1000, Lx_filt, 'b',Time*1000, Ly_filt, 'r');
%                 hold off
%                 subplot(3,1,2)
%                 distanceL = hypot(diff(Lx_filt), diff(Ly_filt));
%                 velocityL = distanceL./diff(Time);
%                 plot(Time(1:end-1)*1000,velocityL)
%                 hold on
%                 plot(7.*ones(1, Time(end-1)*1000));
                
                %finds saccades using velocity threshold
                c = 1;
                saccades = [];
                subplot(3,1,2);
                hold on;
                for m = 1:length(velocityL)
                    if velocityL(m) > 7
                        plot(Time(m)*1000, velocityL(m), 'm*');
                        saccades(c, 1) = m;
                        saccades(c, 2) = velocityL(m);
                        c = c+1;
                    end
                end
               
                
                %xticks(0:200:(ceil(max(Time*1000)/200)*200))
%                 title('velocity')
%                 set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.02, .5, 0.96/2]);
%                 subplot(3,1,3)
%                 hold off
%                 accelerationL = diff(velocityL)./diff(Time(1:end));
%                 plot(Time(1:end-2)*1000,accelerationL);
%                 hold on
%                 plot(threshold(n).*ones(1, Time(end-2)*1000));
%                 plot(-threshold(n).*ones(1, Time(end-2)*1000));
%                 for m = 1:length(accelerationL)
%                     if accelerationL(m) > threshold(n) || accelerationL(m) < -threshold(n)
%                         plot(Time(m)*1000, accelerationL(m), 'm*');
%                     end
%                 end
%                 
                if ~isempty(saccades)
                    [len, wid] = size(saccades);
                    
                    %plots beginning and end point of each saccade on the
                    %velocity graph
                    subplot(3,1,2);
                    hold on;
                    plot(Time(saccades(1,1))*1000, velocityL(saccades(1,1)), 'g*');
                    plot(Time(saccades(end,1))*1000, velocityL(saccades(end,1)), 'k *');
                    
                    
                    %adjusts the start point of the first saccade according to acceleration
                    %threshold and plots on position graph
                    s = saccades(1,1);
                    w=0;
                    while s > 1 && (accelerationL(s) > threshold(n) || accelerationL(s) < -threshold(n))
                        w=1;
                        s = s-1;
                        subplot(3,1,1);
                        hold on;
                        plot(Time(s)*1000, Lx_filt(s), 'm*');
                        plot(Time(s)*1000, Ly_filt(s), 'm*');
                        dataLabels(s) = 1;
                    end
                    subplot(3,1,1);
                    hold on
                    if w == 1
                        plot(Time(s+1)*1000, Lx_filt(s+1), 'g*');
                        plot(Time(s+1)*1000, Ly_filt(s+1), 'g*');
                        dataLabels(s+1) = 1;
                    else
                        plot(Time(s)*1000, Lx_filt(s), 'g*');
                        plot(Time(s)*1000, Ly_filt(s), 'g*');
                        dataLabels(s) = 1;
                    end
                    
                    %adjusts the end point of the last saccade according to acceleration
                    %threshold and plots on position graph
                    s = saccades(end,1);
                    w =0;
                    while s < len && (accelerationL(s) > threshold(n) || accelerationL(s) < -threshold(n))
                        w=1;
                        s = s+1;
                        subplot(3,1,1);
                        plot(Time(s)*1000, Lx_filt(s), 'm*');
                        plot(Time(s)*1000, Ly_filt(s), 'm*');
                        dataLabels(s) = 1;
                    end
                    subplot(3,1,1);
                    hold on
                    if w == 1
                        plot(Time(s-1)*1000, Lx_filt(s-1), 'k*');
                        plot(Time(s-1)*1000, Ly_filt(s-1), 'k*');
                        dataLabels(s-1) = 1;
                    else
                        plot(Time(s)*1000, Lx_filt(s), 'k*');
                        plot(Time(s)*1000, Ly_filt(s), 'k*');
                        dataLabels(s) = 1;
                    end
                    
                    %adjusts rest of the start and end points according to
                    %acceleration
                    for r = 1:len
                        if r ~= 1 && r~= len
                            plot(Time(saccades(r,1))*1000, Lx_filt(saccades(r,1)), 'm*');
                            plot(Time(saccades(r,1))*1000, Ly_filt(saccades(r,1)), 'm*');
                            dataLabels(saccades(r,1)) = 1;
                        end
                        %end point adjustments
                        if r < len && (saccades(r, 1) ~= saccades(r+1) - 1)
                            subplot(3,1,2);
                            plot(Time(saccades(r,1))*1000, velocityL(saccades(r,1)), 'k*');
                            s = saccades(r,1);
                            w=0;
                            while accelerationL(s) > threshold(n) || accelerationL(s) < -threshold(n)
                                w=1;
                                s = s+1;
                                subplot(3,1,1);
                                plot(Time(s)*1000, Lx_filt(s), 'm*');
                                plot(Time(s)*1000, Ly_filt(s), 'm*');
                                dataLabels(s) = 1;
                            end
                            subplot(3,1,1);
                            hold on
                            if w == 1
                                plot(Time(s-1)*1000, Lx_filt(s-1), 'k*');
                                plot(Time(s-1)*1000, Ly_filt(s-1), 'k*');
                                dataLabels(s-1) = 1;
                            else
                                plot(Time(s)*1000, Lx_filt(s), 'k*');
                                plot(Time(s)*1000, Ly_filt(s), 'k*');
                                dataLabels(s) = 1;
                            end
                        end
                        %start point adjustments
                        if r > 1 && (saccades(r, 1) ~= saccades(r-1) + 1)
                            subplot(3,1,2);
                            plot(Time(saccades(r,1))*1000, velocityL(saccades(r,1)), 'g*');
                            s = saccades(r,1);
                            w=0;
                            while s > 1 && accelerationL(s) > threshold(n) || accelerationL(s) < -threshold(n)
                                w=1;
                                s = s-1;
                                subplot(3,1,1);
                                plot(Time(s)*1000, Lx_filt(s), 'm*');
                                plot(Time(s)*1000, Ly_filt(s), 'm*');
                                dataLabels(s) = 1;
                            end
                            subplot(3,1,1);
                            hold on
                            if w == 1
                                plot(Time(s+1)*1000, Lx_filt(s+1), 'g*');
                                plot(Time(s+1)*1000, Ly_filt(s+1), 'g*');
                                dataLabels(s+1) = 1;
                            else
                                plot(Time(s)*1000, Lx_filt(s), 'g*');
                                plot(Time(s)*1000, Ly_filt(s), 'g*');
                                dataLabels(s) = 1;
                            end
                        end
                    end
                end
                
                %gui to edit saccades
                fig = figure(fig);
                h1 = subplot(3,1,1);
                delete = 1;
                dataLabels_1 = dataLabels;
                while delete ==1
                    [dataLabels_1, delete] = edit_saccades(Lx_filt, Ly_filt, Time, dataLabels_1, dataLabels);
                    cla(h1);
                    subplot(3,1,1);
                    hold on
                    plot(Time*1000, Lx_filt, 'b');
                    plot(Time*1000, Ly_filt, 'r');
                    for z = 1:length(dataLabels_1)
                        if dataLabels_1(z) == 1
                            plot(Time(z)*1000, Lx_filt(z), 'm*');
                            plot(Time(z)*1000, Ly_filt(z), 'm*');
                        end
                    end
                end
                dataLabels = dataLabels_1;
                
                %pause;
                %savefig(strcat(num2str(trial), 'figure_', num2str(threshold(n))));
                
            end
            
            %saves png of graphs
            if f == 1
                type = 'median';
            elseif f == 2
                type = 'gaussian';
            else
                type = 'smoothing';
            end
            %saveas(fig, strcat('trial ', num2str(trial), '_', type, '.png'));
            end
            close all;
%           
% 
%             hold off
%             figure(2);
%             plot(Time(2:end), velocityL);
%             pause;
%             
% 
%            
%             figure(1);
%             plot(Time, Rx, 'b');
%             hold on
%             plot(Time, Ry, 'r');
%             hold off
%             figure(2);
%             plot(Time(2:end), velocityR);
%             pause;
%             
%             G=fspecial('gaussian',[15,1],15/6)% what should sigma be?
%             ker = [1 0 -1];
%             f = conv(G,ker,'valid')
%             x_filt = filter(f,1,Lx); % lin to filter code
%             y_filt = filter(f,1,Ly); % lin to filter code
%             
%             for k = 2:length(x_filt)
%                 distanceL_filt(k-1) = sqrt((x_filt(k) - x_filt(k-1))^2 + (y_filt(k) - y_filt(k-1))^2);
%                 velocityL_filt(k-1)= distanceL_filt(k-1)/.002;
%             end
%             
%             figure(1);
%             plot(Time, x_filt, 'b');
%             hold on
%             plot(Time, y_filt, 'r');
%             hold off
%             figure(2);
%             plot(Time(2:end), velocityL_filt);
%             pause;
%             
%             G=fspecial('gaussian',[15,1],15/6)% what should sigma be?
%             ker = [1 0 -1];
%             f = conv(G,ker,'valid')
%             x_filt = filter(f,1,Rx); % lin to filter code
%             y_filt = filter(f,1,Ry); % lin to filter code
%             
%             for k = 2:length(x_filt)
%                 distanceR_filt(k-1) = sqrt((x_filt(k) - x_filt(k-1))^2 + (y_filt(k) - y_filt(k-1))^2);
%                 velocityR_filt(k-1)= distanceR_filt(k-1)/.002;
%             end
%             
%             figure(1);
%             plot(Time, x_filt, 'b');
%             hold on
%             plot(Time, y_filt, 'r');
%             hold off
%             figure(2);
%             plot(Time(2:end), velocityR_filt);
%             pause;
            
            %saves data
            [fileParsed,matches] = split(filename,'_');
            sid = fileParsed{3,1};
            
            fid = fopen(strcat('Lx_filt_', sid, '.txt'), 'a');
            Lx_pix = currTrial(12:end-12,2);
            fs = repmat('%f\t', 1, length(Lx_pix));
            fs = [fs, '\n'];
            fprintf(fid, fs, Lx_pix);
            fclose(fid);
            
            fid = fopen(strcat('Ly_filt_', sid, '.txt'), 'a');
            Ly_pix = currTrial(12:end-12,3);
            fs = repmat('%f\t', 1, length(Ly_pix));
            fs = [fs, '\n'];
            fprintf(fid, fs, Ly_pix);
            fclose(fid);
            
            fid = fopen(strcat('dataLabels_', sid, '.txt'), 'a');
            fs = repmat('%d\t', 1, length(dataLabels));
            fs = [fs, '\n'];
            dataL=transpose(dataLabels);
            fprintf(fid, fs, dataL);
            fclose(fid);
        end
    end
end
% %         
% %         %     %  samples(:,1)     timestamps of the recording in miliseconds
% %         %     %  samples(:,2)     horizontal position of the left eye in degrees
% %         %     %  samples(:,3)     vertical position of the left eye in degrees
% %         %     %  samples(:,4)     horizontal position of the right eye in degrees
% %         %     %  samples(:,5)     vertical position of the right eye in degrees
% %         %
% %         %     blinks = [];
% %         %     %  blinks           binary vector indicating for each sample if it
% %         %                   belons to a blink (1) or not (0)
% %         % =========================================================================
% %         % END EDIT THIS ===========================================================
% %         % =========================================================================
% %         
% %         % Loads the recording and prepares it for processing
% %         recording = ClusterDetection.EyeMovRecording.Create(folder, session, samples, blinks, samplerate);
% %         
% %         % Runs the saccade detection
% %         [saccades stats] = recording.FindSaccades();% breaking here
% %         
% %         if p == 1
% %             % Plots a main sequence
% %             enum = ClusterDetection.SaccadeDetector.GetEnum;
% %             figure
% %             subplot(2,2,1)
% %             plot(saccades(:,enum.amplitude),saccades(:,enum.peakVelocity),'o')
% %             set(gca,'xlim',[0 1.2*max(saccades(:,enum.amplitude))],'ylim',[0 1.2*max(saccades(:,enum.peakVelocity))]); % set axes
% %             xlabel('Saccade amplitude (deg)');
% %             ylabel('Saccade peak velocity (deg/s)');
% %             
% %             trialResults = ones(length(saccades(:,enum.amplitude)), 1);
% %             trialResults = trialResults * trial;
% %             trialResults = [trialResults, saccades(:,enum.amplitude), saccades(:,enum.peakVelocity)];
% %             results = [results; trialResults];
% %             
% %             % Plots the traces with the labeled microsaccades
% %             subplot(2,2,[3:4])
% %             plot(samples(:,1), samples(:,2:end));
% %             hold
% %             yl = get(gca,'ylim');
% %             u1= zeros(size(samples(:,1)))+yl(1);
% %             u2= zeros(size(samples(:,1)))+yl(1);
% %             u1((saccades(:,enum.startIndex))) = yl(2);
% %             u2(saccades(:,enum.endIndex)) = yl(2);
% %             u = cumsum(u1)-cumsum(u2);
% %             plot(samples(:,1), u,'k')
% %             
% %             xlabel('Time (ms)');
% %             ylabel('Eye Position (deg)');
% %             
% %             legend({'Left Horiz', 'Left Vert', 'Right Horiz' , 'Right Vert', 'Microsaccades'})
% %         end
% %         %pause()
% %         close all
% %     end
% % end
% % end

