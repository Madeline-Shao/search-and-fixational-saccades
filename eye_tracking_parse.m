clearvars;

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
            data = importdata(char(filename(1,1)));
            for trial = 1:max(data(:,4))
                indx = find(data(:,4)==trial)   
            end
            
            dataByTrial = [];
            pixelSizeMm = .38;
            displayDistanceMm = 570; 
            pix2deg = atand(pixelSizeMm/displayDistanceMm);
            resolutionWidth = 1024;
            resolutionHeight = 768;
            cx = resolutionWidth/2;
            cy = resolutionHeight/2;
            
            %parses eye tracking data
            for i = 3:length(data)       
                trialPrev = data(i-1,4);
                trial = data(i, 4);
                if trial == trialPrev
                    dataByTrial(end+1, 1) = data(i-1, 2);
                    dataByTrial(end, 2) = data(i-1,3);
                    dataByTrial(end, 3) = data(i-1,1);
                else
                    dataByTrial(end+1, 1) = data(i-1, 2);
                    dataByTrial(end, 2) = data(i-1,3);
                    dataByTrial(end, 3) = data(i-1,1);
                    Lx = dataByTrial(:,1);
                    Ly = dataByTrial(:, 2);
                    Time = dataByTrial(:,3);
                    
                    %FIX THIS
                    Rx = Lx;
                    Ry= Ly;
                    Lx = pix2deg*(Lx - cx);
                    Rx = pix2deg*(Rx - cx);
                    Ly = pix2deg*(-Ly + cy);
                    Ry = pix2deg*(-Ry + cy);
                    Time = (Time - Time(1))/1000; % from milliseconds to seconds
                    samplingRate = round(1/diff(Time(1:2)));
                    
                    %removing blink artifacts
                    Lv = [0; diff(sqrt(Lx.^2 + Ly.^2))./diff(Time)];
                    Rv = [0; diff(sqrt(Rx.^2 + Ry.^2))./diff(Time)];
                    positionThreshold = 20;  % 30 deg

                    velocityThreshold = 1000; % deg/sec % do not make it smaller than 1000 
                    toremoveL = abs(Lx)>positionThreshold | abs(Ly)>positionThreshold | Lv>velocityThreshold;

                    toremoveR = abs(Rx)>positionThreshold | abs(Ry)>positionThreshold | Rv>velocityThreshold;

                    dangerZone = 300; % ms, Before and After blink
                    toremoveL = conv(double(toremoveL),ones(round(dangerZone*samplingRate/1000),1),'same')>0;
                    toremoveR = conv(double(toremoveR),ones(round(dangerZone*samplingRate/1000),1),'same')>0;
                    Lx(toremoveL) = NaN;
                    Ly(toremoveL) = NaN;
                    Rx(toremoveL) = NaN;
                    Ry(toremoveL) = NaN;
                    Lx(toremoveR) = NaN;
                    Ly(toremoveR) = NaN;
                    Rx(toremoveR) = NaN;
                    Ry(toremoveR) = NaN;
                    figure;              
                    plot(Time,Lx,'.r',Time,Ly,'.b',Time,Rx,'.g',Time,Ry,'.m'); hold on; 
                    % until here for saccade code
                    window = 30; % ms
                    n = round(window*samplingRate/1000);
                    Lx = medfilt1(Lx,n);
                    Ly = medfilt1(Ly,n);
                    Rx = medfilt1(Rx,n);
                    Ry = medfilt1(Ry,n);
                    figure(2);
                    plot(Time,Lx,'-r',Time,Ly,'-b',Time,Rx,'-g',Time,Ry,'-m'); hold on;

                %     figure;

                %     plot(Lx,Ly,'.r',Rx,Ry,'.b');

                %     xlim([-10 10]); ylim([-10 10]);
                    %fixationStability = repmat(GetEmptyStruct,1,2);
                    %leftEye = GetEmptyStruct;
                    %rightEye = GetEmptyStruct;
                    GetDensity(Lx, Ly, 1);
                    %blah(dataByTrial);
                    dataByTrial = [];
                end
                if i == length(data)
                    dataByTrial(end+1, 1) = data(i, 2);
                    dataByTrial(end, 2) = data(i,3);
                    dataByTrial(end, 3) = data(i,1);
                    %blah(dataByTrial);
                end
            end
        end
%     catch
%         disp('err');
%     end
end
