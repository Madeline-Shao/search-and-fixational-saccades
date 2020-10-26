clear all;
close all;

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
        [foldernameparse,matches] = split(foldername(1,folders),'/');
        filename = foldernameparse(end,1);
        
        if folders == 1
            dataLabels = importdata(char(filename(1,1)));
        else
            dataX = importdata(char(filename(1,1)));
        end
    end
    %     catch
    %         disp('err');
    %     end
end

%user adds saccades in data
[len, wid] = size(dataX);
results = zeros(len,wid);
for i = 1:len
    saccades = find_saccades(dataX, i);
    for j = 1:length(saccades)
        results(i,saccades(j)) = 1;
    end
end

save('Results.mat', 'results');