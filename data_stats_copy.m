clearvars;

foldername = uipickfiles;
if ~iscell(foldername)
    if foldername == 0 
        fprintf('User cancelled folder selection. Silently exiting...\n');
        return;
    end
end

OU_Data.RT1 = zeros(1, length(foldername));
OU_Data.RT2 = zeros(1, length(foldername));
OU_Data.RT3 = zeros(1, length(foldername));
OU_Data.percentCorrect1 = zeros(1, length(foldername));
OU_Data.percentCorrect2 = zeros(1, length(foldername));
OU_Data.percentCorrect3 = zeros(1, length(foldername));
OU_Data.SID1 = strings(1, length(foldername));
OU_Data.SID2 = strings(1, length(foldername));
OU_Data.SID3 = strings(1, length(foldername));
displacementCounts.bin1 = zeros(1, length(foldername));
displacementCounts.bin2 = zeros(1, length(foldername));
displacementCounts.bin3 = zeros(1, length(foldername));
displacementCounts.bin4 = zeros(1, length(foldername));
displacementCounts.bin5 = zeros(1, length(foldername));
displacementCounts.bin6 = zeros(1, length(foldername));
subCount = 0;

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
    subCount = subCount + 1; %counter for # of subjects
%     try
        % loop over all files in the selected folder
        for i=1:length(indices)
            %finds file name and imports it
            [foldernameparse,matches] = split(foldername(1,folders),'/');
            filename = foldernameparse(length(foldernameparse),1);
            data = importdata(char(filename(1,1)));
            
            [fileParsed,matches] = split(filename,'_');
            newFileName = fileParsed(2,1);
            subName = char(newFileName(1,1)); %subject ID
            IncorrectImage.(strcat('', subName)) = {};
            
            
            data_sorted = sortrows(data.data,3); %use the sortrows command
            [r,c] = size(data_sorted); %r = rows, c = columns
            
            %condition 1
            j1 = 1; %counter for condition 1
            avg = 0;
            numCorrect = 0;
            numIncorrect = 0;
            while data_sorted(j1,3) == 1
                avg = avg + data_sorted(j1,2); %sum of times for condition 1
                if data_sorted(j1, 4) == 1
                    numCorrect = numCorrect + 1;
                else
                    %record incorrect images
                    numIncorrect = numIncorrect + 1;
                    imageNum = data_sorted(j1, 1);
                    IncorrectImage.(strcat('', subName))(1, numIncorrect) = data.textdata(imageNum + 1, 3);
                end
                j1=j1+1;
            end
            OU_Data.SID1(1, folders) = subName;
            OU_Data.RT1(1, folders) = avg/(j1 - 1); %divide sum by # of trials for condition 1
            OU_Data.percentCorrect1(1, folders) = numCorrect / (j1-1); %find percent correct
            
            %condition 2
            j2=j1; %counter for condition 2
            avg = 0;
            numCorrect = 0;
            while data_sorted(j2,3) == 2
                avg = avg + data_sorted(j2,2); %sum of times for condition 2
                if data_sorted(j2, 4) == 1
                    numCorrect = numCorrect +1;
                else
                    numIncorrect = numIncorrect + 1;
                    imageNum = data_sorted(j2, 1);
                    IncorrectImage.(strcat('', subName))(1, numIncorrect) = data.textdata(imageNum + 1, 3);
                end
                j2=j2+1;
            end
            OU_Data.SID2(1, folders) = subName;
            OU_Data.RT2(1, folders) = avg/(j2 - j1); %divide sum by # of trials for condition 2
            OU_Data.percentCorrect2(1, folders) = numCorrect / (j2 - j1);
            
            %condition 3
            j3=j2; %counter for condition 2
            avg = 0;
            numCorrect = 0;
            while j3 < r+1
                avg = avg + data_sorted(j3,2); %sum of times for condition 3
                if data_sorted(j3, 4) == 1
                    numCorrect = numCorrect +1;
                else
                    numIncorrect = numIncorrect + 1;
                    imageNum = data_sorted(j3, 1);
                    IncorrectImage.(strcat('', subName))(1, numIncorrect) = data.textdata(imageNum + 1, 3);
                end
                j3=j3+1;
            end
            OU_Data.SID3(1, folders) = subName;
            OU_Data.RT3(1, folders) = avg/(j3 - j2); %divide sum by # of trials for condition 3
            OU_Data.percentCorrect3(1, folders) = numCorrect / (j3 - j2);
            
            edges = [min(data_sorted(:,37)) -7 -3 0  3 7  max(data_sorted(:,37)) ];
            [N, edges] = histcounts(data_sorted(:,37),edges);
            for k = 1:6
                displacementCounts.(strcat('bin', num2str(k)))(1, folders) = N(1, k);
            end         
        end
%     catch
%         disp('err');
%     end
end

save ('OU_Data.mat', 'OU_Data');
save ('IncorrectImage.mat', 'IncorrectImage');
save ('displacementCounts.mat', 'displacementCounts');
