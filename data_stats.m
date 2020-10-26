clearvars;

foldername = uipickfiles;
if ~iscell(foldername)
    if foldername == 0 
        fprintf('User cancelled folder selection. Silently exiting...\n');
        return;
    end
end

%initializes structures/fields
for i = 1:3
    OU_Data.(strcat('RT', num2str(i))) = [];
    OU_Data.(strcat('percentCorrect', num2str(i))) = [];
    OU_Data.(strcat('SID', num2str(i))) = strings(0);
    AE_Data.(strcat('RT', num2str(i))) = [];
    AE_Data.(strcat('percentCorrect', num2str(i))) = [];
    AE_Data.(strcat('SID', num2str(i))) = strings(0);
    FE_Data.(strcat('RT', num2str(i))) = [];
    FE_Data.(strcat('percentCorrect', num2str(i))) = [];
    FE_Data.(strcat('SID', num2str(i))) = strings(0);
    NDE_Data.(strcat('RT', num2str(i))) = [];
    NDE_Data.(strcat('percentCorrect', num2str(i))) = [];
    NDE_Data.(strcat('SID', num2str(i))) = strings(0);
    DE_Data.(strcat('RT', num2str(i))) = [];
    DE_Data.(strcat('percentCorrect', num2str(i))) = [];
    DE_Data.(strcat('SID', num2str(i))) = strings(0);
end

for i = 1:6
    OU_Data.displacementCounts.(strcat('bin', num2str(i))) = [];
    AE_Data.displacementCounts.(strcat('bin', num2str(i))) = [];
    FE_Data.displacementCounts.(strcat('bin', num2str(i))) = [];
    NDE_Data.displacementCounts.(strcat('bin', num2str(i))) = [];
    DE_Data.displacementCounts.(strcat('bin', num2str(i))) = [];
end

OU_Data.IncorrectImage = struct;
AE_Data.IncorrectImage = struct;
FE_Data.IncorrectImage = struct;
NDE_Data.IncorrectImage = struct;
DE_Data.IncorrectImage = struct;

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
            
            %finds eye type
            [fileParsed,matches] = split(filename,'_');
            suffix = fileParsed(3,1);
            [eyeType,matches] = split(suffix,'.');
            eye = eyeType(1,1);
            eye = convertCharsToStrings(eye(1,1));
            if eye == "OU"
                OU_Data = find_stats(filename, OU_Data);
            elseif eye == "AE"
                AE_Data = find_stats(filename, AE_Data);
            elseif eye == "FE"
                FE_Data = find_stats(filename, FE_Data);
            elseif eye == "NDE"
                NDE_Data = find_stats(filename, NDE_Data);
            elseif eye == "DE"
                DE_Data = find_stats(filename, DE_Data);
            end
        end
%     catch
%         disp('err');
%     end
end

save ('OU_Data.mat', 'OU_Data');
save ('AE_Data.mat', 'AE_Data');
save ('FE_Data.mat', 'FE_Data');
save ('NDE_Data.mat', 'NDE_Data');
save ('DE_Data.mat', 'DE_Data');
