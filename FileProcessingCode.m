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
            % foldername variable has saved all the files you want to load
            % in. Scary .txt file load in here
            % load list
%             fields_per_line_string = 4; 
%             fields_per_line_float = 38;
%             fmt = [repmat('%s',1,fields_per_line_string) repmat('%f',1,fields_per_line_float)];
% %           [filenameparse,matches] = split(foldername(1,folders),'/');
%             file = filenameparse(6,1);
%             fid = fopen([char(file)], 'rt');
%             if folders == 1
%                 filebycolumn1 = textscan(fid, fmt, 'Delimiter', '\t', 'headerlines', 1);
%             else
%                 filebycolumn2 = textscan(fid, fmt, 'Delimiter', '\t','headerlines', 1);
%             end
            
            %imports data
            [foldernameparse,matches] = split(foldername(1,folders),'/');
            filename = foldernameparse(6,1);
            
            [fileParsed,matches] = split(filename,'_');
            newFileName = strcat(fileParsed(1,1), '_', fileParsed(2,1), '_', fileParsed(3,1), '.txt');
            newFileName = char(newFileName(1,1));
            
            data = importdata(char(filename(1,1)));
            fmt = [repmat('%s\t',1,4) repmat('%f\t',1,38) ' \n'];
            [len, wid] = size(data.data);
            
            %uploads data
            if folders == 1
                fid = fopen(newFileName, 'w');
                fprintf(fid,'Name\teyeTested\tscene\ttarget\timageNum\tRT\tcondition\tcorrect\tfixDuration\tcenterX\tcenterY\tsignedDisplacement\timageRect1\timageRect2\timageRect3\timageRect4\ttargetRect1\ttargetRect2\ttargetRect3\ttargetRect4\tdistractor1Rect1\tdistractor1Rect2\tdistractor1Rect3\tdistractor1Rect4\tdistractor2Rect1\tdistractor2Rect2\tdistractor2Rect3\tdistractor2Rect4\tdistractor3Rect1\tdistractor3Rect2\tdistractor3Rect3\tdistractor3Rect4\tdistractor4Rect1\tdistractor4Rect2\tdistractor4Rect3\tdistractor4Rect4\txclick\tyclick\tdisplacements\tjitter\tsignedDisplacementDegrees\tsignedDisplacementPixels\n'); % this is the header
                for jj = 1 : len
                    fprintf(fid,fmt,data.textdata{jj+1,1}, data.textdata{jj+1,2},data.textdata{jj+1,3},data.textdata{jj+1,4},data.data(jj,1:38));
                end
            else
                fid = fopen(newFileName, 'ab');             
                for jj = 1 : len
                    fprintf(fid,fmt,data.textdata{jj+1,1}, data.textdata{jj+1,2},data.textdata{jj+1,3},data.textdata{jj+1,4},data.data(jj,1:38));
                end
            end
        end
%     catch
%         disp('err');
%     end
end

% save (newFileName, 'combinedFile');

% allData1 = [];
% for i = 5:length(filebycolumn1)
%     allData1(:,end+1) = filebycolumn1{1,i}
% end
% 
% allData2 = [];
% for i = 5:length(filebycolumn2)
%     allData2(:,end+1) = filebycolumn2{1,i}
% end
% 
% finalVal = cat(1,allData1,allData2);
%     
