pathTo = '/Users/madeline/Dropbox/Search and fixational saccades';

cd(pathTo)
addpath(pwd);


%read in list
fields_per_line = 25;
fmt = [ repmat('%s',1,2) repmat('%f',1,23)];
fid = fopen(['ImageList','.txt'], 'rt')
filebycolumn = textscan(fid, fmt, 'Delimiter', '\t');
randShuffle = Shuffle(1:length(filebycolumn{1}));
images = filebycolumn{1}(randShuffle);
targets = filebycolumn{2}(randShuffle);
condition = filebycolumn{3}(randShuffle);
imNum = filebycolumn{4}(randShuffle);
targetRect = [filebycolumn{1,6} filebycolumn{1,7} filebycolumn{1,8} filebycolumn{1,9}];
distractorRect = [filebycolumn{1,10},filebycolumn{1,11},filebycolumn{1,12},filebycolumn{1,13},...
    filebycolumn{1,14},filebycolumn{1,15},filebycolumn{1,16},filebycolumn{1,17},...
    filebycolumn{1,18},filebycolumn{1,19},filebycolumn{1,20},filebycolumn{1,21},...
    filebycolumn{1,22},filebycolumn{1,23},filebycolumn{1,24},filebycolumn{1,25}];
targetRectRand = targetRect(randShuffle,:);
distractorRectRand = distractorRect(randShuffle,:);


for i =1:length(targets)
    scene=char(images(i));
    trial=i;
    target_word = char(targets(i));
    trial=trial+1;
    imdata=imread([[pathTo, '/Images/'],scene(find(~isspace(scene)))]); % images are read in
    imshow(imdata);
end