
clearvars;
close all;
% go into folder that holds data 
%cd('Data');
result=zeros(1,3);

data = xlsread(['Results_FixationalSearch_xyz']); % use the xlsread command to read the subjects datasheets...
data_sorted = sortrows(data,3); %use the sortrows command
[r,c] = size(data_sorted); %r = rows, c = columns

%condition 1
j1 = 1; %counter for condition 1
while data_sorted(j1,3) == 1
    result(1,1) = result(1,1) + data_sorted(j1,2); %sum of times for condition 1
    j1=j1+1;
end
result(1,1) = result(1,1)/(j1 - 1); %divide sum by # of trials for condition 1

%condition 2
j2=j1; %counter for condition 2
while data_sorted(j2,3) == 2
    result(1,2) = result(1,2) + data_sorted(j2,2); %sum of times for condition 2
    j2=j2+1;
end
result(1,2) = result(1,2)/(j2 - j1); %divide sum by # of trials for condition 2

%condition 3
%condition 2
j3=j2; %counter for condition 2
while j3 < r+1
    result(1,3) = result(1,3) + data_sorted(j3,2); %sum of times for condition 3
    j3=j3+1;
end
result(1,3) = result(1,3)/(j3 - j2); %divide sum by # of trials for condition 3

c = categorical({'Condition 1','Condition 2','Condition 3'});
bar(c,result)

save Result.mat,result;
