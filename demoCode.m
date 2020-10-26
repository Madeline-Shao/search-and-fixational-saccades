clearvars;

filename = 'FixationSearch_DEN_OU_1.txt';
data = importdata(filename)
[len wid] = size(data.data)
fmt = [repmat('%s\t',1,4) repmat('%f\t',1,38) ' \n']
fid = fopen(['test','.txt'], 'w')
fprintf(fid,'Name\teyeTested\tscene\ttarget\timageNum\tRT\tcondition\tcorrect\tfixDuration\tcenterX\tcenterY\tsignedDisplacement\timageRect1\timageRect2\timageRect3\timageRect4\ttargetRect1\ttargetRect2\ttargetRect3\ttargetRect4\tdistractor1Rect1\tdistractor1Rect2\tdistractor1Rect3\tdistractor1Rect4\tdistractor2Rect1\tdistractor2Rect2\tdistractor2Rect3\tdistractor2Rect4\tdistractor3Rect1\tdistractor3Rect2\tdistractor3Rect3\tdistractor3Rect4\tdistractor4Rect1\tdistractor4Rect2\tdistractor4Rect3\tdistractor4Rect4\txclick\tyclick\tdisplacements\tjitter\tsignedDisplacementDegrees\tsignedDisplacementPixels\n') % this is the header

for jj = 1 : len
    fprintf(fid,fmt,data.textdata{jj+1,1}, data.textdata{jj+1,2},data.textdata{jj+1,3},data.textdata{jj+1,4},data.data(jj,1:38));
end
% 
% for jj = 1 : len
%     fprintf(fid,fmt,data.textdata{jj+1,1}, data.textdata{jj+1,2},data.textdata{jj+1,3},data.textdata{jj+1,4},data.data(1:38));
% end