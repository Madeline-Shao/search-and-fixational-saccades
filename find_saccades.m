function saccades = find_saccades(dataX, i)
plot(1:1000,dataX(i,:));
saccades = [];
[pind,xs,ys] = selectdata('selectionmode','rect')
saccades = [xs, ys];
choice = 1
while choice == 1
    choice = menu('Press next sacc, complete, or undo','next sacc','complete','undo last selection');
    %adds saccade
    if choice == 1
        [pind,xs,ys] = selectdata('selectionmode','rect')
        if iscell(xs)
            xs = xs{length(xs)};
            ys = ys{length(ys)};
        else
        end
        s = [xs,ys];
        saccades = [saccades; s];
    end
end
%undos selections
if choice == 3
    saccades = [];
    clear xy sy pind
    saccades = find_saccades(dataX, i);
end

% plot data selected as a sanity check
figure(2)
plot(1:1000,dataX(i,:));
hold on
plot(saccades(:,1),saccades(:,2),'m*')
pause()
close all
end