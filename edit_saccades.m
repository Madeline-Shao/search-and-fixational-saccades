 function [dataLabels, delete] = edit_saccades(Lx_filt, Ly_filt, Time, dataLabels, orig_dataLabels)
%figure(413);
choice = 1;
while choice ~= 3
    choice = menu('Press add saccade, delete saccade, done, or undo','add saccade','delete saccade','done', 'undo');
    if choice == 1
        %adds selected saccades
        [pind,xs,ys] = selectdata('selectionmode','rect');
        %         plot(xs{end}, ys{end}, 'm*');
        %         plot(xs{end-1}, ys{end-1}, 'm*');
        if ~isempty(xs{end})
            start = xs{end}(1,1);
            [pind,xs,ys] = selectdata('selectionmode','rect');
            if ~isempty(xs{end}) 
                e = xs{end}(end, 1);
                t = 1;
                %switches start and end if needed
                if start > e
                    v = start;
                    start = e;
                    e = v;
                end
                %marks saccades in dataLabels
                for d = start:e
                    if mod(d,2) == 0
                        indices(t) = find(Time==d/1000);
                        t=t+1;
                        dataLabels(find(Time==d/1000)) = 1;
                    end
                end
                %plots saccades
                hold on
                for f = 1:length(indices)
                    plot(Time(indices(f))*1000, Lx_filt(indices(f)), 'm*');
                    plot(Time(indices(f))*1000, Ly_filt(indices(f)), 'm*');
                end
                %         for i = 1:length(xs{end})
                %             dataLabels(trial, find(Time==xs{end}(i,1)/1000)) = 1;
                %         end
                %         for i = 1:length(xs{end-1})
                %             dataLabels(trial, find(Time==xs{end-1}(i,1)/1000)) = 1;
                %         end
            else
                %if no points selected
                disp('Error. Please select a point');
            end
        else
            %if no points selected
            disp('Error. Please select a point');
        end
    elseif choice == 2
        %deletes selected saccades
        [pind,xs,ys] = selectdata('selectionmode','rect');
        for i = 1:length(xs{end})
            dataLabels(find(Time==xs{end}(i,1)/1000)) = 0;
        end
%         for i = 1:length(xs{end-1})
%             dataLabels(trial, find(Time==xs{end-1}(i,1)/1000)) = 0;
%         end
        delete = 1;
        return; 
    elseif choice == 4
        %undos selections by reverting to original data labels
        h1 = subplot(3,1,1);
        cla(h1);
        subplot(3,1,1);
        hold on
        plot(Time*1000, Lx_filt, 'b');
        plot(Time*1000, Ly_filt, 'r');
        for z = 1:length(orig_dataLabels)
            if orig_dataLabels(z) == 1
                plot(Time(z)*1000, Lx_filt(z), 'm*');
                plot(Time(z)*1000, Ly_filt(z), 'm*');
            end
        end
        [dataLabels, delete] = edit_saccades(Lx_filt, Ly_filt, Time, orig_dataLabels,orig_dataLabels);
        return;
    end
end
delete = 0;
