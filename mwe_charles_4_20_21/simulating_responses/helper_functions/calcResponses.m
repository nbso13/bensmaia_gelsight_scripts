function [FRs_ts, FRs_gel, loc, r, figure_handles] = calcResponses(aff_pop_in, ts, gel, ...
     ppm, speed, len, samp_freq, top_neuron_number, ...
    texture_rates, neuron_selection_modes, plot_flag)
%calcResponses given parameters, scans textures in ts structs across hand
%and measures response in aff_pop. FRs are returned.
r = {};
s = {};
figure_handles = {};
ts_structs = [ts, gel];
response_fig = figure;
loc = zeros(2, 4);
for i = 1:length(ts_structs)
    aff_pop = aff_pop_in;
    s{i} = stim_scan_shape(ts_structs(i).shape, ts_structs(i).offset, ppm, ...
        len, samp_freq, speed, ts_structs(i).gel_flag);

    resp = aff_pop.response(s{i});
    
    min_x = min(ts_structs(i).shape(:, 1)); %grabbing extent of pins on hand coordinate system
    max_x = max(ts_structs(i).shape(:, 1));
    min_y = min(ts_structs(i).shape(:, 2));
    max_y = max(ts_structs(i).shape(:, 2));
    
    loc(i,:) = [min_x, max_x, min_y, max_y];
    
%     [resp_new] = chooseNeurons(resp, neuron_selection_modes, ...
%         texture_rates, top_neuron_number, aff_pop, loc);
%    
    resp_new = resp;
    r{i} = resp_new;
    
    if plot_flag
        figure(response_fig);
        subplot(2,3,i);
        plot(resp_new);
        title(ts_structs(i).name);
        
        ax = gca;
        ax.FontSize = 12;
        ax.FontWeight = 'bold';
        if i == 2
            ax.YAxis.Visible = 'off'; % remove y-axis
        end
    end
    %calculate mean frs and sd
    
    FRs{i,1} = resp_new.rate(resp_new.affpop.iPC);
    FRs{i,2} = resp_new.rate(resp_new.affpop.iRA);
    FRs{i,3} = resp_new.rate(resp_new.affpop.iSA1);
    means = zeros(3,1);
    sem = means;
    
    for j = 1:3
        means(j)  = mean(FRs{i,j});
        sem(j) = std(FRs{i,j})/sqrt(length(FRs{i,j}));
    end
    
    FRs{i,4} = means;
    FRs{i, 5} = sem;
    x = [1,2,3];
    
    if plot_flag
        subplot(2,3,i+3);
        bar(x, means);
        hold on
        b = bar(x,means);
        b.FaceColor = 'flat';
        b.CData(1,:) =  [255 127 0]/255;
        b.CData(2,:) =  [30 120 180]/255;
        b.CData(3,:) = [50 160 40]/255;
        er = errorbar(x,means,sem);
        er.Color = [0 0 0];
        er.LineStyle = 'none';
%         title(strcat(ts_structs(i).name, " firing rate"));
        xticks(x)
        xticklabels({'PCs','RAs','SAs'})
        ylabel("Hz")
        ax = gca;
        ax.FontSize = 12;
        ax.FontWeight = 'bold';
    end
    total_force = sum(s{i}.profile(1,:));
    total_area = ts_structs(i).area;
%     total_area ts_structs(i).x_axis(end) * ts_structs(i).y_axis(end)
%     total_pressure
    disp(strcat("Total force: ", num2str(total_force)));
    disp(strcat("Total area: ", num2str(total_area), " sq mm"));
    disp(strcat("Total Pressure: ", num2str(total_force/(total_area/1000000)), " Pa"));
  
end
for i = 1:length(s)
    figure;
    plot(s{i});
%     sgtitle(ts_structs(i).name);
end
FRs_ts = {};
FRs_gel = {};

for i = 1:length(FRs)
    FRs_ts{i} = FRs{1,i};
    FRs_gel{i} = FRs{2,i};
end

end

