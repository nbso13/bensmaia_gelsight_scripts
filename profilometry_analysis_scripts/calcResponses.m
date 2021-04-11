function [FRs_ts, FRs_gel, r, a, figure_handles] = calcResponses(ts, gel, ...
    aff_density, ppm, speed, len, samp_freq, top_neuron_number, ...
    texture_rates, neuron_selection_modes, plot_flag)
%calcResponses given parameters, scans textures in ts structs across hand
%and measures response in aff_pop. FRs are returned.
r = {};
s = {};
a = {};
figure_handles = {};
ts_structs = [ts, gel];

for i = 1:length(ts_structs)
    aff_pop = affpop_hand('D2d',aff_density);
    s{i} = stim_scan_shape(ts_structs(i).shape, ts_structs(i).offset, ppm, ...
        len, samp_freq, ts_structs(i).amp, speed, ts_structs(i).gel_flag);
    if plot_flag
        stim_fig = figure;
        plot(s{i})
        aff_fig = figure;
        plot(aff_pop)
        title("aff pop old")
    end
    
    resp = aff_pop.response(s{i});
    
    min_x = min(ts_structs(i).shape(:, 1)); %grabbing extent of pins on hand coordinate system
    max_x = max(ts_structs(i).shape(:, 1));
    min_y = min(ts_structs(i).shape(:, 2));
    max_y = max(ts_structs(i).shape(:, 2));
    
    loc = [min_x, max_x, min_y, max_y];
    
    [resp_new, aff_pop_new] = chooseNeurons(resp, neuron_selection_modes, ...
        texture_rates, top_neuron_number, aff_pop, loc);
    r{i} = resp_new;
    a{i} = aff_pop_new;
    
    if plot_flag
        aff_pop_new_fig = figure;
        plot(aff_pop_new)
        title("aff pop new")
        response_fig = figure;
        subplot(2,1,1);
        plot(resp_new)
        title(ts_structs(i).name);
    end
    %calculate mean frs and sd
    
    FRs{i,1} = resp_new.rate(aff_pop_new.iPC);
    FRs{i,2} = resp_new.rate(aff_pop_new.iRA);
    FRs{i,3} = resp_new.rate(aff_pop_new.iSA1);
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
        subplot(2,1,2)
        bar(x, means);
        hold on
        b = bar(x,means);
        b.FaceColor = 'flat';
        b.CData(1,:) = [1 0.62 0];
        b.CData(2,:) = [0.2 0.2 1];
        b.CData(3,:) = [0.2 0.85 0.2];
        er = errorbar(x,means,sem);
        er.Color = [0 0 0];
        er.LineStyle = 'none';
        title(strcat(ts_structs(i).name, "firing rate"));
        xticks(x)
        xticklabels({'PCs','RAs','SAs'})
        ylabel("Hz")
    end
    total_force = sum(s{i}.profile(1,:));
    total_area = ts_structs(i).area;
%     total_area ts_structs(i).x_axis(end) * ts_structs(i).y_axis(end)
%     total_pressure
    disp(strcat("Total force: ", num2str(total_force)));
    disp(strcat("Total area: ", num2str(total_area), " sq mm"));
    disp(strcat("Total Pressure: ", num2str(total_force/(total_area/1000000)), " Pa"));
    
    force_profile = shape2profilometry(ts_structs(i).shape, ...
        s{i}.profile(1,:), ts_structs(i).pins_per_mm);
    trace_profile = shape2profilometry(ts_structs(i).shape, ...
        s{i}.trace(1,:), ts_structs(i).pins_per_mm);
    
    if plot_flag
%         subplot(2,2,3); visualizeProfile(force_profile);
%         title(strcat(ts_structs(i).name, " force profile"));
%         subplot(2,2,4); visualizeProfile(trace_profile);
%         title(strcat(ts_structs(i).name, " trace profile"));
%         sgtitle(strcat(ts_structs(i).name, " Responses"));
        
        figure_handles{i, 1} = stim_fig;
        figure_handles{i, 2} = response_fig;
    end
end
FRs_ts = {};
FRs_gel = {};

for i = 1:length(FRs)
    FRs_ts{i} = FRs{1,i};
    FRs_gel{i} = FRs{2,i};
end

end

