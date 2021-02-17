function [FRs, figure_handles] = calcResponses(ts_structs, aff_pop, ppm, speed, len, loc, samp_freq, ramp_len)
%calcResponses given parameters, scans textures in ts structs across hand
%and measures response in aff_pop. FRs are returned.
r = {};
s = {};
FRs = {};
figure_handles = {};
for i = 1:length(ts_structs)
    s{i} = stim_scan_shape(ts_structs(i).shape, ts_structs(i).offset, ppm, ...
        len, samp_freq, ts_structs(i).amp, speed, ts_structs(i).gel_flag);
    stim_fig = figure;
    plot(s{i})
    resp = aff_pop.response(s{i});
%     fr{i} = calcFR(r{i}, plot_flag)
%     %take out neurons that fire less than 2 spikes per second
%     r_new = excludeNeurons(r, 2);
    r{i} = resp;
    response_fig = figure;
    plot(resp)
    title(ts_structs(i).name);
    %calculate mean frs and sd
    
    FRs{i,1} = resp.rate(aff_pop.iPC);
    FRs{i,2} = resp.rate(aff_pop.iRA);
    FRs{i,3} = resp.rate(aff_pop.iSA1);
    means = zeros(3,1);
    sem = means;
    for j = 1:3
        means(j)  = mean(FRs{i,j});
        sem(j) = std(FRs{i,j})/sqrt(length(FRs{i,j}));
    end
    FRs{i,4} = means;
    FRs{i, 5} = sem;
    x = [1,2,3];
    
    mean_response_fig = figure;
    bar(x, means);
    hold on
    b = bar(x,means);
    b.FaceColor = 'flat';
    b.CData(1,:) = [1 0.62 0];
    b.CData(2,:) = [0.2 0.2 1];
    b.CData(3,:) = [0.2 0.85 0.2];
    er = errorbar(x,means,means-sem,means+sem);
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
    title(strcat(ts_structs(i).name, "firing rate"));
    xticks(x)
    xticklabels({'PCs','RAs','SAs'})
    ylabel("Hz")
    
    force_profile = shape2profilometry(ts_structs(i).shape, ...
    s{i}.profile(1,:), ts_structs(i).pins_per_mm);
    force_profile_fig = figure; visualizeProfile(force_profile);
    title(strcat(ts_structs(i).name, " force profile"));
    disp(strcat("Total forces: ", num2str(sum(sum(force_profile.profile)))));
    
    trace_profile = shape2profilometry(ts_structs(i).shape, ...
    s{i}.trace(1,:), ts_structs(i).pins_per_mm);
    trace_profile_fig = figure; visualizeProfile(trace_profile);
    title(strcat(ts_structs(i).name, " trace profile"));
    
    figure_handles{i, 1} = stim_fig;
    figure_handles{i, 2} = response_fig;
    figure_handles{i, 3} = mean_response_fig;
    figure_handles{i, 4} = force_profile_fig;
    figure_handles{i, 5} = trace_profile_fig;
end

end

