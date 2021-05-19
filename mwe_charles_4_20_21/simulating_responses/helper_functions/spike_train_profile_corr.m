function [corrs, spectra_figure] = spike_train_profile_corr(texture_name, ...
    gel_psd, f_gel, no_gel_psd, f_no_gel, aff_psds, f_rate)
%spike_train_profile_corr calculates the corrcoef of each spike train PSD
%with the gel and then the no_gel PSD.

%% interpolate so gel and no_gel are onthe same frequency axis as the rates
gel_psd_interp = interp1(f_gel, gel_psd, f_rate);
norm_gel_psd_interp = gel_psd_interp./max(gel_psd_interp);
no_gel_psd_interp = interp1(f_no_gel, no_gel_psd, f_rate);
norm_no_gel_psd_interp = no_gel_psd_interp./max(no_gel_psd_interp);

%now everything uses f_rate.

%% plot the average PSD for each afferent class with the PSD of the gel and no gel
aff_names = [" PCs", " RAs", " SAs"];
line_thick = 1.5;
transparency = 0.6;
spectra_figure = figure;
for i = 1:3 % for each afferent class
    trial_psds = aff_psds{i}; %array with trial psd for each row
    mean_aff_psd = mean(trial_psds, 1);
    norm_mean_aff_psd = mean_aff_psd./max(mean_aff_psd);
    subplot(3,1, i);
    hold on
    p1 = plot(f_rate, norm_gel_psd_interp, 'cyan', 'LineWidth', line_thick);
    p1.Color(4) = transparency;
    p2 = plot(f_rate, norm_no_gel_psd_interp, 'red', 'LineWidth', line_thick);
    p2.Color(4) = transparency;
    p3 = plot(f_rate, norm_mean_aff_psd, 'black', 'LineWidth', line_thick);
    p3.Color(4) = transparency;
    
    if i == 1
        title(strcat(texture_name, aff_names(i)));
    else
        title(aff_names(i));
    end
    
    if i == 3
        xlabel("Frequency (1/mm)")
    end
    if i == 2
        ylabel("Normalized Power")
    end
    if i ==1
        strs = {'Gel', 'No gel', 'Recorded Spike Trains'}';
        colors = [[0 1 1]; [1 0 0]; [0 0 0]];
        leg = legend(color_legend(strs, colors));
        leg.Box = 0;
    end
    
    yticks([]);
    yticklabels({});
    
    ax = gca;
    ax.FontSize = 12;
    ax.FontWeight = 'bold';
end

corrs = cell(1,3);
%% compute the correlations for each trial PSD with the gel and no gel PSD
for i = 1:3 % for each afferent
    trial_psds = aff_psds{i};
    corrs_trials = zeros(size(trial_psds, 1), 2); %first row gel, second row no gel correlations
    for j = 1:size(trial_psds, 1)
    	R_gel = corrcoef(trial_psds(j, :), gel_psd_interp);
    	R_no_gel = corrcoef(trial_psds(j, :), no_gel_psd_interp);
        corrs_trials(j, 1) = R_gel(1,2);
        corrs_trials(j, 2) = R_no_gel(1,2);
    end
    corrs{i} = corrs_trials;
end
end

