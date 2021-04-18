function [correlations, spectra_figure] = spike_train_profile_corr(texture_name, gel_psd, f_gel, no_gel_psd, f_no_gel, aff_psds, f_rate)
%spike_train_profile_corr calculates the corrcoef after interpolating as
%necessary between the power spectra of the three afferent spike trains and
%the gel profile and the profile without a gel.

%first: interpolate so that everything is on the same frequency axis.

%do each of the six correlations
correlations = zeros(1,6);
aff_names = [" PCs", " RAs", " SAs"];
line_thick = 1.5;
transparency = 0.6;
spectra_figure = figure;
for i = 1:3
    aff_psd = aff_psds(:,i);
    largest_rate = f_rate(end);
    
    % crop off frequency axis and gel power level vector where freq axis is
    % greater than largest rate. TO DO
    
    aff_psd_interp = interp1(f_rate, aff_psd, f_gel);
    known_freq = ~isnan(aff_psd_interp);
    aff_psd_interp = aff_psd_interp(known_freq);
    f_gel = f_gel(known_freq);
    gel_psd = gel_psd(known_freq);
    
    R_gel = corrcoef(aff_psd_interp, gel_psd);
    
    % crop off frequency axis and no gel power level vector where freq axis is
    % greater than largest rate. TO DO
    
    aff_psd_interp_no_gel = interp1(f_rate, aff_psd, f_no_gel);
    
    known_freq = ~isnan(aff_psd_interp_no_gel);
    aff_psd_interp_no_gel = aff_psd_interp_no_gel(known_freq);
    f_no_gel = f_no_gel(known_freq);
    no_gel_psd = no_gel_psd(known_freq);
    
    R_no_gel = corrcoef(aff_psd_interp_no_gel, no_gel_psd);
    
    correlations(i) = R_gel(1,2);
    correlations(i+3) = R_no_gel(1,2);
    
    
    
    
    subplot(3,1, i);
    hold on
    p1 = plot(f_gel, gel_psd, 'cyan', 'LineWidth', line_thick);
    p1.Color(4) = transparency;
    p2 = plot(f_no_gel, no_gel_psd, 'red', 'LineWidth', line_thick);
    p2.Color(4) = transparency;
    p3 = plot(f_rate, aff_psd, 'black', 'LineWidth', line_thick);
    p3.Color(4) = transparency;
    
    title(strcat(texture_name, aff_names(i)));
    xlim([0 4]);
    xlabel("Frequency (1/mm)")
    ylabel("Power")
    yticks([]);
    yticklabels({});
    if i ==1
        legend(["Gel", "No gel", "Recorded Spike Train"])
    end
    
end
sgtitle("Comparing Power Spectra")

end

