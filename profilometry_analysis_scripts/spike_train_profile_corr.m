function [correlations] = spike_train_profile_corr(gel_psd, f_gel, no_gel_psd, f_no_gel, aff_psds, f_rate)
%spike_train_profile_corr calculates the corrcoef after interpolating as
%necessary between the power spectra of the three afferent spike trains and
%the gel profile and the profile without a gel.

%first: interpolate so that everything is on the same frequency axis.

%do each of the six correlations
correlations = zeros(1,6);
for i = 1:3
    aff_psd = aff_psds(:,i);
    largest_rate = f_rate(end);
    
    % crop off frequency axis and gel power level vector where freq axis is
    % greater than largest rate. TO DO
    
    aff_psd_interp = interp1(f_rate, aff_psd, f_gel);
    
    R_gel = corrcoef(aff_psd_interp, gel_psd);
    
    % crop off frequency axis and no gel power level vector where freq axis is
    % greater than largest rate. TO DO
    
    aff_psd_interp_no_gel = interp1(f_rate, aff_psd, f_no_gel);
    
    R_no_gel = corrcoef(aff_psd_interp_no_gel, no_gel_psd);
    
    correlations(i) = R_gel(1,2);
    correlations(i+3) = R_no_gel(1,2);
    
    
end
end

