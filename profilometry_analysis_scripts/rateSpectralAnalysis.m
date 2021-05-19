function [spike_psds, f_rate] = rateSpectralAnalysis(aff_spike_trains, space_axis)
%computes the psd for each afferent PSTH

figure;

%raster
% subplot(
% demean
spike_psds = [];
for i = 1:3
    spike_train = aff_spike_trains(i, :) - mean(aff_spike_trains(i, :)); %demean
    spike_train = spike_train';
    samp_period = space_axis(2);
    fs = 1/samp_period; % samples per mm
    % Welch's implementation
    n_windows = 2;
    win_length = round(length(spike_train) / n_windows);
    [temp_fft, f] = pwelch(spike_train, win_length,[],[], fs);
    
    %normalize 
    temp_fft = temp_fft./(max(temp_fft));
    spike_psds = horzcat(spike_psds, temp_fft);
    subplot(2,2, i+1)
    plot(f, temp_fft)
    xlabel("Frequency (1/mm)")
    ylabel("Power")
    yticks([]);
    yticklabels({});
    xlim([0 4]);
end
f_rate = f;
end

