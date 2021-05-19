function [spike_psds, f_rate] = rateSpectralAnalysis(spike_trains, space_axis)
%computes the psd for each texture for each afferent for each neuron for
%each trial.

spike_psds = spike_trains;

for i = 1:size(spike_trains, 1) % for every texture
    for j = 1:size(spike_trains, 2) % for every aff class
        trials = spike_trains{i,j}; % grab all trials
        
         % DO THE FIRST ONE
        train = trials(1, :); %grab train
        train = train - mean(train); %demean
        train = train';
        samp_period = space_axis(1);
        fs = 1/samp_period; % samples per mm

        % Welch's implementation
        n_windows = 2;
        win_length = round(length(train) / n_windows);
        [temp_fft, f] = pwelch(train, win_length,[],[], fs);

        %normalize 
        temp_fft = temp_fft./(max(temp_fft));
        psd_trials= temp_fft';
        
        % now do the rest, vert catting.
        
        for k = 2:size(trials, 1) % for each trial run
            train = trials(k, :); %grab train
            train = train - mean(train); %demean
            train = train';
            samp_period = space_axis(1);
            fs = 1/samp_period; % samples per mm
            
            % Welch's implementation
            n_windows = 2;
            win_length = round(length(train) / n_windows);
            [temp_fft, f] = pwelch(train, win_length,[],[], fs);

            %normalize 
            temp_fft = temp_fft./(max(temp_fft));
            psd_trials = vertcat(psd_trials, temp_fft');
        end
        spike_psds{i,j} = psd_trials;
    end
end
f_rate = f;
end

