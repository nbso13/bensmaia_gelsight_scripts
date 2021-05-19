function process_noise_pos

try
    load('pma_chubNoise_desiredpositions.mat')
    load('pma_noise.mat')
catch
    try
        load('\\bsl-somsrv1\data\processed\per\pma\desired\pma_chubNoise_desiredpositions.mat')
        load('\\bsl-somsrv1\data\processed\per\pma\desired\pma_noise.mat')
    catch
        error('Unable to find noise position files locally or on server.')
    end
end

stim_noise = L2_str.stim_pos{1}';

for n=1:length(desiredpositions)
    desiredpositions{n} = desiredpositions{n}(:,10001:30000);
end

desirednoisepos = horzcat(stim_noise,desiredpositions);

save desirednoisepos.mat desirednoisepos
