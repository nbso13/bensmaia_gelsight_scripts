function [outStruct] = removeLowFreq(profStruct, stopBand, method)
%removeLowFreq takes out frequencies lower than stopBand in profile.


%charles method
%demean
profStruct.profile = profStruct.profile - mean(profStruct.profile(:));

fs = 1 ./ [profStruct.x_res, profStruct.y_res];
nyq_f = fs ./ 2;


if strcmp(method, 'charles')
    
    freq_ax_x = linspace(0,2*nyq_f(1), length(profStruct.x_axis));
    freq_ax_y = linspace(0,2*nyq_f(2), length(profStruct.y_axis));
    
    f = fft2(profStruct.profile);
    mag = abs(f);
    % how do I correct for frequency decay?
    phase = atan2(imag(f), real(f));
    % How do I frequency-wise divide the magnitude?
    
    [X, Y] = meshgrid(freq_ax_x, freq_ax_y);
    mask = (X < stopBand) & (Y < stopBand);
    mag(mask) = 0; %filter
    f_filtered = mag .* exp(1i*phase); %recombine phase
    filtered_image = ifft2(f_filtered);
    outStruct = profStruct;
    outStruct.profile = real(filtered_image);
    outStruct.profile = outStruct.profile - min(outStruct.profile(:));
    
elseif strcmp(method, "new_charles")
    fs = 1 ./ profStruct.x_res;
    nyq_f = fs ./ 2;
    % In the case that the profile isn't square get the appropriate segment length for each axis
    f_L_row = 2^(nextpow2(size(profStruct.profile,1)) - 1);
    f_ax_row = fs*(0:(f_L_row/2))/f_L_row;
    f_ax_row_2s = [f_ax_row,fliplr(f_ax_row(2:end-1))];
    
    f_L_col = 2^(nextpow2(size(profStruct.profile,2)) - 1);
    f_ax_col = fs*(0:(f_L_col/2))/f_L_col;
    f_ax_col_2s = [f_ax_col,fliplr(f_ax_col(2:end-1))];
    
    % Truncate to axis lengths and fft
    prof_trun = profStruct.profile(1:f_L_row, 1:f_L_col);
    prof_trun_fft = fft2(prof_trun);
    prof_trun_fft_mag = abs(prof_trun_fft);
    prof_trun_fft_mag_shifted = fftshift(prof_trun_fft_mag);
    x_vec = linspace(-nyq_f,nyq_f, size(prof_trun,1));
    y_vec = linspace(-nyq_f,nyq_f, size(prof_trun,2));
    figure; subplot(1,2,1)
    imagesc(x_vec, y_vec, prof_trun_fft_mag_shifted)
    xlim([-10, 10]); ylim([-10, 10])
    title("Truncated Profile fft")
    
    % Masking
    f_ax_row_2d = repmat(f_ax_row_2s',[1,size(prof_trun,2)]);
    f_ax_col_2d = repmat(f_ax_col_2s,[size(prof_trun,1),1]);
    f_t = stopBand; %stopBand? in example it's 0.05
    f_t_logit = f_ax_row_2d < f_t & f_ax_col_2d < f_t;
    prof_trun_fft_filt = prof_trun_fft;
    prof_trun_fft_filt(f_t_logit) = 0 + 1i;
    
    subplot(1,2,2);
    imagesc(linspace(-nyq_f,nyq_f, size(prof_trun,1)), linspace(-nyq_f,nyq_f, size(prof_trun,2)), abs(prof_trun_fft_filt))
    xlim([-10, 10]); ylim([-10, 10])
    title("Truncated Profile fft")
    
    
    prof_trun_filt = abs(ifft2(ifftshift(prof_trun_fft_filt)));
    outStruct = profStruct;
    outStruct.profile = prof_trun_filt; %-min(prof_trun_filt(:));
    
elseif strcmp(method, 'disk')
    fs = fs(1);
    nyq_f = fs / 2;
    spacing = 2 *nyq_f / (size(profStruct.profile,2)-1); % frequency difference between array spots
    disc_radius = round(stopBand/spacing); %how many indices we need the radius of the disc to be
    f = fftshift(fft2(profStruct.profile));
    [x, y] = size(f);
    D = disc_radius;
    mask = fspecial('disk', D) == 0;
    mask = imresize(padarray(mask, [floor((x/2)-D) floor((y/2)-D)], 1, 'both'), [x y]);
    masked_ft = f .* mask;
    filtered_image = ifft2(ifftshift(masked_ft), 'symmetric');
    outStruct = profStruct;
    outStruct.profile = filtered_image-min(filtered_image(:));
else
    error("method str unrecognized")
end



end

