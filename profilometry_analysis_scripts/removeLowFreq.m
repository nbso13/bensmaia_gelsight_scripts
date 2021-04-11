function [outStruct] = removeLowFreq(profStruct, stopBand)
%removeLowFreq takes out frequencies lower than stopBand in profile.


%charles method
%demean
profStruct.profile = profStruct.profile - mean(profStruct.profile(:));

fs = 1 ./ [profStruct.x_res, profStruct.y_res];
nyq_f = fs ./ 2;
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







% disc method
% fs = 1 / profStruct.x_res;
% nyq_f = fs / 2;
% freq_ax = linspace(-nyq_f,nyq_f, size(profStruct.profile,2));
% spacing = 2 *nyq_f / (size(profStruct.profile,2)-1); % frequency difference between array spots
% disc_radius = round(stopBand/spacing); %how many indices we need the radius of the disc to be
% 
% f = fftshift(fft2(profStruct.profile));
% [x, y] = size(f);
% 
% D = disc_radius;
% mask = fspecial('disk', D) == 0;
% mask = imresize(padarray(mask, [floor((x/2)-D) floor((y/2)-D)], 1, 'both'), [x y]);
% masked_ft = f .* mask;
% filtered_image = ifft2(ifftshift(masked_ft), 'symmetric');
% outStruct = profStruct;
% outStruct.profile = filtered_image-min(min(filtered_image));
end

