function [temp] = calcPlot2dFFT(profile_struct)
%calcPlot2dFFT calculates and plots the 2dim fourier transform. - thanks to charles greenspon for code - 

temp = abs(fftshift(fft2(profile_struct.profile)));
fs = 1 ./ [profile_struct.x_res profile_struct.y_res];
nyq_f = fs ./ 2;
r_freq_ax = linspace(-nyq_f(2),nyq_f(2), size(profile_struct.profile,1));
c_freq_ax = linspace(-nyq_f(1),nyq_f(1), size(profile_struct.profile,2));

imagesc(r_freq_ax, c_freq_ax, temp); pbaspect([1 1 1])
xlabel('Frequency (Hz)'); ylabel('Frequency (Hz)');
xlim([-10 10]); ylim([-10 10])
end

