function [amp_mat, fx, fy] = freqAnalysisFull(data)
%freqAnalysisFull takes the full fourier transform

samp_freq_y = 1/(data.y_res); % once every x microns, -> 1/mm
samp_freq_x = 1/(data.x_res); 
[size_y, size_x] = size(data.profile);

%demeaning
data.profile(0 == data.profile) = nan;
data_demeaned = data.profile - nanmean(data.profile, 'all');
data_demeaned(isnan(data_demeaned)) = 0;

% dft 2d
NFFTY = 2^nextpow2(size_y);
NFFTX = 2^nextpow2(size_x);

% Find X and Y frequency spaces
fx = samp_freq_x/2*linspace(0,1,NFFTX/2+1);
fy = samp_freq_y/2*linspace(0,1,NFFTY/2+1);
amp_mat = fft2(data_demeaned, NFFTY,NFFTX);
amp_mat = abs(amp_mat(1:NFFTY/2+1, 1:NFFTX/2+1));


% 
% Ly = size_y;
% Lx = size_x;
% fx = (samp_freq_x/(Lx))*(0:(Lx-1));
% fy = (samp_freq_y/(Ly))*(0:(Ly-1));
end