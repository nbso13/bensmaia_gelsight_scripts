function [mean_y_trace, fy] = freqAnalysisFullNew(data)
%freqAnalysisFull takes the 1d F transform down every column. Averages all
%columns of the result. Returns a 1 d result.
samp_freq_y = 1/(data.y_res); 
size_y = size(data.profile, 1);

%demeaning
data.profile(0 == data.profile) = nan;
data_demeaned = data.profile - nanmean(data.profile, 'all');
data_demeaned(isnan(data_demeaned)) = 0;

% dft 1d
NFFTY = 2^nextpow2(size_y);

% Find Y freq space
fy = samp_freq_y/2*linspace(0,1,NFFTY/2+1);
amp_mat = fft(data_demeaned, NFFTY);
P2 = abs(amp_mat/size_y);
P1 = P2(1:NFFTY/2+1, :);
P1(:,2:end-1) = 2*P1(:,2:end-1);
mean_y_trace = mean(P1, 2); %single sided




% 
% Ly = size_y;
% Lx = size_x;
% fx = (samp_freq_x/(Lx))*(0:(Lx-1));
% fy = (samp_freq_y/(Ly))*(0:(Ly-1));
end