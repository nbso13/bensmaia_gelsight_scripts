%%% 2D FFT Testing
fp = '\\BENSMAIA-LAB\LabSharing\Nick\Profilometry Analysis\10x_Gel_Dot_200305.csv';
profilometry = csvread(fp, 19, 1);
resolutions_xy = round(dlmread(fp, ',', [3 1 4 1]) / 1000, 4); % Convert to mm
if resolutions_xy(1) == resolutions_xy(2); resolutions_xy = resolutions_xy(1); end
if sum(profilometry(:,end)) == 0; profilometry(:,end) = []; end

% 2D FFT
fs = 1 ./ resolutions_xy;
nyq_f = fs / 2;

% In the case that the profile isn't square get the appropriate segment length for each axis
f_L_row = 2^(nextpow2(size(profilometry,1)) - 1);
f_ax_row = fs*(0:(f_L_row/2))/f_L_row;
f_ax_row_2s = [f_ax_row,fliplr(f_ax_row(2:end-1))];

f_L_col = 2^(nextpow2(size(profilometry,2)) - 1);
f_ax_col = fs*(0:(f_L_col/2))/f_L_col;
f_ax_col_2s = [f_ax_col,fliplr(f_ax_col(2:end-1))]; 

% Truncate to axis lengths and fft
prof_trun = profilometry(1:f_L_row, 1:f_L_col);
prof_trun_fft = fft2(prof_trun);
prof_trun_fft_mag = abs(prof_trun_fft);
prof_trun_fft_mag_shifted = fftshift(prof_trun_fft_mag);

figure;
imagesc(linspace(-nyq_f,nyq_f, size(prof_trun,1)), linspace(-nyq_f,nyq_f, size(prof_trun,2)), prof_trun_fft_mag_shifted)
xlim([-2 2]); ylim([-2 2])

% Masking
f_ax_row_2d = repmat(f_ax_row_2s',[1,size(prof_trun,2)]);
f_ax_col_2d = repmat(f_ax_col_2s,[size(prof_trun,1),1]);
f_t = 0.05;
f_t_logit = f_ax_row_2d < f_t & f_ax_col_2d < f_t;
prof_trun_fft_filt = prof_trun_fft;
prof_trun_fft_filt(f_t_logit) = 0 + 1i;
prof_trun_filt = abs(ifft2(prof_trun_fft_filt));