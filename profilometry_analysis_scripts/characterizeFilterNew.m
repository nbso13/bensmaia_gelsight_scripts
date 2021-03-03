function [amp_ratio_vec, amp_fy, fig_handles] = characterizeFilterNew(data_unfilt, data_filt, ...
    vert_line, horiz_line, plotflag,  one_dim)
%characterizeFilterFull visualizes data filtered and unfiltered, conducts
%fourier analysis and produces frequency response. vertline and horiz line
%specify where to draw sample traces for visualization verically and
%horizontally from the data. c_ax specifies the colormap range (either
%"max" or a value) and log specificies if log scale (1 if yes).
% max_freq gives the maximum frequency to be displayed (i.e., 4 dots per
% mm) variable "one dim" indicates this is a regular texture in one
% dimension, like a grating. if it is one dim, second number is which dim
% (1 is down, 2 is across)

%max value in data
max_filt = max(data_filt.profile, [], 'all');
max_unfilt = max(data_unfilt.profile, [], 'all');
max_val = max([max_filt, max_unfilt]);
filt_color = 'r';
unfilt_color = 'b';
fig_handles = [];
%% plot in subplot array A, B, A-B
if plotflag
    fig_handles(1) = figure;
    subplot(3,3,2);
    title_str = "Filtered Data";
    imagesc(data_filt.x_axis, data_filt.y_axis, data_filt.profile)
    if ~one_dim
        xline(data_filt.x_axis(vert_line));
        yline(data_filt.y_axis(horiz_line));
    end
    c = colorbar;
    ylabel(c, 'mm');
    caxis([0, max_unfilt]);
    title(title_str);
    xlabel('mm'); ylabel('mm');
    
    subplot(3,3,1);
    title_str = "Unfiltered Data";
    imagesc(data_unfilt.x_axis, data_unfilt.y_axis, data_unfilt.profile)
    if ~one_dim
        xline(data_unfilt.x_axis(vert_line));
        yline(data_unfilt.y_axis(horiz_line));
    end
    
    c = colorbar;
    ylabel(c, 'mm');
    caxis([0, max_unfilt]);
    title(title_str);
    xlabel('mm'); ylabel('mm');
    
    
    %align and calculate difference
    subplot(3,3,3);
    prof_diff = max_unfilt-max_filt;
    filt_to_subtract = data_unfilt.profile(1:size(data_filt.profile, 1), 1:size(data_filt.profile,2));
    diff = data_filt.profile+prof_diff - filt_to_subtract;
    title_str = "Difference, max aligned";
    imagesc(data_filt.x_axis, data_filt.y_axis, diff)
    c = colorbar;
    ylabel(c, 'mm');
    caxis([0, max_val]);
    title(title_str);
    xlabel('mm'); ylabel('mm');
    
    %% pick out vertical trace
    data_filt_vline = data_filt.profile(:, vert_line);
    data_unfilt_vline = data_unfilt.profile(:,vert_line);
    common_y_axis = data_unfilt.y_axis;
    
    if one_dim
        if length(one_dim) < 2
            error("one_dim var needs to be a size 2 array - first index indicates one dim and second indicates direction, 1 is horizontal 2 is vertical");
        end
        if one_dim(2) == 2 % vertical grating averaging (right to left)
            data_filt_vline = mean(data_filt.profile, 2);
            data_unfilt_vline = mean(data_unfilt.profile, 2);
        end
    end
    
    subplot(3,3,5);
    line = data_filt_vline;
    plot(common_y_axis, line, filt_color)
    title("V Trace from Filtered")
    xlabel("mm")
    ylabel("mm")
    
    subplot(3,3,4);
    line = data_unfilt_vline;
    plot(common_y_axis, line, unfilt_color)
    title("V Trace from Unfiltered")
    xlabel("mm")
    ylabel("mm")
    
    %adjust max
    data_filt_vmax = max(data_filt_vline);
    data_unfilt_vmax = max(data_unfilt_vline);
    vdifference = data_unfilt_vmax - data_filt_vmax;
    
    data_filt_vline = data_filt_vline+vdifference;
    subplot(3,3,6);
    plot(common_y_axis, data_filt_vline, filt_color)
    hold on;
    plot(common_y_axis, data_unfilt_vline, unfilt_color)
    title("V Filt and Unfilt Traces Adjusted")
    xlabel("mm")
    ylabel("mm")
    
    %% pick out horizontal trace
    data_filt_hline = data_filt.profile(horiz_line, :);
    data_unfilt_hline = data_unfilt.profile(horiz_line, :);
    common_x_axis = data_unfilt.x_axis;
    
    if one_dim
        if one_dim(2) == 1 % vertical grating averaging (up and down)
            data_filt_hline = mean(data_filt.profile, 1);
            data_unfilt_hline = mean(data_unfilt.profile, 1);
        end
    end
    
    subplot(3,3,8);
    line = data_filt_hline;
    plot(common_x_axis, line, filt_color)
    title("H Trace from Filtered")
    xlabel("mm")
    ylabel("mm")
    
    subplot(3,3,7);
    line = data_unfilt_hline;
    plot(common_x_axis, line, unfilt_color)
    title("H Trace from Unfiltered")
    xlabel("mm")
    ylabel("mm")
    
    %adjust max
    data_filt_hmax = max(data_filt_hline);
    data_unfilt_hmax = max(data_unfilt_hline);
    hdifference = data_unfilt_hmax - data_filt_hmax;
    
    data_filt_hline = data_filt_hline+hdifference;
    subplot(3,3,9);
    plot(common_x_axis, data_filt_hline, filt_color)
    hold on;
    plot(common_x_axis, data_unfilt_hline, unfilt_color)
    title("H Filt and Unfilt Traces Adjusted")
    xlabel("mm")
    ylabel("mm")
    set(gcf, 'position', [100 100 1100 800]);
end

% freq analysis on sim
% max_freq_x=max_freq;
% max_freq_y =max_freq;
disp("Fourier transform down columns. Please ensure scanning direction is down.")
[unfilt_amp_vec, unfilt_fy] = freqAnalysisFullNew(data_unfilt);
[filt_amp_vec, filt_fy] = freqAnalysisFullNew(data_filt);

%limiting viewing so that we don't have to see super high freqs we don't
%care about
% view_x_filt = filt_fx(filt_fx<max_freq_x);
% view_y_filt = filt_fy(filt_fy<max_freq_y);
% view_x_unfilt = unfilt_fx(unfilt_fx<max_freq_x);
% view_y_unfilt = unfilt_fy(unfilt_fy<max_freq_y);
% view_unfilt_mat = unfilt_amp_mat(1:length(view_y_unfilt), 1:length(view_x_unfilt));
% view_filt_mat = filt_amp_mat(1:length(view_y_filt), 1:length(view_x_filt));

amp_ratio_vec = filt_amp_vec./unfilt_amp_vec;
amp_fy = unfilt_fy;

% %finding max between them to set colorbar sclae
% maxo = max([max(filt_amp_mat, [], 'all') max(unfilt_amp_mat, [], 'all')]);

if plotflag
    fig_handles(2) = figure;
    subplot(2,2,1)
    plot(unfilt_fy,unfilt_amp_vec);
    title("Amplitude Spectrum, Unfiltered")
    ylabel("Amplitude")
    xlabel("Frequency (1/mm)")
    ax = gca;
    ax.YDir = 'normal';
    set(gcf, 'position', [100 100 1100 800]);
%     colorbar;
%     if isstring(cax)
%         if cax == "max"
%             caxis([0 maxo]);
%         else
%             error("Enter a proper value for cax limit.");
%         end
%     else
%         caxis([0 cax]);
%     end
%     if log == 1
%         set(gca,'ColorScale','log')
%     end
    
    subplot(2,2,2)
    plot(filt_fy, filt_amp_vec);
    title("Amplitude Spectrum, Filtered")
    ylabel("amplitude")
    xlabel("Frequency (1/mm)")
    ax = gca;
    ax.YDir = 'normal';
%     
%     colorbar;
%     set(gcf, 'position', [100 100 1100 800]);
%     if isstring(cax)
%         if cax == "max"
%             caxis([0 maxo]);
%         else
%             error("Enter a proper value for cax limit.");
%         end
%     else
%         caxis([0 cax]);
%     end
%     if log == 1
%         set(gca,'ColorScale','log')
%     end
%     
    
    %max_amp = max(amp_ratio_mat, [], 'all');
    subplot(2,2,3)
    plot(unfilt_fy, amp_ratio_vec)
    title("Amplitude Ratio For Filt vs Unfilt")
    ylabel("Ratio, Filtered to Unfiltered")
    xlabel("Frequency (1/mm)")
    ax = gca;
    ax.YDir = 'normal';
%     colorbar;
%     caxis([0 1]);
    set(gcf, 'position', [100 100 1100 800]);
end

% %crop mat
% s_un_x = length(unfilt_fx);
% s_un_y = length(unfilt_fy);
% s_f_x = length(filt_fx);
% s_f_y = length(filt_fy);
% s_x_min = min(s_un_x, s_f_x);
% s_y_min = min(s_un_y, s_f_y);
% 
% cropped_unfilt_amp = unfilt_amp_mat(1:s_y_min, 1:s_x_min);
% cropped_filt_amp = filt_amp_mat(1:s_y_min, 1:s_x_min);
% cropped_unfx = unfilt_fx(1:s_x_min);
% cropped_unfy = unfilt_fy(1:s_y_min);
% cropped_fx = filt_fx(1:s_x_min);
% cropped_fy = filt_fy(1:s_y_min);
% 
% % interpolate unfilt amp_mat to filt frequencies to calculate ratio
% cropped_unfilt_amp_interp = interp2(cropped_unfx, cropped_unfy, cropped_unfilt_amp, cropped_fx, cropped_fy);
% amp_ratio_mat = cropped_unfilt_amp_interp./cropped_filt_amp; %calculate amp ratio
% amp_fx = cropped_fx;
% amp_fy = cropped_fy;

% 
% figure;
% surf(view_x_unfilt, view_y_unfilt, view_unfilt_mat);
% title("Amplitude Spectrum, Unfiltered")
% ylabel("Frequency (1/mm)")
% xlabel("Frequency (1/mm)")
% set(gca, 'ZScale', 'log');
% 
% figure;
% surf(view_x_filt, view_y_filt, view_filt_mat);
% title("Amplitude Spectrum, Filtered")
% ylabel("Frequency (1/mm)")
% xlabel("Frequency (1/mm)")
% set(gca, 'ZScale', 'log');


end

