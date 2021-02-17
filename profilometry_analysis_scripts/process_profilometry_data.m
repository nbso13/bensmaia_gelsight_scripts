%% Processing Data

%% Load data
file_names = {"210212_3_mm_grating_no_gel", "210216_3mm_grating_gel_11", "210216_3mm_grating_gel_7", "210216_wool_blend_no_gel"};
gel_id = [0, 1, 1, 0]; %1 if gel 0 if nah
for index = 1:length(file_names)
    clearvars -except file_names gel_id index
    
    %macros
    DETILT = 0;
    LEXT = 1;
    GEL = gel_id(index);
    file_name = file_names{index};
    %% THIS IS FOR LEXT FILE
    cd ../../csv_data;
    
    if LEXT == 1
        file_list = dir;
        starting_row = 19;
        x_res_row = 3;
        y_res_row = 4;
        z_res_row = 5;
        res_col = 1;
        title_str = file_name;
        
        %picking out pertainent files
        for file = 1:size(file_list,1)
            file_logit(file) = (contains(file_list(file).name, 'csv'))&(contains(file_list(file).name, file_name));
        end
        target_file = file_list(logical(file_logit));
        
        %reading data
        temp_data  = csvread(target_file(1).name, starting_row, 1);
        %read in resolutions from spreadsheet, put them into resolution cell -
        %x, y, z.
        temp_data_res = csvread(target_file(1).name, x_res_row, res_col, [x_res_row res_col z_res_row res_col]);
        
        %samp_freq = 1/(x_res/1000); % once every 2.5 microns
        temp_data_filtered = temp_data;%filter2(fir1(10,0.6), temp_data);
        
        x_res = temp_data_res(1,1)/1000;
        y_res = temp_data_res(2,1)/1000;
        z_res = temp_data_res(3,1)/1000; % all in mms
        x_axis = linspace(0, x_res*(size(temp_data,1)-1), size(temp_data_filtered,1));
        y_axis = linspace(0, y_res*(size(temp_data,2)-1), size(temp_data_filtered,2));
        new_window = temp_data_filtered;
        new_window = new_window./1000;
    end
    %% THIS IS FOR GWEN FILE
    if LEXT == 0
        clear;
        close all;
        file_list = dir;
        starting_row = 5;
        for file = 1:size(file_list,1)
            file_logit(file) = (contains(file_list(file).name, 'csv'));
        end
        target_file = file_list(logical(file_logit));
        
        %reading data
        temp_data  = csvread(target_file(1).name, starting_row, 1);
        temp_data = (temp_data').*1000; %now in mms
        new_window = temp_data;
        load(file_name)
        x_res = gel_dot_200524.x_res;
        y_res = gel_dot_200524.x_res;
        z_res = gel_dot_200524.x_res;
        
        x_axis = linspace(0, x_res*(size(temp_data,1)-1)/1000, size(temp_data,1));
        y_axis = linspace(0, y_res*(size(temp_data,2)-1)/1000, size(temp_data,2));
    end
    
    cd ../bensmaia_gelsight_scripts/profilometry_analysis_scripts
    
    
    %% turning csv data into 3xn matrix of x,y,z points
    y_axis = y_axis'; x_axis = x_axis';
    y_size = size(y_axis, 1); x_size = size(x_axis, 1);
    %only do this if we are detilting.
    if DETILT
        N = zeros(y_size*x_size, 3);
        count = 1;
        reverseStr = '';
        %go through height map and make a 3 column matrix with 3D points for each
        %point
        for x_val = 1:x_size
            
            x_amount = x_axis(x_val);
            msg = sprintf('    fitting line %d of %d', count, x_size);
            fprintf([reverseStr, msg]);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));
            count = count + 1;
            for y_val = 1:y_size
                matrix_val = (x_val-1)*y_size+y_val;
                N(matrix_val, 1) = x_amount;
                N(matrix_val, 2) = y_axis(y_val);
                N(matrix_val, 3) = new_window(x_val, y_val);
            end
        end
        fprintf('\n');
        %calculate plane and fit
        [plane, gof] = fit([N(:,1), N(:,2)], N(:,3), 'poly11');
        params = coeffvalues(plane);
        
        %% Filling in New Plane and Subtracting
        grid_form_plane = zeros(size(new_window));
        %fill a new grid height map with the trend plane
        count = 1;
        reverseStr = '';
        for x_val = 1:x_size
            msg = sprintf('    fitting line %d of %d', count, x_size);
            fprintf([reverseStr, msg]);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));
            count = count + 1;
            for y_val = 1:y_size
                grid_form_plane(x_val, y_val) = x_axis(x_val)*params(2) + y_axis(y_val)*params(3);
            end
        end
        fprintf('\n');
    end
    
    figure;
    [l_lim, u_lim] = bounds(new_window(:));
    imagesc(x_axis', y_axis', new_window')
    c = colorbar;
    ylabel(c, 'mm');
    title([title_str, " Before Processing"]);
    xlabel('mm'); ylabel('mm');
    caxis([l_lim, u_lim]); %ylim([1240 1270]); xlim([2810 2930])
    
    if DETILT
        new_window = new_window - grid_form_plane; % Subtract out
    end
    y_axis = y_axis'; x_axis = x_axis';
    
    [l_lim, u_lim] = bounds(new_window(:));
    figure;
    imagesc(x_axis, y_axis, new_window')
    c = colorbar;
    ylabel(c, 'mm');
    title([title_str, " After DeTilt"]);
    xlabel('mm'); ylabel('mm');
    caxis([l_lim, u_lim]); %ylim([1240 1270]); xlim([2810 2930])
    
    %
    % figure;
    % hist(new_window)
    % if strcmp(file_name, '201118_corduroy_35_gel')
    %     new_window = new_window + 0.002;
    %     new_window(new_window<0) = 0;
    % end
    % if strcmp(file_name, '201118_corduroy_no_gel')
    %     disp("corduroy no gel");
    %     new_window(new_window<0.1) = 0;
    % end
    
    gel_dot_200305 = struct;
    gel_dot_200305.profile = new_window';
    gel_dot_200305.x_res = x_res;
    gel_dot_200305.y_res = y_res;
    gel_dot_200305.z_res = z_res;
    gel_dot_200305.x_axis = x_axis';
    gel_dot_200305.y_axis = y_axis';
    gel_dot_200305.type = "dots";
    filename = strcat(file_name, "_processed.mat");
    visualizeProfile(gel_dot_200305);
    cd ../../mat_files
    if GEL
        gel = gel_dot_200305;
        save(filename, "gel");
    else
        no_gel = gel_dot_200305;
        save(filename, "no_gel");
    end
    
    cd ../bensmaia_gelsight_scripts/profilometry_analysis_scripts
end

%% processsing - crop and rotate if necessary
%
% crop_vector_dots = [300 300 4300 4300];
% crop_vector_sine = [0 800 5000 5000];
% crop_vector_white = [400 0 6000 4400];
% crop_vector_uphol = [0 0 4200 4200];
% crop_vector_uphol2 = [0 0 4600 4600];
% crop_vector_rope = [0 0 4400 4400];
% rotate_amount = -4; %degrees
%
% processed_out = imcrop(imrotate(new_window, rotate_amount), crop_vector_dots);
% x_axis_new = x_axis(1:4301);
% x_axis = 0;
% x_axis = x_axis_new;
% processed_out = processed_out./1000;
% x_axis = x_axis/1000;
% y_axis = y_axis/1000;
%
% texture = strsplit(file_name, ' ');
%
% min_sine = 0.3;
% min_dots = 0.38;
% max_sine = 0.55;
% max_dots = 0.7;
% max_white = 0.65;
% min_white = 0.36;
% min_uphol = 0.2;
% max_uphol = 1.1;
% min_uphol_g = 0.25;
% max_uphol_g = 0.45;
% min_uphol2 = 0;
% max_uphol2 = 3;
% min_uphol2_g = 0.37;
% max_uphol2_g = 0.6;
% min_rope = 0;
% max_rope = 2;
% min_rope_g = 0.05;
% max_rope_g = 0.25;
% if(contains(texture{1}, 'rope'))
%     if(contains(texture{2}, 'gelsight'))
%         mini = min_rope_g;
%         maxi = max_rope_g;
%     else
%         mini = min_rope;
%         maxi = max_rope;
%     end
% elseif(contains(texture{1}, 'upholstry2'))
%     if(contains(texture{2}, 'gelsight'))
%         mini = min_uphol2_g;
%         maxi = max_uphol2_g;
%     else
%         mini = min_uphol2;
%         maxi = max_uphol2;
%     end
% elseif(contains(texture{1}, 'upholstry'))
%     if(contains(texture{2}, 'gelsight'))
%         mini = min_uphol_g;
%         maxi = max_uphol_g;
%     else
%         mini = min_uphol;
%         maxi = max_uphol;
%     end
% elseif(contains(texture{1}, 'whitenoise'))
%      mini = min_white;
%      maxi = max_white;
% elseif(contains(texture{1}, 'sine'))
%     mini = 0;
%     maxi = 10000;
% elseif(contains(texture{1}, 'dots'))
%     mini = min_dots;
%     maxi = max_dots;
% end
% processed_out(processed_out<mini) = mini;
% processed_out(processed_out>maxi) = maxi;


% %% write to a csv
% temp_data_res = temp_data_res./1000;
% new_name = strcat(file_name, name_addition);
% dlmwrite(new_name, temp_data_res);
% dlmwrite(new_name, processed_out, '-append');

%% plotting
% [l_lim, u_lim] = bounds(processed_out(:));
% imagesc(x_axis, y_axis, processed_out')
% c = colorbar;
% ylabel(c, 'mm');
% title(title_str);
% xlabel('mm'); ylabel('mm');
% caxis([l_lim, u_lim]); %ylim([1240 1270]); xlim([2810 2930])


%% Plot Cross Section

% processed_out = processed_out';
% line = processed_out(2000, :);
% line = smoothdata(line);
% plot(x_axis, line, 'LineWidth', 1.5);
% title('Row 1000 of Dots Profile');
% xlabel('mm');
% ylabel('mm');



%
% % Frequency analysis
% samp_freq = 1/(x_res/1000); % once every 2.5 microns
% samp_period = 1/samp_freq;
% L = length(new_window);
% lines_total = size(new_window,2);
% f_ax = samp_freq*(0:(L/2))/L;
%
% reverseStr = '';
% count = 1;
% [size_x, size_y] = size(new_window);
% transform = zeros(size_x, size_y);
% for line = 1:lines_total
%     temp_line = new_window(:,line);%(1:(end-1));
%     %subtract out mean
%     temp_line = temp_line - mean(temp_line);
%     %do fft and normalize by period
%     temp_line = fft(temp_line)*samp_period;
%     %magnitude of complex components
%     two_sidedfft = abs(temp_line);
%     %spectrum is symmetrical around nyquist freq so only first half
%     one_sidedfft = two_sidedfft(1:floor(L/2+1));
%     %double the magnitude since we only took the first half
%     one_sidedfft(2:end-1) = 2*one_sidedfft(2:end-1);
%     transform(line, :) = one_sidedfft;
%     msg = sprintf('    Analyzing trace %d of %d', count, lines_total);
%     fprintf([reverseStr, msg]);
%     reverseStr = repmat(sprintf('\b'), 1, length(msg));
%     count = count + 1;
% end
% fprintf('\n')
% fourier_average = mean(transform, 1);
% plot(f_ax, fourier_average)
%
% title('Amplitude Spectrum of Profilometry Data')
% xlabel('frequency (1/mm)')
% ylabel('Amplitude')
%
%
%
% %
% % test = temp_data{file} - temp_data_filtered{file};
% % [l_lim, u_lim] = bounds(test(:));
% % imagesc(test);
% % % %
% % [l_lim, u_lim] = bounds(temp_data{2}(:));
% % subplot(2,2,3);
% % imagesc(temp_data{2}')
% % caxis([l_lim, u_lim]); %ylim([1240 1270]); xlim([2810 2930])
% % subplot(2,2,4);
% % imagesc(temp_data_filtered{2}')
%
%
% %%
% %daspect([1 1 1]);
% %surf((temp_data{1})', 'LineStyle', 'none');
