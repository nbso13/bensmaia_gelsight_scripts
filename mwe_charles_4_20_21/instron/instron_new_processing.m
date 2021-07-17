%% New Instron Processing
% Read and Process files from ratio test in july

close all
clear
addpath('helper_functions');
load('colorscheme.mat');

%%

%data structure - cell array of structs. Each struct is one gel. Each
%struct has a field "traces" with a cell array of trace recordings. Each cell has a nx2 array
%(mm and forces). For each entry in that cell array there's an entry in the
%array in the 'dates' field of the 
gel_structs = {};
date_made = 210709;
date_measured = 210715;
for i = 1:10
    if i > 5
        x = 2;
    else
        x = 1;
    end
    gel_structs{1,i} = struct;
    gel_id_str = strcat('B', num2str(x), '_', num2str(i));
    
    gel_structs{1,i}.batch_num = x;
    gel_structs{1,i}.id = i;
    gel_structs{2,i} = gel_id_str;
    gel_structs{1,i}.date_made = date_made;
    gel_structs{1,i}.dates_measured = [];
    gel_structs{1,i}.dates_measured = [gel_structs{1,i}.dates_measured, date_measured];
    filename = strcat(num2str(date_measured), '_', gel_id_str, '.xls');
    cd '../../../mwe_data/instron_new_data'
    gel_structs{1,i}.data = {};
    gel_structs{1,i}.data{length(gel_structs{1,i}.data)+1} = xlsread(filename, gel_id_str, 'A:B');
    cd ../../bensmaia_gelsight_scripts/mwe_charles_4_20_21/instron
    
end
%% examine data
figure;
hold on
for i = 1:10
    temp_data = gel_structs{1,i}.data{1};
    if gel_structs{1,i}.batch_num == 1
        plot(temp_data(:, 1), temp_data(:,2), 'r')
    else
        plot(temp_data(:, 1), temp_data(:,2), 'b')
    end
end

%% determine contact
initial_indentation = 0.05; %mm
mult = 3; %sd's to exceed
resolution_x = linspace(0, 2, 2000);

figure;
hold on
for i = 1:10
    
    temp_data = gel_structs{1,i}.data{1};
    % get average force in first n = 0.05 mm
    initial = temp_data((temp_data(:,1) < initial_indentation), :);
    mean_force = mean(initial(:, 2));
    sd_force = std(initial(:,2));
    inds_more = find(temp_data(:, 2) > (mean_force + mult*sd_force));
    ind_contact = inds_more(1); %first reading above that value
    
    processed_data = temp_data(ind_contact:end, :);
    
    processed_data = processed_data - processed_data(1,:); %begin at 0 everything
    
    %only examine 2mm indentation
    inds_more = find(processed_data(:, 1) > 2);
    if ~isempty(inds_more)
        two_mm_ind = inds_more(1);
        processed_data = processed_data(1:two_mm_ind, :);
    end
    
    %interpolate
    [force, indentation] = instronInterp(processed_data(:,2) , processed_data(:,1), resolution_x);
    processed_data = [indentation', force'];
    
    gel_structs{1,i}.data{1} = processed_data;
    
    
    if i == 1 || i == 6
        hand_vis = 'on';
    else
        hand_vis = 'off';
    end
    if gel_structs{1,i}.batch_num == 1
        plot(processed_data(:, 1), processed_data(:,2), 'Color', colorscheme(3,:), 'LineWidth', ...
        1.5, 'HandleVisibility', hand_vis)
    elseif gel_structs{1,i}.batch_num == 2
        plot(processed_data(:, 1), processed_data(:,2), 'Color', colorscheme(4,:), 'LineWidth', ...
        1.5, 'HandleVisibility', hand_vis)
    end
    
end

ylabel("Force (N)")
xlabel("indentation (mm)")
title("Stiffness")

strs = ["Batch 1", "Batch 2"]';
colors = colorscheme(3:4, :);
leg = legend([color_legend(strs, colors)]);
leg.Box = 0;

figure;
hold on
for i = 1:10
    processed_data = gel_structs{1,i}.data{1};
    
    
    if i == 1 || i == 6
        hand_vis = 'on';
    else
        hand_vis = 'off';
    end
    if gel_structs{1,i}.batch_num == 1
        plot(processed_data(:, 1), processed_data(:,2), 'Color', colorscheme(3,:), 'LineWidth', ...
        1.5, 'HandleVisibility', hand_vis)
    elseif gel_structs{1,i}.batch_num == 2
        plot(processed_data(:, 1), processed_data(:,2), 'Color', colorscheme(4,:), 'LineWidth', ...
        1.5, 'HandleVisibility', hand_vis)
    end
end
ylabel("Force (N)")
xlabel("indentation (mm)")
title("Stiffness")

strs = ["Batch 1", "Batch 2"]';
colors = colorscheme(3:4, :);
leg = legend([color_legend(strs, colors)]);
leg.Box = 0;

        
        