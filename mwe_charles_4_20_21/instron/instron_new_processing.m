%% New Instron Processing
% Read and Process files from ratio test in july

close all
clear
addpath('helper_functions');
load('colorscheme.mat');

%% 210715

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
trace_number = 1;
[fig] = plot_all_together(gel_structs, trace_number, colorscheme);

%% determine contact

trace_number = 1;
for i = 1:10
    processed_data = contact_process(i, trace_number, gel_structs);
    gel_structs{1,i}.data{trace_number} = processed_data;
end

[fig] = plot_all_together(gel_structs, trace_number, colorscheme);
title("210715")

%% 210720 
date_measured = 210720;
for i = 1:10
    if i > 5
        x = 2;
        gel_id_str = strcat('B', num2str(x), '_', num2str(i));
        gel_structs{1,i}.dates_measured = [gel_structs{1,i}.dates_measured, date_measured];
        filename = strcat(num2str(date_measured), '_', gel_id_str, '.xls');
        cd '../../../mwe_data/instron_new_data'
        gel_structs{1,i}.data{length(gel_structs{1,i}.data)+1} = xlsread(filename, gel_id_str, 'A:B');
        cd ../../bensmaia_gelsight_scripts/mwe_charles_4_20_21/instron
    else
        x = 1;
        for trial = 1:3
            gel_id_str = strcat('B', num2str(x), '_', num2str(i), '_', num2str(trial));
            gel_structs{1,i}.dates_measured = [gel_structs{1,i}.dates_measured, date_measured];
            filename = strcat(num2str(date_measured), '_', gel_id_str, '.xls');
            cd '../../../mwe_data/instron_new_data'
            gel_structs{1,i}.data{length(gel_structs{1,i}.data)+1} = xlsread(filename, gel_id_str, 'A:B');
            cd ../../bensmaia_gelsight_scripts/mwe_charles_4_20_21/instron
        end
    end
end
% examine data
figure;
hold on
for i = 1:10
    if i < 6
        for j = 3:4
            temp_data = gel_structs{1,i}.data{j};
            if gel_structs{1,i}.batch_num == 1
                plot(temp_data(:, 1), temp_data(:,2), 'r')
            else
                plot(temp_data(:, 1), temp_data(:,2), 'b')
            end
        end
    end
    temp_data = gel_structs{1,i}.data{2};
    if gel_structs{1,i}.batch_num == 1
        plot(temp_data(:, 1), temp_data(:,2), 'r')
    else
        plot(temp_data(:, 1), temp_data(:,2), 'b')
    end
end

figure;
hold on
for i = 1:10
    if i < 6
        for trace_number = 3:4
            processed_data = contact_process(i, trace_number, gel_structs);
            gel_structs{1,i}.data{trace_number} = processed_data;
            hand_vis = 'off';
            plot(processed_data(:, 1), processed_data(:,2), 'Color', colorscheme(3,:), 'LineWidth', ...
                1.5, 'HandleVisibility', hand_vis)
        end
    end
    trace_number = 2;
    processed_data = contact_process(i, trace_number, gel_structs);
    gel_structs{1,i}.data{trace_number} = processed_data;
    
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
title("210720")

strs = ["Batch 1", "Batch 2"]';
colors = colorscheme(3:4, :);
leg = legend([color_legend(strs, colors)]);
leg.Box = 0;



%% 210721
date_measured = 210721;
for i = 1:10
    if i > 5
        x = 2;
        
    else
        x = 1;
    end
    gel_id_str = strcat('B', num2str(x), '_', num2str(i));
    gel_structs{1,i}.dates_measured = [gel_structs{1,i}.dates_measured, date_measured];
    filename = strcat(num2str(date_measured), '_', gel_id_str, '.xls');
    cd '../../../mwe_data/instron_new_data'
    gel_structs{1,i}.data{length(gel_structs{1,i}.data)+1} = xlsread(filename, gel_id_str, 'A:B');
    cd ../../bensmaia_gelsight_scripts/mwe_charles_4_20_21/instron
end
% examine data
[fig] = plot_all_together(gel_structs, "last", colorscheme);

figure;
hold on
for i = 1:10
    trace_number = "last";
    processed_data = contact_process(i, trace_number, gel_structs);
    gel_structs{1,i}.data{length(gel_structs{1,i}.data)} = processed_data;
end

[fig] = plot_all_together(gel_structs, "last", colorscheme);
title("210721")

%% processing together
% examine changes
mean_percent_difference = [];
for i = 1:10
    if i <6
        trial_two = gel_structs{1, i}.data{2};
        trial_three= gel_structs{1, i}.data{3};
        trial_four = gel_structs{1, i}.data{4};
        second_trace = mean([trial_two(:,2), trial_three(:,2), trial_four(:,2)], 2);
%         second_trace = trial_two(:,2);
    else
        second_trace = gel_structs{1, i}.data{2};
        second_trace = second_trace(:,2);
    end
    first_trace = gel_structs{1, i}.data{1};
    mean_percent_difference(i) = nanmean(first_trace(:, 2)./second_trace);
end


for i = 1:10
    plot_all_measurements(gel_structs, i, colorscheme);
    title(strcat("Gel ", num2str(i)));
end

%% Fit to power law and to polynomial and try to predict change in coefficients with change in time.
fs = {};
ps = cell(10);
for i = 1:10
    for j = 1:length(gel_structs{1, i}.data)
        fs{i,j} = fit(gel_structs{1, i}.data{j}(:, 1), gel_structs{1, i}.data{j}(:, 2), 'poly2');
        ps{1, i} = [ps{1, i}, fs{i, j}.p1];
        ps{2, i} = [ps{2, i}, fs{i, j}.p2];
        ps{3, i} = [ps{3, i}, fs{i, j}.p3];
    end
end



%plot p1s over time
figure;
hold on
for i = 1:10
    dates = gel_structs{1,i}.dates_measured;
    days_later = subtract_dates(gel_structs{1, i}.date_made, dates);
    if i < 6
        plot(days_later, ps{1,i}, 'Color', colorscheme(3,:));
    else
        plot(days_later, ps{1,i}, 'Color', colorscheme(4,:));
    end
end

%plot p2s over time
figure;
hold on
for i = 1:10
    dates = gel_structs{1,i}.dates_measured;
    days_later = subtract_dates(gel_structs{1, i}.date_made, dates);
    if i < 6
        plot(days_later, ps{2,i}, 'Color', colorscheme(3,:));
    else
        plot(days_later, ps{2,i}, 'Color', colorscheme(4,:));
    end
end

%plot p3s over time
figure;
hold on
for i = 1:10
    dates = gel_structs{1,i}.dates_measured;
    days_later = subtract_dates(gel_structs{1, i}.date_made, dates);
    if i < 6
        plot(days_later, ps{3,i}, 'Color', colorscheme(3,:));
    else
        plot(days_later, ps{3,i}, 'Color', colorscheme(4,:));
    end
end

%% functions

function [processed_data] = contact_process(gel_ID, trace_num, gel_structs)
if isstring(trace_num) && strcmp(trace_num, "last")
    temp_data = gel_structs{1,gel_ID}.data{length(gel_structs{1,gel_ID}.data)};
else
    temp_data = gel_structs{1,gel_ID}.data{trace_num};
end
initial_indentation = 0.05; %mm
mult = 3; %sd's to exceed
resolution_x = linspace(0, 2, 2000);
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
force(isnan(force)) = 0;
indentation(isnan(indentation)) = 0;
processed_data = [indentation', force'];
end


function [fig] = plot_all_measurements(gel_structs, gel_ID, colorscheme)
fig = figure;
hold on

% get unique measurement days
dates = gel_structs{1,gel_ID}.dates_measured;
unique_dates = unique(dates, 'stable');
leg_flags = ones(1, length(unique_dates));
days_later = subtract_dates(gel_structs{1, gel_ID}.date_made, unique_dates);
strs = [];
for i = 1:length(days_later)
    strs = [strs, string(strcat(num2str(days_later(i)), " days later"))];
end
    
data = gel_structs{1,gel_ID}.data;

for i = 1:length(data)
    temp = find(unique_dates == dates(i)); ind = temp(1);
    if leg_flags(ind)
        leg_flags(ind) = 0;
        hand_vis = 'on';
    else
        hand_vis = 'off';
    end
    plot(data{i}(:,1), data{i}(:,2), 'Color', colorscheme(2+ind,:), 'LineWidth', ...
        1.5, 'HandleVisibility', hand_vis)
end

colors = colorscheme(3:2+length(unique_dates), :);
leg = legend([color_legend(strs', colors)]);
leg.Box = 0;

end

function [days_later] = subtract_dates(date_ref, dates_compare)
    year = floor(date_ref/10000) + 2000;
    month = floor((date_ref - floor(date_ref/10000)*10000)/100);
    day = date_ref - floor(date_ref/100)*100;
    dt_ref = datetime(year, month, day);
    year = floor(dates_compare/10000) + 2000;
    month = floor((dates_compare - floor(dates_compare/10000)*10000)/100);
    day = dates_compare - floor(dates_compare/100)*100;
    dt_compare = datetime(year, month, day);
    days_later = dt_compare - dt_ref;
    days_later.Format = 'd';
    days_later = time2num(days_later);
end

function [fig] = plot_all_together(gel_structs, trace_num, colorscheme)
fig = figure;
hold on
for i = 1:10
    if isstring(trace_num) && strcmp(trace_num, "last")
        processed_data = gel_structs{1,i}.data{length(gel_structs{1,i}.data)};
    else
        processed_data = gel_structs{1,i}.data{trace_num};
    end
    
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
end
    

        