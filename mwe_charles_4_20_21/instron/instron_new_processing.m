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
    temp = readcell(filename, 'Sheet', gel_id_str, 'Range', 'A:B');
    gel_structs{1,i}.data{length(gel_structs{1,i}.data)+1} = cell2mat(temp(5:end, :));
    cd ../../bensmaia_gelsight_scripts/mwe_charles_4_20_21/instron
    
end
%% examine data
trace_number = 1;
hand_vis_on_indices = [1,6];
[fig] = plot_all_together(gel_structs, trace_number, colorscheme, hand_vis_on_indices);

%% determine contact

trace_number = 1;
for i = 1:10
    processed_data = contact_process(i, trace_number, gel_structs);
    gel_structs{1,i}.data{trace_number} = processed_data;
end

[fig] = plot_all_together(gel_structs, trace_number, colorscheme, hand_vis_on_indices);
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
        temp = readcell(filename, 'Sheet', gel_id_str, 'Range', 'A:B');
        gel_structs{1,i}.data{length(gel_structs{1,i}.data)+1} = cell2mat(temp(5:end, :));
        cd ../../bensmaia_gelsight_scripts/mwe_charles_4_20_21/instron
    else
        x = 1;
        for trial = 1:3
            gel_id_str = strcat('B', num2str(x), '_', num2str(i), '_', num2str(trial));
            gel_structs{1,i}.dates_measured = [gel_structs{1,i}.dates_measured, date_measured];
            filename = strcat(num2str(date_measured), '_', gel_id_str, '.xls');
            cd '../../../mwe_data/instron_new_data'
            temp = readcell(filename, 'Sheet', gel_id_str, 'Range', 'A:B');
            gel_structs{1,i}.data{length(gel_structs{1,i}.data)+1} = cell2mat(temp(5:end, :));
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
    temp = readcell(filename, 'Sheet', gel_id_str, 'Range', 'A:B');
    gel_structs{1,i}.data{length(gel_structs{1,i}.data)+1} = cell2mat(temp(5:end, :));
    cd ../../bensmaia_gelsight_scripts/mwe_charles_4_20_21/instron
end
% examine data
[fig] = plot_all_together(gel_structs, "last", colorscheme, hand_vis_on_indices);

figure;
hold on
for i = 1:10
    trace_number = "last";
    processed_data = contact_process(i, trace_number, gel_structs);
    gel_structs{1,i}.data{length(gel_structs{1,i}.data)} = processed_data;
end

[fig] = plot_all_together(gel_structs, "last", colorscheme, hand_vis_on_indices);
title("210721")


%% 210721 NEW GELS
date_made = 210719;
date_measured = 210721;
for i = 11:18
    if i > 14
        x = 4;
    else
        x = 3;
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
    temp = readcell(filename, 'Sheet', gel_id_str, 'Range', 'A:B');
    gel_structs{1,i}.data{length(gel_structs{1,i}.data)+1} = cell2mat(temp(5:end, :));
    cd ../../bensmaia_gelsight_scripts/mwe_charles_4_20_21/instron
    
end
% examine data
trace_number = 1;
new_structs = gel_structs(:, 11:18);
hand_vis_on_indices = [1,6, 11, 15];
[fig] = plot_all_together(new_structs, trace_number, colorscheme, hand_vis_on_indices);

% determine contact

for i = 11:18
    processed_data = contact_process(i, trace_number, gel_structs);
    gel_structs{1,i}.data{trace_number} = processed_data;
end

new_structs = gel_structs(:, 11:18);

[fig] = plot_all_together(new_structs, trace_number, colorscheme, hand_vis_on_indices);
title("after finding contact");


[fig] = plot_all_together(gel_structs, trace_number, colorscheme, hand_vis_on_indices);
title("210721 NEW")


for i = 1:10
    plot_all_measurements(gel_structs, i, colorscheme);
    title(strcat("Gel ", num2str(i)));
end


%% 210728 1:10 again!
file_counter = 1;
date_measured = 210728;
for i = 1:10
    if i == 4
        continue;
    end
    if i > 5
        x = 2;
        
    else
        x = 1;
    end
    gel_id_str = strcat('B', num2str(x), '_', num2str(mod(i,5)+5*~logical(mod(i,5))));
    gel_structs{1,i}.dates_measured = [gel_structs{1,i}.dates_measured, date_measured];
    filename = strcat(num2str(date_measured), '_', "B1B2_old_indentor_", num2str(file_counter), '.xls');
    file_counter = file_counter + 1;
    cd '../../../mwe_data/instron_new_data'
    temp = readcell(filename, 'Sheet', gel_id_str, 'Range', 'A:B');
    gel_structs{1,i}.data{length(gel_structs{1,i}.data)+1} = cell2mat(temp(5:end, :));
    cd ../../bensmaia_gelsight_scripts/mwe_charles_4_20_21/instron
end
% examine data
[fig] = plot_all_together(gel_structs, "last", colorscheme, hand_vis_on_indices);

figure;
hold on
for i = 1:10
    trace_number = "last";
    processed_data = contact_process(i, trace_number, gel_structs);
    gel_structs{1,i}.data{length(gel_structs{1,i}.data)} = processed_data;
end

[fig] = plot_all_together(gel_structs, "last", colorscheme, hand_vis_on_indices);
title("210721")


%% get change in slope per change in time for every gel


figure;
hold on
change_slopes = [];
change_slope_controlled = [];
for i = 1:10
    if i == 4
        continue
    end
    % get unique dates
    dates = gel_structs{1,i}.dates_measured;
    [days_later, unique_dates] = get_unique_days_later(dates, ...
        gel_structs{1,i}.date_made);
    data = gel_structs{1,i}.data;
    slopes = zeros(1,length(unique_dates));
    for j = 1:length(unique_dates)
        temp = find(dates == unique_dates(j), 1);  %find the first measurement for each of the unique dates
        run = data(temp); %run from that date
        slopes(j) = fit_run(run); % get slope
    end
    scatter(days_later, slopes, [], colorscheme(i, :));
    best_fit = fit(days_later', slopes', "poly1");
    av_y = feval(best_fit, days_later); 
    plot(days_later, av_y, 'Color', colorscheme(i, :))
    change_slopes = [change_slopes, best_fit.p1];
end
xlabel("Days Later")
ylabel("Slope")
title("Change in Compliance Slope over Time for 2.8 gels")


disp("Mean change in N/mm stiffness slope per day:")
disp(strcat(num2str(mean(change_slopes))), " +/- ", num2str(std(change_slopes)));








%% 210722 finger
date_made = 210722;
date_measured = 210722;
for i = 1:7
    gel_structs{1, length(gel_structs)+1} = struct;
    gel_id_str = strcat('nick_left_finger', '_', num2str(i));
    gel_structs{1,length(gel_structs)}.batch_num = 0;
    gel_structs{1,length(gel_structs)}.id = i;
    gel_structs{2,length(gel_structs)} = gel_id_str;
    worksheet_name = strcat("Nick_left_index_", num2str(i));
    gel_structs{1,length(gel_structs)}.date_made = date_made;
    gel_structs{1,length(gel_structs)}.dates_measured = [];
    gel_structs{1,length(gel_structs)}.dates_measured = [gel_structs{1,length(gel_structs)}.dates_measured, date_measured];
    filename = strcat(num2str(date_measured), '_', "nick_finger_test_", num2str(i), '.xls');
    cd '../../../mwe_data/instron_new_data'
    gel_structs{1,length(gel_structs)}.data = {};
    gel_structs{1,length(gel_structs)}.data{length(gel_structs{1,...
        length(gel_structs)}.data)+1} = xlsread(filename, worksheet_name, 'A:B');
    cd ../../bensmaia_gelsight_scripts/mwe_charles_4_20_21/instron
end
% examine data
trace_number = 1;
new_structs = gel_structs(:, 19:end);
hand_vis_on_indices = [1];
[fig] = plot_all_together(new_structs, trace_number, colorscheme, hand_vis_on_indices);

% determine contact

for i = 19:25
    processed_data = contact_process(i, trace_number, gel_structs);
    gel_structs{1,i}.data{trace_number} = processed_data;
end

% kick 20 and 23
gel_structs = [gel_structs(:, 1:19), gel_structs(:, 21:22), gel_structs(:, 24:end)];
new_structs = gel_structs(:, 19:end);
[fig] = plot_all_together(new_structs, trace_number, colorscheme,  hand_vis_on_indices);
title("after finding contact");

hand_vis_on_indices = [1, 6, 11, 15, 19];
[fig] = plot_all_together(gel_structs, trace_number, colorscheme, hand_vis_on_indices);
title("finger NEW")

hand_vis_on_indices = [1, 6, 11, 15];
[fig] = plot_all_together(gel_structs(:, 1:18), trace_number, colorscheme, hand_vis_on_indices);
plot_av(gel_structs(:, 19:end), "last", "Nick Left Index Finger");


hand_vis_on_indices = [1, 6, 11, 15];
[fig] = plot_all_together(gel_structs(:, 1:5), trace_number, colorscheme, hand_vis_on_indices);
plot_av(gel_structs(:, 19:end), "last", "Nick Left Index Finger");
title("Batch 1")

[fig] = plot_all_together(gel_structs(:, 6:10), trace_number, colorscheme, hand_vis_on_indices);
plot_av(gel_structs(:, 19:end), "last", "Nick Left Index Finger");
title("Batch 2")

[fig] = plot_all_together(gel_structs(:, 11:14), trace_number, colorscheme, hand_vis_on_indices);
plot_av(gel_structs(:, 19:end), "last", "Nick Left Index Finger");
title("Batch 3")

[fig] = plot_all_together(gel_structs(:, 14:18), trace_number, colorscheme, hand_vis_on_indices);
plot_av(gel_structs(:, 19:end), "last", "Nick Left Index Finger");
title("Batch 4")

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


function [days_later, unique_dates, leg_flags] = get_unique_days_later(dates_measured, date_made)
% get unique measurement days
unique_dates = unique(dates_measured, 'stable');
leg_flags = ones(1, length(unique_dates));
days_later = subtract_dates(date_made, unique_dates);
end

function [processed_data] = contact_process(gel_ID, trace_num, gel_structs)
if isstring(trace_num) && strcmp(trace_num, "last")
    temp_data = gel_structs{1,gel_ID}.data{length(gel_structs{1,gel_ID}.data)};
else
    temp_data = gel_structs{1,gel_ID}.data{trace_num};
end
initial_indentation = 0.05; %mm
mult = 3; %sd's to exceed
resolution_x = linspace(0, 1.5, 1500);
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
[days_later, unique_dates, leg_flags] = get_unique_days_later(gel_structs{1,gel_ID}.dates_measured, ...
    gel_structs{1, gel_ID}.date_made);
strs = [];
for i = 1:length(days_later)
    strs = [strs, string(strcat(num2str(days_later(i)), " days later"))];
end
    
data = gel_structs{1,gel_ID}.data;

for i = 1:length(data)
    temp = find(unique_dates == gel_structs{1,gel_ID}.dates_measured(i)); ind = temp(1);
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

function [fig] = plot_all_together(gel_structs, trace_num, colorscheme, hand_vis_on_indices)
fig = figure;
hold on
for i = 1:length(gel_structs)
    if isstring(trace_num) && strcmp(trace_num, "last")
        processed_data = gel_structs{1,i}.data{length(gel_structs{1,i}.data)};
    else
        processed_data = gel_structs{1,i}.data{trace_num};
    end
    
    if ismember(i, hand_vis_on_indices)
        hand_vis = 'on';
    else
        hand_vis = 'off';
    end
    
    plot(processed_data(:, 1), processed_data(:,2), 'Color', ...
        colorscheme(2+gel_structs{1,i}.batch_num,:), 'LineWidth', ...
        1.5, 'HandleVisibility', hand_vis)
end
ylabel("Force (N)")
xlabel("indentation (mm)")
title("Stiffness")

line_list = fig.Children.Children;
num_batches = length(line_list);
strs = [];
for i = 1:num_batches
    strs = [strs, strcat("Batch ", num2str(i))];
end

colors = colorscheme(3:2+num_batches, :);
leg = legend([color_legend(strs', colors)]);
leg.Box = 0;
end

function [av_ind, av_force] = plot_av(gel_structs, trace_num, disp_name)

ind = zeros(length(gel_structs{1, 1}.data{1}), 1);
force = ind;
for i = 1:length(gel_structs)
    if isstring(trace_num) && strcmp(trace_num, "last")
    processed_data = gel_structs{1,i}.data{length(gel_structs{1,i}.data)};
    else
        processed_data = gel_structs{1,i}.data{trace_num};
    end
    ind = [ind, processed_data(:, 1)];
    force = [force, processed_data(:, 2)];
end
ind = ind(:, 2:end);
force = force(:, 2:end);
av_force = mean(force, 2);
av_ind = mean(ind, 2);
std_force = std(force')';
gcf;
hold on
plot(av_ind, av_force, 'Color', ...
    'k', 'LineWidth', ...
    1.7, 'DisplayName', disp_name);

lower = av_force - std_force;
upper = av_force + std_force;
x2 = [av_ind, fliplr(av_ind')'];
inBetween = [lower, fliplr(upper')'];
inBetween(isnan(inBetween)) = 0;
h = fill(x2, inBetween, 'r',  'HandleVisibility','off');
set(h, 'facealpha', 0)
% h(1).EdgeColor = 'none';
% h(2).EdgeColor = 'none';
end

    
function [slope] = fit_run(run)
    run = run{1};
    ps = fit(run(:,1), run(:,2), 'poly1');
    slope = ps.p1;
end

function [av, high, low] = plot_av_sd_fit(gel_structs, trace_num, colorscheme_color, disp_name)

% get average for these structs
ind = zeros(length(gel_structs{1, 1}.data{1}), 1);
force = ind;
for i = 1:length(gel_structs)
    if isstring(trace_num) && strcmp(trace_num, "last")
    	processed_data = gel_structs{1,i}.data{length(gel_structs{1,i}.data)};
    else
        processed_data = gel_structs{1,i}.data{trace_num};
    end
    ind = [ind, processed_data(:, 1)];
    force = [force, processed_data(:, 2)];
end
ind = ind(:, 2:end);
force = force(:, 2:end);
av_force = mean(force, 2);
av_ind = mean(ind, 2);
std_force = std(force')';
lower = av_force - std_force;
upper = av_force + std_force;

%fit a line to average and high and low sds
av = fit(av_ind, av_force, 'poly1');
high = fit(av_ind, upper, 'poly1');
low = fit(av_ind, lower, 'poly1');

av_y = feval(av, av_ind); 
high_y = feval(high, av_ind); 
low_y = feval(low, av_ind); 

%plot
gcf;
hold on
plot(av_ind, av_y, 'Color', ...
    colorscheme_color, 'LineWidth', ...
    1.7, 'DisplayName', disp_name);

x2 = [av_ind', fliplr(av_ind')];
inBetween = [low_y', fliplr(high_y')];
inBetween(isnan(inBetween)) = 0;
h = patch(x2, inBetween, colorscheme_color,  'HandleVisibility','off');
set(h, 'facealpha', 0.2)
h(1).EdgeColor = 'none';
end

        
        