%% New Instron Processing
% Read and Process files from ratio test in july 27th

close all
clear
% addpath('helper_functions');
% cd ../simulating_responses
addpath('helper_functions');
load('colorscheme.mat');

batch_edges = [1, 5, 10, 14, 18, 22, 25, 28, 32, 36, 39, 42, 45, 48, 51, 54, 57, 60];
disp_names = ["batch 1 - 2.8 - ", "batch 2 - 2.8 -", "batch 3 - 2.6 -", ...
    "batch 4 - 2.7 -", "batch 5 - 2.9 -", "batch 6 - 3.1 - ", ...
    "batch 7 - 3.3 - ", "batch 8 - 3.8 - ", "batch 9 - 4.2 - ", ...
    "nick", "drew", "ev", "olivia", "paul", "char", "liza"];
batch_dates = [210709, 210709, 210719, 210719, 210722];
fit_type = "poly2";
%% B1 and B2

%data structure - cell array of structs. Each struct is one gel. Each
%struct has a field "traces" with a cell array of trace recordings. Each cell has a nx2 array
%(mm and forces). For each entry in that cell array there's an entry in the
%array in the 'dates' field of the 
gel_structs = {};
date_made_B1_B2 = 210709;
date_measured = 210727;

cd 'data'
filename = '210726_redo_batches_1.xls';
ssds = spreadsheetDatastore(filename);
measurements = sheetnames(ssds, 1);
cd ..


for i = 1:length(measurements)-2
    if i > 4
        x = 2;
    else
        x = 1;
    end
    gel_structs{1,i} = struct;
    gel_id_str = measurements(i+2);
    gel_structs{1,i}.batch_num = x;
    gel_structs{1,i}.id = i;
    gel_structs{2,i} = gel_id_str;
    gel_structs{1,i}.date_made = date_made_B1_B2;
    gel_structs{1,i}.dates_measured = [];
    gel_structs{1,i}.dates_measured = [gel_structs{1,i}.dates_measured, date_measured];
    cd 'data'
    gel_structs{1,i}.data = {};
    temp = readcell(filename, 'Sheet', gel_id_str, 'Range', 'A:B');
    gel_structs{1,i}.data{length(gel_structs{1,i}.data)+1} = cell2mat(temp(5:end, :));
    cd ..
    
end
% examine data
trace_number = 1;
[fig] = plot_all_together(gel_structs, trace_number, colorscheme, batch_edges, disp_names);

% determine contact

trace_number = 1;
for i = 1:length(measurements)-2
    processed_data = contact_process(i, trace_number, gel_structs);
    gel_structs{1,i}.data{trace_number} = processed_data;
end

[fig] = plot_all_together(gel_structs, trace_number, colorscheme, batch_edges, disp_names);
title("Batches 1 and 2, 210727")

figure;
batch_nums = [1,2];
plot_from_batch_nums(gel_structs, trace_number, colorscheme, batch_nums, batch_edges, disp_names)

%% B3 B4
date_made_B3_B4 = 210719;
date_measured = 210727;

cd 'data'
filename = '210726_redo_batches_2.xls';
ssds = spreadsheetDatastore(filename);
measurements = sheetnames(ssds, 1);
cd ..
current_gel_list_length = length(gel_structs);

for j = 1:length(measurements)-2
    if j > 4
        x = 4;
    else
        x = 3;
    end
    i = j+current_gel_list_length;
    gel_structs{1,i} = struct;
    gel_id_str = measurements(j+2);
    gel_structs{1,i}.batch_num = x;
    gel_structs{1,i}.id = i;
    gel_structs{2,i} = gel_id_str;
    gel_structs{1,i}.date_made = date_made_B3_B4;
    gel_structs{1,i}.dates_measured = [];
    gel_structs{1,i}.dates_measured = [gel_structs{1,i}.dates_measured, date_measured];
    cd 'data'
    gel_structs{1,i}.data = {};
    temp = readcell(filename, 'Sheet', gel_id_str, 'Range', 'A:B');
    gel_structs{1,i}.data{length(gel_structs{1,i}.data)+1} = cell2mat(temp(5:end, :));
    cd ..
    
end



% examine data
trace_number = 1;
[fig] = plot_all_together(gel_structs, trace_number, colorscheme, batch_edges, disp_names);

% determine contact

trace_number = 1;
for i = current_gel_list_length+1:current_gel_list_length+length(measurements)-2
    processed_data = contact_process(i, trace_number, gel_structs);
    gel_structs{1,i}.data{trace_number} = processed_data;
end

[fig] = plot_all_together(gel_structs, trace_number, colorscheme, batch_edges, disp_names);
title("Batches 1, 2, 3, 4, 210727")


figure;
batch_nums = [3,4];
plot_from_batch_nums(gel_structs, trace_number, colorscheme, batch_nums, batch_edges, disp_names)

%% B5
date_made_B5 = 210722;
date_measured = 210727;

cd 'data'
filename = '210727_B5_1.xls';
ssds = spreadsheetDatastore(filename);
measurements = sheetnames(ssds, 1);
cd ..
current_gel_list_length = length(gel_structs);

for j = 1:length(measurements)-2
    x=5;
    i = j+current_gel_list_length;
    gel_structs{1,i} = struct;
    gel_id_str = measurements(j+2);
    gel_structs{1,i}.batch_num = x;
    gel_structs{1,i}.id = i;
    gel_structs{2,i} = gel_id_str;
    gel_structs{1,i}.date_made = date_made_B5;
    gel_structs{1,i}.dates_measured = [];
    gel_structs{1,i}.dates_measured = [gel_structs{1,i}.dates_measured, date_measured];
    cd 'data'
    gel_structs{1,i}.data = {};
    temp = readcell(filename, 'Sheet', gel_id_str, 'Range', 'A:B');
    gel_structs{1,i}.data{length(gel_structs{1,i}.data)+1} = cell2mat(temp(5:end, :));
    cd ..
    
end


% examine data
trace_number = 1;
[fig] = plot_all_together(gel_structs, trace_number, colorscheme, batch_edges, disp_names);

% determine contact

trace_number = 1;
for i = current_gel_list_length+1:current_gel_list_length+length(measurements)-2
    processed_data = contact_process(i, trace_number, gel_structs);
    gel_structs{1,i}.data{trace_number} = processed_data;
end

[fig] = plot_all_together(gel_structs, trace_number, colorscheme, batch_edges, disp_names);
title("Batches 1, 2, 3, 4, 5 210727")


figure;
batch_nums = [1, 2, 5];
plot_from_batch_nums(gel_structs, trace_number, colorscheme, batch_nums, batch_edges, disp_names)


%% B6, 7
date_made_B6B7 = 210726;
date_measured = 210728;

cd 'data'
filename = '210728_B6_B7_new-indentor.xls';
ssds = spreadsheetDatastore(filename);
measurements = sheetnames(ssds, 1);
cd ..
current_gel_list_length = length(gel_structs);

for j = 1:length(measurements)-2
    if j > 3
        x = 7;
    else
        x= 6;
    end
    i = j+current_gel_list_length;
    gel_structs{1,i} = struct;
    gel_id_str = measurements(j+2);
    gel_structs{1,i}.batch_num = x;
    gel_structs{1,i}.id = i;
    gel_structs{2,i} = gel_id_str;
    gel_structs{1,i}.date_made = date_made_B6B7;
    gel_structs{1,i}.dates_measured = [];
    gel_structs{1,i}.dates_measured = [gel_structs{1,i}.dates_measured, date_measured];
    cd 'data'
    gel_structs{1,i}.data = {};
    temp = readcell(filename, 'Sheet', gel_id_str, 'Range', 'A:B');
    gel_structs{1,i}.data{length(gel_structs{1,i}.data)+1} = cell2mat(temp(5:end, :));
    cd ..  
end


% examine data
trace_number = 1;
[fig] = plot_all_together(gel_structs, trace_number, colorscheme, batch_edges, disp_names);

% determine contact

trace_number = 1;
for i = current_gel_list_length+1:current_gel_list_length+length(measurements)-2
    processed_data = contact_process(i, trace_number, gel_structs);
    gel_structs{1,i}.data{trace_number} = processed_data;
end

[fig] = plot_all_together(gel_structs, trace_number, colorscheme, batch_edges, disp_names);
title("Batches 1, 2, 3, 4, 5 210727")


figure;
batch_nums = [6, 7];
plot_from_batch_nums(gel_structs, trace_number, colorscheme, batch_nums, batch_edges, disp_names)


%% B8, 9
date_made_B8B9 = 210728;
date_measured = 210729;

cd 'data'
filename = '210729_B8B9.xls';
ssds = spreadsheetDatastore(filename);
measurements = sheetnames(ssds, 1);
cd ..
current_gel_list_length = length(gel_structs);

for j = 1:length(measurements)-2
    if j > 2
        x = 9;
    else
        x= 8;
    end
    i = j+current_gel_list_length;
    gel_structs{1,i} = struct;
    gel_id_str = measurements(j+2);
    gel_structs{1,i}.batch_num = x;
    gel_structs{1,i}.id = i;
    gel_structs{2,i} = gel_id_str;
    gel_structs{1,i}.date_made = date_made_B8B9;
    gel_structs{1,i}.dates_measured = [];
    gel_structs{1,i}.dates_measured = [gel_structs{1,i}.dates_measured, date_measured];
    cd 'data'
    gel_structs{1,i}.data = {};
    temp = readcell(filename, 'Sheet', gel_id_str, 'Range', 'A:B');
    gel_structs{1,i}.data{length(gel_structs{1,i}.data)+1} = cell2mat(temp(5:end, :));
    cd ..
end


% examine data
trace_number = 1;
[fig] = plot_all_together(gel_structs, trace_number, colorscheme, batch_edges, disp_names);

% determine contact

trace_number = 1;
for i = current_gel_list_length+1:current_gel_list_length+length(measurements)-2
    processed_data = contact_process(i, trace_number, gel_structs);
    gel_structs{1,i}.data{trace_number} = processed_data;
end

chosen_structs = gel_structs(:, batch_edges(8):batch_edges(10)-1);
[fig] = plot_all_together(chosen_structs, trace_number, colorscheme, batch_edges, disp_names);
title("Batches 8 and 9")


figure;
batch_nums = [8, 9];
plot_from_batch_nums(gel_structs, trace_number, colorscheme, batch_nums, batch_edges, disp_names)



%% fingers_1
close all
date_measured = 210726;

cd 'data'
filename = '210726_fingers.xls';
ssds = spreadsheetDatastore(filename);
measurements = sheetnames(ssds, 1);
cd ..
current_gel_list_length = length(gel_structs);

for j = 1:length(measurements)-2
    if j < 4
        x = "nick";
    elseif j < 7
        x = "drew";
    elseif j < 10
        x = "ev";
    elseif j < 13
        x = "olivia";
    else
        error("index not recognized")
    end
    i = j+current_gel_list_length;
    gel_structs{1,i} = struct;
    gel_id_str = measurements(j+2);
    gel_structs{1,i}.batch_num = x;
    gel_structs{1,i}.id = i;
    gel_structs{2,i} = gel_id_str;
    gel_structs{1,i}.dates_measured = [];
    gel_structs{1,i}.dates_measured = [gel_structs{1,i}.dates_measured, date_measured];
    cd 'data'
    gel_structs{1,i}.data = {};
    if ~ (x == "nick")
        temp = readcell(filename, 'Sheet', ...
            strcat(x, "_right_index_", num2str(mod(j,3)+ ~logical(mod(j,3))*3)), ...
            'Range', 'A:B');
        gel_structs{1,i}.data{length(gel_structs{1,i}.data)+1} = cell2mat(temp(5:end, :));
    else
        temp = readcell(filename, 'Sheet', ...
            strcat(x, "_left_index_", num2str(mod(j,3)+ ~logical(mod(j,3))*3)), ...
            'Range', 'A:B');
        gel_structs{1,i}.data{length(gel_structs{1,i}.data)+1} = cell2mat(temp(5:end, :));
    end
    cd ..
    
end


% examine data
trace_number = 1;
[fig] = plot_all_together(gel_structs, trace_number, colorscheme, batch_edges, disp_names);

% determine contact

trace_number = 1;
for i = current_gel_list_length+1:current_gel_list_length+length(measurements)-2
    processed_data = contact_process(i, trace_number, gel_structs);
    gel_structs{1,i}.data{trace_number} = processed_data;
end

[fig] = plot_all_together(gel_structs, trace_number, colorscheme, batch_edges, disp_names);
title("Batches 1, 2, 3, 4, 5 210727")




%% fingers_2
date_measured = 210727;

cd 'data'
filename = '210727_more_fingers_2.xls';
ssds = spreadsheetDatastore(filename);
measurements = sheetnames(ssds, 1);
cd ..
current_gel_list_length = length(gel_structs);

for j = 1:length(measurements)-2
    if j < 4
        x = "paul";
    elseif j < 7
        x = "char";
    elseif j < 10
        x = "liza";
    else
        error("index not recognized")
    end
    i = j+current_gel_list_length;
    gel_structs{1,i} = struct;
    gel_id_str = measurements(j+2);
    gel_structs{1,i}.batch_num = x;
    gel_structs{1,i}.id = i;
    gel_structs{2,i} = gel_id_str;
    gel_structs{1,i}.dates_measured = [];
    gel_structs{1,i}.dates_measured = [gel_structs{1,i}.dates_measured, date_measured];
    cd 'data'
    gel_structs{1,i}.data = {};
    temp = readcell(filename, 'Sheet', ...
            strcat(x, "_right_index_", num2str(mod(j,3)+ ~logical(mod(j,3))*3)), ...
            'Range', 'A:B');
    gel_structs{1,i}.data{length(gel_structs{1,i}.data)+1} = cell2mat(temp(5:end, :));
    cd ..
    
end


% examine data
trace_number = 1;
[fig] = plot_all_together(gel_structs, trace_number, colorscheme, batch_edges, disp_names);

% determine contact

trace_number = 1;
for i = current_gel_list_length+1:current_gel_list_length+length(measurements)-2
    processed_data = contact_process(i, trace_number, gel_structs);
    gel_structs{1,i}.data{trace_number} = processed_data;
end


% all fingers
figure;
subplot(2, 2, 1)
ylim([-0.05, 0.4])
batch_nums = [10:16];
finger_slopes = plot_from_batch_nums(gel_structs, trace_number, colorscheme, batch_nums, batch_edges, disp_names);
title("all fingers")

% all gels
subplot(2, 2, 2);
ylim([-0.05, 0.4])
batch_nums = [1, 5, 6, 7, 8, 9];
plot_from_batch_nums(gel_structs, trace_number, colorscheme, batch_nums, batch_edges, disp_names)
title("all gels")


% plot human finger average
subplot(2, 2, 3)
ylim([-0.05, 0.4])
target_structs = gel_structs(:,28:48);
mean_finger_slope = plot_av_sd_fit(target_structs, trace_number, colorscheme(13, :), "finger average");
title("Average Finger")

subplot(2, 2, 4)
ylim([-0.05, 0.4])
target_structs = gel_structs(:,horzcat(28:33, 40:46));
plot_av_sd_fit(target_structs, trace_number, colorscheme(19, :), "male average");
target_structs = gel_structs(:,horzcat(34:39, 46:48));
plot_av_sd_fit(target_structs, trace_number, colorscheme(17, :), "female average");
title("Male and female averages")


% get different stiffness slopes for each ratio
figure;
batch_nums = [1, 5, 6, 7, 8, 9];
slopes = plot_from_batch_nums(gel_structs, trace_number, colorscheme, batch_nums, batch_edges, disp_names);

%plot stiffness slopes vs gel ratio
ratios = [2.8, 2.9, 3.1, 3.3, 3.8, 4.2];
ratio_strs = string(ratios);
figure; scatter(ratios, slopes);
hold on
xlabel("ratios"); ylabel("slope");
best_fit = fit(ratios', slopes', 'poly1');
ys = feval(best_fit, ratios);
plot(ratios, ys)
x = 5;
% beeswarm finger slopes
bee_swarm(finger_slopes, 5, 'r', length(finger_slopes));
title("Finger Stiffness Slopes Distribution")
ylabel("Stiffness, N/mm")

xticks([ratios, 5]);
xticklabels([ratio_strs, "finger slopes"]);



figure;
batch_nums = [10:16];
finger_slopes = plot_from_batch_nums(gel_structs, trace_number, colorscheme, batch_nums, batch_edges, disp_names);
target_structs = gel_structs(:,horzcat(34:39, 46:48));
plot_av_sd_fit(target_structs, trace_number, colorscheme(17, :), "female average");
title("Gels and Female average")

%% functions

function [processed_data] = contact_process(gel_ID, trace_num, gel_structs)
if isstring(trace_num) && strcmp(trace_num, "last")
    temp_data = gel_structs{1,gel_ID}.data{length(gel_structs{1,gel_ID}.data)};
else
    temp_data = gel_structs{1,gel_ID}.data{trace_num};
end
initial_indentation = 0.08; %mm
mult = 4.1; %sd's to exceed
resolution_x = linspace(0, 1.5, 1500);
% get average force in first n = 0.05 mm
% get data from before in that initial indentation value
initial = temp_data((temp_data(:,1) < initial_indentation), :);
temp_force = temp_data(:, 2);
mean_force = mean(initial(:, 2)); % get the mean force in that initial descent
sd_force = std(initial(:,2));
% find first index where indentation is greater than init ind value
more_than_indent = find(temp_data(:,1) > initial_indentation, 1);
% find the first index, beyond the init indent value, where force exceeds
% mean force + mult * sd
ind_contact = find(temp_force(more_than_indent:end) > (mean_force + mult*sd_force), 1);
% add index for first indent to be able to index original data
ind_contact = ind_contact + more_than_indent;
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

function [fig] = plot_all_together(gel_structs, trace_num, colorscheme, batch_edges, disp_names)
fig = figure;
hold on
for i = 1:length(gel_structs)
    if isstring(trace_num) && strcmp(trace_num, "last")
        processed_data = gel_structs{1,i}.data{length(gel_structs{1,i}.data)};
    else
        processed_data = gel_structs{1,i}.data{trace_num};
    end
    
    if ismember(i, batch_edges)
        hand_vis = 'on';
    else
        hand_vis = 'off';
    end
    
    if isstring(gel_structs{1,i}.batch_num)
        color_ind = find(disp_names == gel_structs{1,i}.batch_num);
        color_ind = color_ind(1);
    else
        color_ind = gel_structs{1,i}.batch_num;
    end
    
    plot(processed_data(:, 1), processed_data(:,2), 'Color', ...
        colorscheme(2+color_ind,:), 'LineWidth', ...
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

x2 = [av_ind', fliplr(av_force')];
inBetween = [low_y', fliplr(high_y')];
inBetween(isnan(inBetween)) = 0;
h = patch(x2, inBetween, colorscheme_color,  'HandleVisibility','off');
set(h, 'facealpha', 0.2)
h(1).EdgeColor = 'none';
% h(1).EdgeColor = 'none';
% h(2).EdgeColor = 'none';
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

function [slopes] = plot_from_batch_nums(gel_structs, trace_num, colorscheme, batch_nums, batch_edges, disp_names)
%wrapper for plot_av_sd_fit
slopes = [];
for i = batch_nums
    target_structs = gel_structs(:,batch_edges(i):batch_edges(i+1)-1);
    slope = plot_av_sd_fit(target_structs, trace_num, colorscheme(i+2, :), disp_names(i));
    slopes = [slopes, slope.p1];
end


colors = colorscheme(2+batch_nums, :);
leg = legend([color_legend(disp_names(batch_nums)', colors)]);
leg.Box = 0; 

end
        