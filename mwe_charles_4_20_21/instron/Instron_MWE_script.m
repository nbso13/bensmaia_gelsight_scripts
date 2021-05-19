%% INSTRON min working example 4/20/21

clear
close all

addpath('data')
addpath('helper_functions')

%% Load Person Data
%run this function if you need to re-process person files.
% processPersonData('person_data');
load('person_files')

%% Visualize Data
load('toaster_tests.mat');
load('ratio_tests.mat');
load('colorscheme.mat');

f = figure;
f.Position = [100 100 700 500];
hold on
indentation = p_array_ratio{2,1}.runs{1,1}.indentation; %indentation vector as standardized
disp_names = {};
% add P10
% len_ratio_array = length(ratio_array);
% ratio_array{2,len_ratio_array+1}.mean_trace = p_array{2,4}.mean_trace(1:1200);
% ratio_array{2,len_ratio_array+1}.std_trace = p_array{2,4}.std_trace(1:1200);
% ratio_array{1, len_ratio_array+1} = 30;
for i = 1:length(ratio_array)
    mean_run = struct; %make mean_struct
    mean_run.force = ratio_array{2,i}.mean_trace;
    mean_run.indentation = indentation;
    std_run = struct;
    std_run.force = ratio_array{2,i}.std_trace;
    std_run.indentation = indentation;
    [x2, inBetween] = stdevSpan(mean_run, std_run); %find the fill parameters for std
    h = fill(x2, inBetween, colorscheme(i+3,:),  'HandleVisibility','off');
    h.EdgeColor = 'none';
    set(h,'facealpha',.4)
    plot(indentation, mean_run.force, 'Color', colorscheme(i+3,:), 'LineWidth', ...
        1.5);
    disp_names{i, 1} = char(strcat(num2str(round(100/(ratio_array{1,i}+1), 2)), "%"));
end

%plot 1:30 from best toaster test
% mean_run = struct; %make mean_struct
% mean_run.force = p_array{2,4}.mean_trace(1:1200); %this is P10, 1:30 baked for 3 hrs
% mean_run.indentation = indentation;
% std_run = struct;
% std_run.force = p_array{2,4}.std_trace(1:1200); %this is P10, 1:30 baked for 3 hrs
% std_run.indentation = indentation;

% [x2, inBetween] = stdevSpan(mean_run, std_run); %find the fill parameters for std
% h = fill(x2, inBetween, colorscheme(i+4,:),  'HandleVisibility','off');
% h.EdgeColor = 'none';
% set(h,'facealpha',.4)
% plot(indentation, mean_run.force, 'Color', colorscheme(i+4,:), 'LineWidth', 1.5, 'DisplayName', '3.23%');
% disp_names{i+1, 1}  = '3.23%';

%plot person
[x2, inBetween] = stdevSpan(person_files.mean_run, person_files.std_run);
h = fill(x2, inBetween, 'r',  'HandleVisibility','off');
set(h,'facealpha',.4)
h.EdgeColor = 'none';
plot(person_files.mean_run.indentation, person_files.mean_run.force, 'Color', 'k', 'LineWidth', 2.2, 'DisplayName', 'Human');
disp_names{i+1, 1}  = 'Human';
leg_text = color_legend(disp_names, [colorscheme(4:length(ratio_array)+3, :); [0 0 0]]);
leg = legend(leg_text);
leg.Box = 0;
ax = gca;
ax.FontSize = 12;
ax.FontWeight = 'bold';

xlabel("Indentation Depth (mm)")
ylabel("Force (N)");

%% COMPARING MEAN RESIDUAL to HUMAN for each gel

resids = zeros(2,length(ratio_array));
ratios = zeros(1,length(ratio_array));
% gradient_diffs = resids;
human_run = person_files.mean_run.force;
% indentation_grad = gradient(person_files.mean_run.indentation(:));
% human_gradients = gradient(human_run(:)) ./ indentation_grad;
for i = 1:length(ratio_array)
    residuals = human_run./ratio_array{2,i}.mean_trace;
%     norm_residuals = residuals./ human_run;
%     norm_residuals(1) = 0;
    mean_norm_resid = mean(residuals);
    sd_norm_resid = std(residuals);
    resids(1, i) = mean_norm_resid;
    resids(2, i) = sd_norm_resid;
    ratios(i) = ratio_array{1,i};
%     gradient_resid = abs(human_gradients - ...
%         (gradient(ratio_array{2,i}.mean_trace(:)) ./ indentation_grad));
%     gradient_diffs(1, i) = mean(gradient_resid);
%     gradient_diffs(2, i) = std(gradient_resid);
end

ratios = 100./(ratios + 1);

figure;
hold on
bar(ratios, resids(1,:));
er = errorbar(ratios, resids(1,:), resids(2,:));
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
xlabel("Percentage Activator")
ylabel("Mean Ratio to Human Trace");
% title("Mean Ratio of Gel to Human Trace as a Function of Activator Content")
yline(1)
ax = gca;
ax.FontSize = 12;
ax.FontWeight = 'bold';