%% Ratio Test Processing
% Read and Process files from ratio test in january

close all
clear
addpath('../')        
cd('../data_raw/Testing New Molds_1_28\all tests\')

%% Read csvs, store as structs
% prototype struct fields:
%   name, ratio, gel_height (mm), trials (list of trial structs),
%   all_runs (list of run_structs)
% trial struct fields:
%    time, indentor_radius (in mm), indentor_type 
% run struct fields:
%   time, indentor_radius (in mm), indentor_type (1 is spherical, 2 
% is cylindrical), run_num, force (array), indentation (array)

folder_list = dir;
folder_list = folder_list(3:end);
p_array_ratio = {};
max_ind = 2-(2400-2325)*(2/2400);
resolution_x = linspace(0, max_ind, 2325); %querry points for interpolation

for i = 1:length(folder_list)
     name = folder_list(i).name; %get prototype name (P10)
     prototype = struct;
     cd(name) %into prototype folder
     name = strsplit(name, '.');
     prototype.name = name{1}; %name of gel is ratio_IDnumber
     ratio = strsplit(prototype.name, '_');
     prototype.ratio = str2num(ratio{1}); %ratio silocon gel parts A:1 B
     prototype.gel_height = 5; %gel height in mm
     prototype.runs = {}; %cell array of trial structs
     run_list = dir; %get list of trial directories
     run_list = run_list(3:end);
     all_forces = zeros(1,length(resolution_x)); % a big matrix of all the forces
     for j = 1:length(run_list)
         run = struct;
         run.gel_loc = nan;
         if j < 4
             run.gel_loc = "center";
         elseif j<7
             run.gel_loc = "middle";
         elseif j <10
             %run.gel_loc = "edge";
             continue %THROWING OUT EDGE READINGS
         end
         [force, indentation] = readInstronFile(run_list(j).name); %read csv
         [force, indentation] = instronInterp(force, indentation, resolution_x);
         
         %limit to linear part at less than one..
         less_than_one_logit = (indentation<=1);
         clipped_ind = indentation(less_than_one_logit);
         clipped_force = force(less_than_one_logit);
         
         all_forces = [all_forces; force];
         run.force = clipped_force;
         run.indentation = clipped_ind;
         prototype.runs{j} = run;
     end
     all_forces = all_forces(2:end, :); %kill that first row of zeros
     all_forces = all_forces(:, less_than_one_logit); % kill measured forces more than one
     prototype.all_forces = all_forces;
     cd ..
     p_array_ratio{1,i} = prototype.name;
     p_array_ratio{2,i} = prototype;
end

%% Tabulate Together Runs For Each Ratio, calculate averages and stds
% write a cell array for each ratio with first row, ratio, second row,
% matrix of rows of forces from each run for that ratio
ratio_array = {};
for i = 1:length(p_array_ratio) % for every prototype
    forces = p_array_ratio{2,i}.all_forces;
    ratio = p_array_ratio{2,i}.ratio;
    match = 0;
    for j = 1:size(ratio_array,2) %check every known ratio
        if ratio == ratio_array{1,j} %if it matches this one
            ratio_array{2,j} = [ratio_array{2,j}; forces]; % add all the forces
            match = 1; %set match flag
            continue %skip the rest
        end
    end
    if not(match) %if no match found
        ind = size(ratio_array,2)+1;
        ratio_array{1, ind} = ratio; %label first row with ratio
        ratio_array{2, ind} = forces;
    end
end

% make mean_run and std_run structs
for i = 1:length(ratio_array)
    prototype = struct;
    prototype.mean_trace = nanmean(ratio_array{2,i});
    prototype.std_trace = nanstd(ratio_array{2,i});
    slope = regress(prototype.mean_trace',clipped_ind');
    prototype.slope =slope;
    ratio_array{2,i} = prototype;
end

cd ('../../../data')
save('ratio_tests', 'p_array_ratio', 'ratio_array')
disp("p_array_ratio and ratio_array saved as ratio_tests.mat")    
        
        
        
        