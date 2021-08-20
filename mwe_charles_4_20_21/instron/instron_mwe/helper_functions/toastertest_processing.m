%% Process Instron Data
% Read and process instron data from toaster tests in december
close all
clear
addpath('../')        
cd('../data_raw/Toaster Tests')

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
p_array = {};
max_ind = 2-(2400-2325)*(2/2400);
resolution_x = linspace(0, max_ind, 2325); %querry points for interpolation

for i = 1:length(folder_list)
    name = folder_list(i).name; %get prototype name (P10)
    prototype = struct;
    cd(name) %into prototype trial folder
    prototype.name = name; %name field in struct
    prototype.ratio = 30; %ratio silocon gel parts A:1 B
    prototype.gel_height = 5; %gel height in mm
    prototype.trials = {}; %cell array of trial structs
    trial_list = dir; %get list of trial directories
    trial_list = trial_list(3:end);
    
    for j = 1:length(trial_list) %for every trial
        trial_dir = trial_list(j).name; %get trial directory
        trial = struct;
        trial.prototype = name;
        trial.time = getTime(trial_dir); %pull trial time from name
        trial.indentor_radius = 0; %indentor radius in mm
        trial.indentor_type = 'sphere'; 
        trial.runs = {}; % cell array of runs
        cd(trial_dir) %pop into trial directory
        run_list = dir;
        run_list = run_list(3:end);
        
        for k = 1:length(run_list) %for every run
            run = struct;
            run.num = k;
            [force, indentation] = readInstronFile(run_list(k).name); %read csv
            [force, indentation] = instronInterp(force, indentation, resolution_x);
            run.force = force;
            run.prototype = name;
            run.indentation = indentation;
            run.time = trial.time;
            run.indentor_radius = trial.indentor_radius;
            run.indentor_type = trial.indentor_type;
            trial.runs{k} = run; %store run struct in runs cell array for this trial
        end
        cd ..
        prototype.trials{2, j} = trial; %store trial struct
        prototype.trials{1, j} = trial.time;
    end
    cd ..
    p_array{2,i} = prototype; %store prototype struct
    p_array{1,i} = prototype.name;
end

%% Compile All Runs In Structs
complete_runs = {}; % cell array that stores all the runs
proto_len = length(p_array);
run_counter = 0;
for i = 1:proto_len % for each prototype
    all_runs = {};% a cell array that stores all runs per prototype
    forces = zeros(1,length(run.force));
    trial_num = length(p_array{2,i}.trials);
    for j = 1:trial_num % for each trial
        run_num = length(p_array{2,i}.trials{2,j}.runs);
        for k = 1:run_num % go thru runs
            ind = run_num*(j-1)+k; % all run index
            run = p_array{2,i}.trials{2,j}.runs{k}; %current run
            forces = [forces; run.force];
            all_runs{1, ind} = run.time;% trial time
            all_runs{2, ind} = run.num;% run num
            all_runs{3,ind} = run; %copy run 
            run_counter = run_counter+1;
            complete_runs{1, run_counter} = run.prototype; %write run prototype
            complete_runs{2, run_counter} = all_runs{1, ind}; %copy from all runs
            complete_runs{3, run_counter} = all_runs{2, ind}; %copy from all runs
            complete_runs{4, run_counter} = all_runs{3, ind}; %copy from all runs
        end
    end
    p_array{2,i}.all_runs = all_runs; %save all runs as struct field
    forces = forces(2:end, :); %kill that first row of zeros
    %now this is a matrix where rows are runs and cols are forces at a
    %given aligned indentation along the shared x axis (in mm)
    mean_trace = nanmean(forces);
    std_trace = nanstd(forces);
    p_array{2,i}.mean_trace = mean_trace;
    p_array{2,i}.std_trace = std_trace; 
end
p_array{1, end+1} = 'complete runs';
p_array{2, end} = complete_runs;

cd ../data
save('toaster_tests', 'p_array')
disp("p_array saved as toaster_tests.mat")    