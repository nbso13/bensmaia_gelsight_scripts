function [] = processPersonData(folder_str)
%processPersonData cd's into dir and treats files in that dir as instron
%csv's, outputting a .mat file after interpolation and mean calculation.
cd(folder_str)
file_list = dir;
file_list = file_list(3:end);
person_runs = {};
%making person_runs a cell array with entries for each run
for file = 1:length(file_list)
    person_runs{file} = file_list(file).name;
end
max_ind = 2-(2400-2325)*(2/2400);
resolution_x = linspace(0, max_ind, 2325); %querry points for interpolation
%read and interpolate person runs
for i = 1:length(person_runs)
    run = parseRun(person_runs{1,i}); %parse each person run
    %interpolate each run at resolution _x
    [force, indentation] = instronInterp(run.force, run.indentation, resolution_x);
    run.force = force(1:1200);
    run.indentation = indentation(1:1200);
    run.indentor_type = 'spherical';
    run.indentor_radius = 0;
    person_runs{2,i} = run;
end

% get average and standard dev at each point
mean_trace = zeros(size(person_runs{2,1}.force));%we're trying to write the average trace
std_trace = mean_trace; %and find the std at each point
for i = 1:length(mean_trace)
    forces = zeros(length(person_runs),1); %an array at each indentation querry point of the person forces
    for j = 1:length(person_runs)
        forces(j) = person_runs{2,j}.force(i);
    end
    mean_trace(i) = nanmean(forces);
    std_trace(i) = nanstd(forces);
end
person_files = struct;
person_files.runs = person_runs;
%save mean and std data to person_files struc
mean_run = struct;
mean_run.force = mean_trace - mean_trace(1); %starts at 0
mean_run.indentation = resolution_x(1:1200);
mean_run.name = 'mean run';
person_files.mean_run = mean_run;

std_run = struct;
std_run.force = std_trace;
std_run.indentation = resolution_x(1:1200);
std_run.name = 'std trace';
person_files.std_run = std_run;

cd data
save('person_files','person_files');
disp("Human instron runs interpolated and saved as person_files.mat");
end

