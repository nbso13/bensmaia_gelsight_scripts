function [run] = parseRun(filename)
%parseRun, when pointed at a csv in the directory, returns run struct. A
%wrapper for readInstronFile
[force, indentation] = readInstronFile(filename);
run = struct;
run.force = force;
run.indentation = indentation;
run.name = filename;
end

