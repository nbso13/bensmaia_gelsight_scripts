function [force, indentation] = readInstronFile(filename)
%readInstronFile parses the csv with 'filename' and returns two arrays
%of force and indentation traces
starting_row = 9;
final_row = 2336;
temp = csvread(filename, starting_row, 1, [starting_row 1 final_row 2]);
force = temp(:, 2);
indentation = temp(:,1);
end

