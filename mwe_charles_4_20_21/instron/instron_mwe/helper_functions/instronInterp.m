function [force, indentation] = instronInterp(force, indentation, res_x)
%instronInterp takes in force and indentation values and interpolates
%evenly and recalculates those values, returning two vectors
temp_ax = indentation(2:end);
temp_f = force(2:end);
[~,IA,~] = unique(temp_ax);
temp_ax = temp_ax(IA);
temp_f = temp_f(IA);
force = interp1(temp_ax,temp_f, res_x);
indentation = res_x;
end

