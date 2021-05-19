function [bool] = checkSizeMatch(gel, no_gel)
%checkSizeMatch takes in two gel structs and checks to make sure size
%matches, so they can go through charactFilter
[gel_y_size, gel_x_size] = size(gel.profile);
[no_gel_y_size, no_gel_x_size] = size(no_gel.profile);


if no_gel_y_size == gel_y_size && no_gel_x_size == gel_x_size
    bool = 1;
else
    disp(strcat("Gel/NoGel y len: ", num2str(gel_y_size), " / ", num2str(no_gel_y_size)));
    disp(strcat("Gel/NoGel x len: ", num2str(gel_x_size), " / ", num2str(no_gel_x_size)));
    bool = 0;
end
end

