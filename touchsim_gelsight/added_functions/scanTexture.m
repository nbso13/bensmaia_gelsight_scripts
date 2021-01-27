function [trace] = scanTexture(offset, trace, x_len, pinspermm, speed, samp_freq)
%scanTexture takes in a vector of pinoffsets (the
%same size as the x length of the stimulus trace). It then takes the offset
%and rewrites the trace matrix by wrapping it around forward for each
%timepoint for each  TIME RESOLUTION SHOULD BE MUCH HIGHER THAN SPATIAL otherwise
%it's just one pin per time increment. 
offset_len = length(offset);
pinspersec = pinspermm*speed; %pins per second (speed)
timeperpin = ceil(samp_freq/pinspersec); % (time_intervals per pin move)
if ne(size(trace,2), length(offset))
    error("trace width and offset length are unequal and do not match to the same pins");
end
trace(1,:) = offset;
speed_counter = 1; %keeps track of how many time intervals since last move
num_x_swaths = offset_len/x_len;
for i = 2:size(trace,1)
    if speed_counter < timeperpin %if its not yet time to move up
        speed_counter = speed_counter+1; %increment counter
        trace(i,:) = trace(i-1,:); %trace does not move
    
    else %rotate within each x swath
        for j= 1:num_x_swaths %for every x swath
            ind_1 = 1+(j-1)*x_len;
            ind_2 = ind_1+x_len-1;
            x_swath = offset(ind_1:ind_2);
            x_swath_rotated = [x_swath(end), x_swath(1:end-1)];
            offset(ind_1:ind_2) = x_swath_rotated;
        end
        trace(i, :) = offset; 
        speed_counter = 1; %reset counter
    end
end
end

