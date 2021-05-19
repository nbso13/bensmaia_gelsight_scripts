function leg_text = color_legend(input_text, colors)
    % input text = cell(n,1) with chars in each cell
    % colors = double(n,3) in RGB format with range 0-1
    
    if size(input_text,1) ~= size(colors,1)
        display('Error: uneven number of text strings and colors.')
        return
    end

    leg_text = [];
    for i = 1:size(input_text,1)
        leg_text = [leg_text, {['\color[rgb]{', num2str(colors(i,:)), '}',  input_text{i,1}]}];
    end

end