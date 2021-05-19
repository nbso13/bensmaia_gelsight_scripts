function [] = plotWelch(pxx, f, color)
%plotWelch plots the normalized power spectrum as passed.
% - thanks to charles greenspon for code - 
plot(f, pxx, color);
ylabel('Power'); 
xlabel('Spatial Frequency (1/mm)');
xlim([0 7]);
yticks([]);
yticklabels({});
end

