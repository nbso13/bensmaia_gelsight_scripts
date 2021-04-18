function [] = plotWelch(pxx, f)
%plotWelch plots the normalized power spectrum as passed.
% - thanks to charles greenspon for code - 
plot(f, pxx, 'k');
ylabel('Power'); 
xlabel('Spatial Frequency (mm)');
xlim([0.25 10]);
yticks([]);
yticklabels({});
end

