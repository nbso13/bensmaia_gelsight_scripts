function [pxx_gel, f_gel, pxx_no_gel, f_no_gel, pxx_gel_to_nogel_ratio, ...
    pxx_gel_to_nogel_ratio_norm] = profSpectralAnalysis(gel, no_gel)
%pullResponses: given struct filenames and other hyperparams, calc firing
%rates. filenames indicate mat file name. ppm is pins per millimeter for
%touchsim model. ts amplitude indicates how much of the texture to input to
%skin mechanics for touchsim (must use entire gelsight profile). if figure
%dir is a string, save figures there.


%% Load data and process data

% save_figures = 0;
% if isstring(figure_dir)
%     dir_list = split(figure_dir, "/");
%     texture_name = dir_list(end);
%     save_figures = 1;
%     disp(strcat("Saving figures to ", figure_dir));
% end

% freq = 0.1; %2/mm
% amp = 1500;
% res = 10;
% window_size = 5000;
% gel = generate_texture("grating", freq, amp, 10, res, window_size);


% disp("Adjusting crop and sample rate...")
% if ~checkSizeMatch(gel, no_gel)
%         [gel, no_gel] = resampleToMin(gel, no_gel); %resamples to the min resolution
%     [gel, no_gel] = bruteCropFit(gel, no_gel); %crops to same size
% end


f = figure;
f.Position = [366,740,987,200];
% subplot(3,2,1)
% visualizeProfile(gel);
% title(gel.name)
% subplot(3,2,2)
% visualizeProfile(no_gel);
% title(no_gel.name)

[pxx_gel, f_gel] = welchProfile(gel);
[pxx_no_gel, f_no_gel] = welchProfile(no_gel);

%normalize

pxx_gel_norm = pxx_gel./(max(pxx_gel));
pxx_no_gel_norm = pxx_no_gel./(max(pxx_no_gel));


subplot(1, 3, 1)

plotWelch(pxx_gel, f_gel, 'c')

hold on;

plotWelch(pxx_no_gel, f_no_gel, 'r')
title(strcat(no_gel.name, " Raw PSD"));
ax = gca;
ax.FontSize = 12;
ax.FontWeight = 'bold';
% legend(["Gel", "Raw"]); 

subplot(1, 3, 2)

plotWelch(pxx_no_gel_norm, f_no_gel, 'r')
hold on;
plotWelch(pxx_gel_norm, f_gel, 'c')
title("Normalized PSD");
xlabel(""); ylabel("");
ax = gca;
ax.FontSize = 12;
ax.FontWeight = 'bold';


subplot(1, 3, 3)
interp_gel = interp1(f_gel, pxx_gel, f_no_gel);
pxx_gel_to_nogel_ratio = interp_gel./pxx_no_gel;
plot(f_no_gel, pxx_gel_to_nogel_ratio, 'k');
ylabel('Raw Ratio'); 
ax = gca;
ax.FontSize = 12;
ax.FontWeight = 'bold';

yyaxis right
interp_gel_norm = interp1(f_gel, pxx_gel_norm, f_no_gel);
pxx_gel_to_nogel_ratio_norm = interp_gel_norm./pxx_no_gel_norm;
plt = plot(f_no_gel, pxx_gel_to_nogel_ratio_norm);
% delete(plt)
ylabel('Normalized Ratio'); 
xlim([0 7]);
title("Gel : NoGel PSD Ratio");
ax = gca;
ax.FontSize = 12;
ax.FontWeight = 'bold';






% 
% subplot(2, 2, 4)
% interp_gel_norm = interp1(f_gel, pxx_gel_norm, f_no_gel);
% pxx_gel_to_nogel_ratio_norm = interp_gel_norm./pxx_no_gel_norm;
% plot(f_no_gel, pxx_gel_to_nogel_ratio_norm, 'k');
% ylabel('Ratio'); 
% xlabel('Spatial Frequency (1/mm)');
% xlim([0 10]);
% title(strcat(no_gel.name, " Normalized Gel:NoGel Ratio"));
% sgtitle(strcat(no_gel.name, " Spectral Analysis"));

% 
% subplot(6,3,7)
% visualizeProfile(new_gel);
% title("Gel after charles filtering")
% subplot(6,3,8)
% visualizeProfile(new_no_gel);
% title("No Gel after charles filtering")
% 
% [pxx_gel, f_gel] = welchProfile(new_gel);
% [pxx_no_gel, f_no_gel] = welchProfile(new_no_gel);
% 
% subplot(6, 3, 10)
% plotWelch(pxx_gel, f_gel)
% title("Gel after charles filtering")
% 
% subplot(6, 3, 11)
% plotWelch(pxx_no_gel, f_no_gel)
% title("No Gel after charles filtering")
% 
% 
% subplot(6,3,12)
% interp_gel = interp1(f_gel, pxx_gel, f_no_gel);
% plotWelch(interp_gel./pxx_no_gel, f_no_gel)
% title("Gel : No Gel ratio")
% ylabel("Ratio")
% 
% % Now evaluating disk method
% 
% 
% new_gel = removeLowFreq(gel, stopBand, "disk");
% 
% new_no_gel = removeLowFreq(no_gel, stopBand, "disk");
% 
% subplot(6,3,13)
% visualizeProfile(new_gel);
% title("Gel after disk filtering")
% subplot(6,3,14)
% visualizeProfile(new_no_gel);
% title("No Gel after disk filtering")
% 
% [pxx_gel, f_gel] = welchProfile(new_gel);
% [pxx_no_gel, f_no_gel] = welchProfile(new_no_gel);
% 
% subplot(6, 3, 16)
% plotWelch(pxx_gel, f_gel)
% title("Gel after disk filtering")
% 
% subplot(6, 3, 17)
% plotWelch(pxx_no_gel, f_no_gel)
% title("No Gel after disk filtering")
% 
% 
% subplot(6,3,18)
% interp_gel = interp1(f_gel, pxx_gel, f_no_gel);
% plotWelch(interp_gel./pxx_no_gel, f_no_gel)
% title("Gel : No Gel ratio")
% ylabel("Ratio")

end
