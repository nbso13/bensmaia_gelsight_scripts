function [gain_factor, r] = gain_process(gel, no_gel, location_mat, bottom_mat, color_vec)
%gain_process takes in the gel and no gel structs containing the aligned
%profiles of the gain stimulus. Location mat gives the locations of the
%feature peaks. Bottom mat gives the swaths at which to calculate the
%bottom of the gel.
gel_vals = zeros(1,size(location_mat, 1));
no_gel_vals = zeros(1, size(location_mat, 1));

f = figure;
f.Position = [100 100 700 700];
subplot(2,2,1)
title(gel.name)
visualizeProfile(gel);
max_no_gel = max(no_gel.profile(:));
ax = gca;
ax.FontSize = 12;
ax.FontWeight = 'bold';

subplot(2,2,2)
title("Raw")
visualizeProfile(no_gel);
ax = gca;
ax.FontSize = 12;
ax.FontWeight = 'bold';
axis off;
for i = 1:size(location_mat, 1)
    bottom_gel_vals(i) = mean(sectionbyMm(gel, bottom_mat(i,:)), 'all');
    bottom_no_gel_vals(i) = mean(sectionbyMm(no_gel, bottom_mat(i,:)), 'all');
    gel_vals(i) = mean(sectionbyMm(gel, location_mat(i,:)), 'all') - bottom_gel_vals(i);
    no_gel_vals(i) = mean(sectionbyMm(no_gel, location_mat(i,:)), 'all') - bottom_no_gel_vals(i);
    
    loc_vec = location_mat(i,:);
    new_loc = [loc_vec(1), loc_vec(3), loc_vec(2)-loc_vec(1), loc_vec(4)-loc_vec(3)];
    loc_vec = bottom_mat(i,:);
    new_bottom_loc = [loc_vec(1), loc_vec(3), loc_vec(2)-loc_vec(1), loc_vec(4)-loc_vec(3)];
    
    
    subplot(2,2,1); hold on; rectangle('Position', new_loc, 'FaceColor', color_vec(i), 'EdgeColor', color_vec(i));
    rectangle('Position', new_bottom_loc, 'FaceColor', color_vec(i), 'EdgeColor', color_vec(i));
    subplot(2,2,2); hold on; rectangle('Position', new_loc, 'FaceColor', color_vec(i), 'EdgeColor', color_vec(i));
    rectangle('Position', new_bottom_loc, 'FaceColor', color_vec(i), 'EdgeColor', color_vec(i));
    
    subplot(2,2,3)
    hold on
    scatter(no_gel_vals(i), gel_vals(i), color_vec(i), 'filled');
end

x = no_gel_vals;
y = gel_vals;

subplot(2,2,3)
hold on

fit_ob = fit(x',y', 'poly1');
lin_mod = fitlm(x', y');
m = fit_ob.p1;
n = fit_ob.p2;
r = lin_mod.Rsquared.Ordinary;
ax = plot(fit_ob);
legend(ax, sprintf('y=%f*x + %f, R^2=%f',m, n, r))
xlabel("Dot Height No Gel (mm)")
ylabel("Dot Height Gel (mm)")
title("Measuring Gel Gain")
ax = gca;
ax.FontSize = 12;
ax.FontWeight = 'bold';

subplot(2,2,4);
diff = no_gel.profile - (1/m).*gel.profile;
gel.profile = diff;
min_diff = min(diff(:));
visualizeProfile(gel);
title("Scaled Difference")
ax = gca;
ax.FontSize = 12;
ax.FontWeight = 'bold';
% sgtitle(strcat(gel.name, " Gain Analysis"))

subplot(2,2,1); caxis([0 max_no_gel]);
subplot(2,2,2); caxis([0 max_no_gel]);
subplot(2,2,4); caxis([min_diff max_no_gel+min_diff]);

gain_factor = 1/m;
end

