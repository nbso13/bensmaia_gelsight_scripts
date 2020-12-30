function [] = gainAnalysis(gel, no_gel, resolution, predicted_ratio)
%gainAnalysis takes in a gel and no gel of the same image and compares the
%ratio of gel to no gel in terms of height at a given resolution
ratio = zeros(floor((1/resolution)*gel.x_axis(end)), floor((1/resolution)*gel.y_axis(end)));
no_gello = zeros(floor((1/resolution)*gel.x_axis(end)), floor((1/resolution)*gel.y_axis(end)));
gello = no_gello;
for i=resolution:resolution:gel.x_axis(end)-resolution
    for j = resolution:resolution:gel.y_axis(end)-resolution
        no_gel_cross = sectionbyMm(no_gel, [i, i+resolution, j, j+resolution]);
        gel_cross = sectionbyMm(gel, [i, i+resolution, j, j+resolution]);
        no_gel_cross_mean = mean(no_gel_cross, 'all');
        gel_cross_mean = mean(gel_cross, 'all');
        no_gello(floor((1/resolution)*i),floor((1/resolution)*j)) = no_gel_cross_mean;
        gello(floor((1/resolution)*i),floor((1/resolution)*j)) = gel_cross_mean;
        ratio(floor((1/resolution)*i),floor((1/resolution)*j)) = gel_cross_mean/no_gel_cross_mean;
    end
end

visualizeProfile(no_gel);
title("Gain Stim No Gel")

visualizeProfile(gel);
title("Gain Stim with Gel")

figure
x = no_gello(:);
y = gello(:);
scatter(x, y, 1)
title("Gel to No Gel Ratio")
xlabel("height, no gel, mm")
ylabel("height, gel, mm")
xlim([0 inf])
ylim([0 inf])
hold on
x= 0:0.2:1.2;
y = x*predicted_ratio;
ax = plot(x,y, 'r');
legend(ax, sprintf('y=%f*x', predicted_ratio))

figure
scatter(no_gello(:), ratio(:), 1)
title("Feature Height by Gel/No Gel Ratio")
xlabel("height, no gel, mm")
ylabel("ratio, gel : no gel")
xlim([0 inf])
ylim([0 5])
hold on
x= 0:0.2:1.2;
y = ones(1,length(x)) * predicted_ratio;
ax = plot(x,y, 'r');
legend(ax, sprintf('y=%f*x', predicted_ratio))

diff = gel;
diff.profile = no_gel.profile - gel.profile./predicted_ratio;
visualizeProfile(diff)


end

