names = ["wool blend", "velvet", "hucktowel", "sueded cuddle", "corduroy", "blizzard fleece"];
close all
figure;

subplot(1,2,1)
hold on
for i = 1:size(forces_depths_GS, 2)
    scatter(forces_depths_TS(1,i), forces_depths_TS(2,i), 'filled')
end
xlabel("indentation max depth, mm")
ylabel("total forces")
legend( names)
title("TouchSim")
ylim([0, 0.6])
xlim([0 1.3])

subplot(1,2,2)
hold on
for i = 1:size(forces_depths_GS, 2)
    scatter(forces_depths_GS(1,i), forces_depths_GS(2,i), 'filled')
end

xlabel("indentation max depth, mm")
ylabel("total forces")
title("Gels")
ylim([0, 0.6])
xlim([0 1.3])

sgtitle("Amplitudes and Forces Across Measurements for 200 gram profiles")