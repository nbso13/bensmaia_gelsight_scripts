clear
close all
cd activities
load("compliant_activities_31-Mar-2021_200_0_20_0_7_80_0.8.mat")
cd ..
motherOfAllPlotsFunc(activities)
sgtitle("200 Grams, top 20 neurons")

cd activities
load("compliant_activities_31-Mar-2021_200_0_30_0_7_80_0.8.mat")
cd ..
motherOfAllPlotsFunc(activities)
sgtitle("200 Grams, top 30 neurons")

cd activities
load("compliant_textures_200_grams.mat")
cd ..
motherOfAllPlotsFunc(activities)
sgtitle("200 Grams, top 10 neurons")

cd activities
load("march_10_final_activities_100_grams")
cd ..
activities.names = activities.names(1:5); %only looking at compliant ones
activities.real = activities.real(1:5,:);
activities.ts = activities.ts(1:5, :);
activities.gel = activities.gel(1:5, :);
motherOfAllPlotsFunc(activities)
sgtitle("100 Grams, top 10 neurons")

cd activities
load("noncompliant_activities_03-Apr-2021_200_0_20_0_7_80_0.8.mat")
cd ..
motherOfAllPlotsFunc(activities)
sgtitle("Noncompliant, 200 Grams, top 20 neurons")

cd activities
load("march_10_final_activities_100_grams")
cd ..
activities.names = activities.names(6:8); %only looking at compliant ones
activities.real = activities.real(6:8,:);
activities.ts = activities.ts(6:8, :);
activities.gel = activities.gel(6:8, :);
motherOfAllPlotsFunc(activities)
sgtitle("Noncompliant, 100 Grams, top 10 neurons")
