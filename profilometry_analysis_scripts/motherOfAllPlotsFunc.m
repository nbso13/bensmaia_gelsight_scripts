function [rmses] = motherOfAllPlotsFunc(activities)
%motherOfAllPlots takes in activity statistics from real afferents and
%simulated afferents from touchsim and gelsight and compares the simulation
%to the real values.

num_textures = length(activities.names);

rmses = zeros(2, 3);

aff_class = ["PCs", "RAs", "SA1s"];
figure;
scatter_list = [];
for i = 1:3 %for each afferent
    subplot(2, 3, i)
    hold on;
    title(strcat(aff_class(i), ", Gel"));
    for j = 1:num_textures
        scatter_list(j) = scatter(activities.real(j,i),activities.gel(j,i), 'filled');
        er = errorbar(activities.real(j,i), activities.gel(j,i), activities.real(j,i+3), 'horizontal'); %real
        er.Color = [0 0 0];
        er.LineStyle = 'none';
        er = errorbar(activities.real(j,i),activities.gel(j,i), activities.gel(j,i+3));
        er.Color = [0 0 0];
        er.LineStyle = 'none';
    end
    x = 1:max(max(activities.real(:,i)), max(activities.gel(:, i)));
    y = x;
    plot(x,y);
%     mdl = fitlm(activities.real(:,i),activities.gel(:,i));
%     [r, p] = corrcoef(activities.real(:,i),activities.gel(:,i));
    rmses(1,i) = rmse_calc(activities.real(:,i),activities.gel(:,i));
    r_sq_str = {strcat("RMSE = ", num2str(round(rmses(1,i))))};
    %strcat("p = ", num2str(round(p(1,2), 3)))
    text(x(round(end/2)), x(end)+5, r_sq_str, 'HorizontalAlignment','right');
    xlabel("Real Recorded Activity (Hz)");
    ylabel("Gel Simulated Activity (Hz)");
    if (i == 1)
    	legend(scatter_list, activities.names, 'Interpreter', 'none');
    end
end

scatter_list = [];
for i = 1:3
    subplot(2, 3, i+3)
    hold on;
    title(strcat(aff_class(i), ", TouchSim"));
    for j = 1:num_textures
        scatter_list(j) = scatter(activities.real(j,i), activities.ts(j,i), 'filled');
        er = errorbar(activities.real(j,i),activities.ts(j,i), activities.real(j,i+3), 'horizontal');
        er.Color = [0 0 0];
        er.LineStyle = 'none';
        er = errorbar(activities.real(j,i),activities.ts(j,i), activities.ts(j,i+3));
        er.Color = [0 0 0];
        er.LineStyle = 'none';
    end
    x = 1:max(max(activities.real(:,i)), max(activities.ts(:, i)));
    y = x;
    plot(x,y);
    mdl = fitlm(activities.real(:,i),activities.ts(:,i));
    [r, p] = corrcoef(activities.real(:,i),activities.ts(:,i));
    rmses(2, i) = rmse_calc(activities.real(:,i),activities.ts(:,i));
    r_sq_str = {strcat("RMSE = ", num2str(round(rmses(2, i))))};
    %strcat("p = ", num2str(round(p(1,2), 3)))
    text(x(round(end/2)), x(end)+5, r_sq_str, 'HorizontalAlignment','right');
    xlabel("Real Recorded Activity (Hz)");
    ylabel("TouchSim Simulated Activity (Hz)");
end

end

function [rmse] = rmse_calc(real_rates, sim_rates)
    rmse = sqrt(mean((real_rates-sim_rates).^2));
end
