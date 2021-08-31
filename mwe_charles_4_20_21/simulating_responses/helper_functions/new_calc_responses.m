function [FRs_ts, FRs_gel] = new_calc_responses(r_in, plot_flag,...
    names)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%FLIP r_in for this function - goes ts then gel
r_in = {r_in{2}, r_in{1}};
response_fig = figure;
for i = 1:2
    resp = r_in{i};
    if plot_flag
        figure(response_fig);
        subplot(2,3,i);
        plot(resp);
        title(names(i));
        
        ax = gca;
        ax.FontSize = 12;
        ax.FontWeight = 'bold';
        if i == 2
            ax.YAxis.Visible = 'off'; % remove y-axis
        end
    end
    
    %calculate mean frs and sd
    
    FRs{i,1} = resp.rate(resp.affpop.iPC);
    FRs{i,2} = resp.rate(resp.affpop.iRA);
    FRs{i,3} = resp.rate(resp.affpop.iSA1);
    means = zeros(3,1);
    sem = means;
    
    for j = 1:3
        means(j)  = mean(FRs{i,j});
        sem(j) = std(FRs{i,j})/sqrt(length(FRs{i,j}));
    end
    
    FRs{i,4} = means;
    FRs{i, 5} = sem;
    x = [1,2,3];
    
    if plot_flag
        subplot(2,3,i+3);
        bar(x, means);
        hold on
        b = bar(x,means);
        b.FaceColor = 'flat';
        b.CData(1,:) =  [255 127 0]/255;
        b.CData(2,:) =  [30 120 180]/255;
        b.CData(3,:) = [50 160 40]/255;
        er = errorbar(x,means,sem);
        er.Color = [0 0 0];
        er.LineStyle = 'none';
        %         title(strcat(ts_structs(i).name, " firing rate"));
        xticks(x)
        xticklabels({'PCs','RAs','SAs'})
        ylabel("Hz")
        ax = gca;
        ax.FontSize = 12;
        ax.FontWeight = 'bold';
    end
end

FRs_ts = {};
FRs_gel = {};

for i = 1:length(FRs)
    FRs_ts{i} = FRs{1,i};
    FRs_gel{i} = FRs{2,i};
end



end

