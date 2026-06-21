function [fignum] = Eckland_HistogramFig(data5all,calc5all,fignum)

% EVALUATE CAP LOSS FOR DIFFERENT SIZE RESERVOIRS

calc=calc5all;
data=data5all;
canSID=500000;


%% Define masks

now=find(calc.t==2025);

mask_rattes = (data5all.inMLR == 1 | (data5all.issite == 1 & data5all.PermStorag == 1)); %rattes is us only w/ perm storage (inMLR does not include canada)
mask_usa=(data.SID<canSID);
mask_notremoved= (data.yrr > 2050); %already confirmed no future dam removal dates past 2025

rattes25=find(mask_rattes & mask_usa & mask_notremoved);

%% ASSIGN SIZE CLASS
size_class = strings(height(rattes25), 1);
size_class(data.capcALL(rattes25) < 6.17e6) = "Small";
size_class(data.capcALL(rattes25) >= 6.17e6 & data.capcALL(rattes25) < 1e8) = "Medium";
size_class(data.capcALL(rattes25) >= 1e8) = "Large";
sizes=size_class;

%% PREP HISTOGRAM DATA FOR 2025 AND 2050
loss_2025 = calc.sed(rattes25,now)./calc.capdesign(rattes25,now)*100;
loss_2050 = calc.sed(rattes25,end)./calc.capdesign(rattes25,end)*100;

plot_df = table(data.SID(rattes25), sizes, loss_2025, loss_2050, ...
    'VariableNames', {'SID', 'size_class', 'pct_cap_loss_2025', 'pct_cap_loss_2050'});
%% PLOT HISTOGRAMS
figure(fignum)
tiledlayout(2, 2, "TileSpacing", "compact", "Padding", "loose");
size_order = ["All", "Small", "Medium", "Large"];
% colors (match your R plots)
color_2025 = [27, 158, 119] / 255;   % green
color_2050 = [217, 95, 2] / 255;     % orange
bin_edges = 0:5:100;

figure(fignum)
set(gcf, ...
    'Color', 'w', ...
    'Units', 'inches', ...
    'Position', [1 1 6.5 4.5]);   % 6.5 inches wide

for i = 1:length(size_order)
    nexttile;
    sz = size_order(i);
    if sz == "All"
        n_res = numel(plot_df.SID);
        idx = 1:n_res;
        idx=idx';
    else
        idx = find(size_class == sz);
        n_res=numel(idx);
    end
    % extract values
    vals_2025 = plot_df.pct_cap_loss_2025(idx);
    vals_2050 = plot_df.pct_cap_loss_2050(idx);
    
    % plot histograms
    h2025=histogram(vals_2025, bin_edges, ...
        "FaceColor", color_2025, ...
        "FaceAlpha", 0.5, ...
        "EdgeColor", "black");
    hold on;
    h2050=histogram(vals_2050, bin_edges, ...
        "FaceColor", color_2050, ...
        "FaceAlpha", 0.5, ...
        "EdgeColor", "black");
  
    % medians
    med_2025 = median(vals_2025);
    med_2050 = median(vals_2050);
    h2025med=xline(med_2025, "--", "Color", color_2025, "LineWidth", 2);
    h2050med=xline(med_2050, "--", "Color", color_2050, "LineWidth", 2);
    % subtitle
    nstr = regexprep(num2str(n_res,'%d'),'\d(?=(\d{3})+$)', '$0,');

    if sz == "All"
        ttl = sprintf('All (n = %s)', nstr);
    else
        ttl = sprintf('%s (n = %s)', char(sz), nstr);
    end

    title(ttl, ...
    'FontName', 'Arial', ...
    'FontSize', 11, ...
    'FontWeight', 'bold', ...
    'Units', 'normalized', ...
    'Position', [0.5 1.04 0]);

    % panel label
    panel_labels = ["A", "B", "C", "D"];
    if i <= 2
        label_y = 1.09;   % A, B
    else
        label_y = 1.07;   % C, D moved down slightly
    end
    text(-0.08, label_y, panel_labels(i), ...
        'Units', 'normalized', ...
        'FontWeight', 'bold', ...
        'FontSize', 12);
    ylabel("Count", ...
        'FontName','Arial', ...
        'FontSize',10);
    if i >= 3
        xlabel("Capacity loss (%)", ...
        'FontName','Arial', ...
        'FontSize',10);
    else
        set(gca, "XTickLabel", []);
    end
    xlim([0 100]);
    xticks(0:25:100);
    grid on;
    ax = gca;
    ax.TickDir = 'out';
    ax.FontName = 'Arial';
    ax.FontSize = 9;    % choose 9 or 10 as desired
    ax.Color = 'w';     % white axes background

    if i == 4
        overlap_color = [179 151 95]/255;
        hOverlap = patch(NaN,NaN,overlap_color, ...
            'EdgeColor','black');
        lgd = legend([h2025 h2050 hOverlap h2025med h2050med], ...
            {'2025','2050','2025 & 2050','2025 median','2050 median'}, ...
            'Location','northeast');
    end
    hold off;
end
set(gcf,'Color','w');

%% SAVE FIGURE
exportgraphics(gcf, fullfile('Figures', ...
    'Eckland_Historgram.png'), ...
    'Resolution', 600)

savefig(gcf, fullfile('Figures', ...
    'Eckland_Histogram.fig'))

fignum=fignum+1;
end