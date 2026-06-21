function [fignum] = Eckland_FigureBrune(data5all,calc5all,fignum)
% figure Brune


data=data5all;
calc=calc5all;
numdam=numel(data.SID);
idxbrune=[];

for i=1:numdam
    if (max(calc.trapbrune(i,:))>0)
    idxbrune=vertcat(idxbrune,i);
    end
end

TEbr=calc.trap(idxbrune,:);
TEbu=calc.trapbrune(idxbrune,:);
TEbr(isnan(TEbu)) = NaN; %makes TEbr nan where TEbu is nan so the means are using the same years for comparison
tauactive=calc.tauactive(idxbrune,:);

TEbr_mean = mean(TEbr, 2, 'omitnan') * 100;
TEbu_mean = mean(TEbu, 2, 'omitnan') * 100;


%% plot figure and generate table
figure(fignum); clf

% ---- figure appearance ----
set(gcf, 'Units', 'inches', 'Position', [1 1 6.5 3], ...
    'Color', 'w')   % white background

set(groot, 'DefaultAxesFontName', 'Arial')
set(groot, 'DefaultTextFontName', 'Arial')

tiledlayout(1,2, 'TileSpacing','compact', 'Padding','compact')

x = TEbr_mean;
y = TEbu_mean;

% classify using mean tau per site
tau_meanactive = mean(tauactive, 2, 'omitnan');
a=isnan(tau_meanactive); b=find(a==1);
tau_meanactive(a)=1; %there are 2 missing live storage.  Assume seasonal.
clear a b


idx_holdover = tau_meanactive > 1;
idx_seasonal = tau_meanactive <= 1;

%% ------------------------------------------------------------------------
% CREATE OUTPUT TABLE
% -------------------------------------------------------------------------

% percent difference: positive means Brune > Brown
diff_pct = TEbu_mean - TEbr_mean;

ReservoirType = strings(size(diff_pct));
ReservoirType(idx_holdover) = "Holdover";
ReservoirType(idx_seasonal) = "Seasonal";

ResultsTable = table( ...
    TEbr_mean, ...
    TEbu_mean, ...
    diff_pct, ...
    tau_meanactive, ...
    ReservoirType, ...
    'VariableNames', ...
    {'TE_Brown_pct','TE_Brune_pct','Diff_BruneMinusBrown_pct', ...
     'TauActive_mean','ReservoirType'});

disp(ResultsTable)


% ------------------------------------------------------------------------
% CATEGORY DEFINITIONS
% -------------------------------------------------------------------------

category_names = { ...
    '>25% higher', ...
    '10-25% higher', ...
    '5-10% higher', ...
    '0-5% higher', ...
    '0-5% lower', ...
    '5-10% lower', ...
    '10-25% lower', ...
    '>25% lower'};

count_categories = @(d) [ ...
    sum(d > 25), ...
    sum(d > 10 & d <= 25), ...
    sum(d > 5  & d <= 10), ...
    sum(d >= 0 & d <= 5), ...
    sum(d < 0  & d >= -5), ...
    sum(d < -5 & d >= -10), ...
    sum(d < -10 & d >= -25), ...
    sum(d < -25) ];


% ------------------------------------------------------------------------
% HOLDOVER COUNTS
% -------------------------------------------------------------------------

diff_holdover = diff_pct(idx_holdover);

holdover_counts = count_categories(diff_holdover);

HoldoverSummary = table( ...
    category_names', ...
    holdover_counts', ...
    'VariableNames', {'DifferenceCategory','Count'});

disp(' ')
disp('HOLDOVER RESERVOIRS')
disp(HoldoverSummary)


% ------------------------------------------------------------------------
% SEASONAL COUNTS
% -------------------------------------------------------------------------

diff_seasonal = diff_pct(idx_seasonal);

seasonal_counts = count_categories(diff_seasonal);

SeasonalSummary = table( ...
    category_names', ...
    seasonal_counts', ...
    'VariableNames', {'DifferenceCategory','Count'});

disp(' ')
disp('SEASONAL RESERVOIRS')
disp(SeasonalSummary)


%% back to table

% ---- LEFT: FULL VIEW ----
nexttile; hold on

plot(x(idx_holdover), y(idx_holdover), ...
    's', 'Color', 'b', 'MarkerFaceColor', 'none')

plot(x(idx_seasonal), y(idx_seasonal), ...
    '^', 'Color', 'k', 'MarkerFaceColor', 'none')

plot([0 100], [0 100], 'r--', 'LineWidth', 1)


xlim([30 100])
ylim([30 100])
axis square

xlabel('Brown Trap Efficiency (%)')
ylabel('Brune Trap Efficiency (%)')

title('Full Range')

grid on
box on

legend('Holdover (\tau > 1)', 'Seasonal (\tau \leq 1)', '1:1 line', ...
    'Location', 'southeast')

hold off


% ---- RIGHT: ZOOMED (80–100%) ----
nexttile; hold on

plot(x(idx_holdover), y(idx_holdover), ...
    's', 'Color', 'b', 'MarkerFaceColor', 'none')

plot(x(idx_seasonal), y(idx_seasonal), ...
    '^', 'Color', 'k', 'MarkerFaceColor', 'none')

plot([90 98], [90 98], 'r--', 'LineWidth', 1)

xlim([90 98])
ylim([90 98])
axis square

xlabel('Brown Trap Efficiency (%)')
ylabel('Brune Trap Efficiency (%)')

title('Zoom: 80–100%')

grid on
box on

hold off

exportgraphics(gcf, fullfile('Figures', ...
    'Eckland_BruneBrown.png'), ...
    'Resolution', 600)

savefig(gcf, fullfile('Figures', ...
    'Eckland_BruneBrown.fig'))

fignum=fignum+1;

end