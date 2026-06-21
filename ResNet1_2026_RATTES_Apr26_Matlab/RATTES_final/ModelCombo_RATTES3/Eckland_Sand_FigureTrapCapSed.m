function [fignum] = Eckland_Sand_FigureTrapCapSed(data5sand,calc5sand,fignum)


%% Eckland et al. figures
data = data5sand;
calc = calc5sand;

mask_rattes = (data.inMLR == 1 | (data.issite == 1 & data.PermStorag == 1));
rattes = find(mask_rattes);

%% Figure 2 first pane; set up
t1900 = find(calc.t == 1900);
numt  = numel(calc.t);

tfig1 = calc.t(t1900:numt);
nt    = numel(tfig1);

%%
%capacity bands
capup = zeros(nt,1);
caplo = zeros(nt,1);
cap   = zeros(nt,1);

%sed bands
sedup = zeros(nt,1);
sedlo = zeros(nt,1);
sed   = zeros(nt,1);

design   = zeros(nt,1);

for i = 1:nt
    n = i + t1900 - 1;
    capup(i) = sum(calc.capup(rattes,n))/1e9;
    caplo(i) = sum(calc.caplo(rattes,n))/1e9;
    cap(i)   = sum(calc.cap(rattes,n))/1e9;

    sedup(i) = sum(calc.sedup(rattes,n))/1e9;
    sedlo(i) = sum(calc.sedlo(rattes,n))/1e9;
    sed(i)   = sum(calc.sed(rattes,n))/1e9;

    design(i)=sum(calc.capdesign(rattes,n))/1e9;
end


%% Fig 2 new panel
% identify small and large reservoirs to subset trap efficiency

sm=6.17e6;
mid=100e6;
canSID=500000;

mask_usa=(data.SID<canSID); %identifys dams in USA
mask_rattes = (data.inMLR == 1 | (data.issite == 1 & data.PermStorag == 1)); %rattes is us only w/ perm storage (inMLR does not include canada)- is either in MLR model or is a site & also has perm storage
mask_sm= (data.capcALL<sm); % find small reservoirs
mask_mid= (data.capcALL>=sm & data.capcALL<mid); % find mid-size reservoirs
mask_lg= (data.capcALL>=mid);

trap_sm=zeros(nt,1); 
trap_mid=zeros(nt,1);
trap_lg=zeros(nt,1);

count_sm=zeros(nt,1);
count_mid=zeros(nt,1);
count_lg=zeros(nt,1);

trap_smup=zeros(nt,1); 
trap_midup=zeros(nt,1);
trap_lgup=zeros(nt,1);

trap_smlo=zeros(nt,1); 
trap_midlo=zeros(nt,1);
trap_lglo=zeros(nt,1);

% get time series of median trap efficiency
for i=1:nt
    n = i + t1900 - 1;
    mask_exists = (data.yrc<=calc.t(n) & data.yrr>calc.t(n)); % make sure the dam exists
    smtemp=find(mask_rattes & mask_usa & mask_sm & mask_exists);
    midtemp=find(mask_rattes & mask_usa & mask_mid & mask_exists);
    lgtemp=find(mask_rattes & mask_usa & mask_lg & mask_exists);
    
    if isempty(smtemp)==0
        trap_sm(i)=mean(calc.trap(smtemp,n));
        trap_smup(i)=mean(calc.trapup(smtemp,n));
        trap_smlo(i)=mean(calc.traplo(smtemp,n));
        count_sm(i)=numel(smtemp);
    end
    if isempty(midtemp)==0
        trap_mid(i)=mean(calc.trap(midtemp,n));
        trap_midup(i)=mean(calc.trapup(midtemp,n));
        trap_midlo(i)=mean(calc.traplo(midtemp,n));
        count_mid(i)=numel(midtemp);
    end
    if isempty(lgtemp)==0
        trap_lg(i)=mean(calc.trap(lgtemp,n));
        trap_lgup(i)=mean(calc.trapup(lgtemp,n));
        trap_lglo(i)=mean(calc.traplo(lgtemp,n));
        count_lg(i)=numel(lgtemp);
    end
end


%% combine into a two panel figure.  
% Combined Figure: Volume Trends (Left) and Trap Efficiency (Right)

% Combined Figure: Volume Trends (Left) and Trap Efficiency (Right)
% This script generates a 2-panel figure for publication at 6.5 inches wide.
% Tighter fit between panels achieved via TileSpacing and Padding adjustments.

fig = figure(fignum); clf;
set(fig, 'Units', 'inches', 'Position', [1, 1, 6.5, 4.5], 'Color', 'w');

% Use 'compact' padding to keep left/right tight
tlo = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- PANEL A: Volume Trends (Left) ---
nexttile
hold on;
XData1 = [tfig1(:); flipud(tfig1(:))];
YData1 = [capup(:); flipud(caplo(:))];
YData2 = [sedup(:); flipud(sedlo(:))];

% Background bands (Confidence Intervals)
patch('XData', XData1, 'YData', YData1, 'FaceAlpha', 0.3, ...
      'FaceColor', [0.25, 0.41, 0.88], 'EdgeColor', 'none');
patch('XData', XData1, 'YData', YData2, 'FaceAlpha', 0.3, ...
      'FaceColor', [0.82, 0.71, 0.55], 'EdgeColor', 'none');

% Plot main lines (Design, Model Capacity, Model Sediment)
hDesign = plot(tfig1, design, 'Color', [0 0 0], 'LineWidth', 1.5); 
hModelCap = plot(tfig1, cap, 'Color', [0 0 0.5], 'LineStyle', '--', 'LineWidth', 1.5);
hModelSed = plot(tfig1, sed, 'Color', [0.54, 0.27, 0.07], 'LineStyle', '--', 'LineWidth', 1.5);

% Volume Labels (Using specified manual coordinates)
text(1997.98, 725.07, 'water capacity', 'Color', [0.04, 0.27, 0.42], ...
     'FontWeight', 'bold', 'FontName', 'Arial', 'FontSize', 8, 'BackgroundColor', 'none');
text(2009.75, 137.39, 'sediment', 'Color', [0.51, 0.41, 0.19], ...
     'FontWeight', 'bold', 'FontName', 'Arial', 'FontSize', 8, 'BackgroundColor', 'none');

ylabel('Volume, km^3', 'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold');
xlabel('Year', 'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold');
xlim([1900 2050]); ylim([0 850]); box on; grid on;
set(gca, 'FontName', 'Arial', 'FontSize', 8, 'GridLineStyle', ':', 'TickDir', 'out', 'Layer', 'top');

% Left Legend (Now reduced to 3 primary items)
hConfDummy = patch(nan, nan, [0.7 0.7 0.7], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
leg1 = legend([hDesign, hModelCap, hConfDummy], ...
    {'Design', 'Model', '95% Conf'}, ...
    'Location', 'northwest', 'FontSize', 8, 'FontName', 'Arial', 'Box', 'on');
set(leg1, 'ItemTokenSize', [18, 8]);

%% --- PANEL B: Trap Efficiency (Right) ---
nexttile
hold on;

annot_years = [1900 1950 2000 2050];

colors_trap = [0.10, 0.62, 0.46;   % small
               0.46, 0.44, 0.70;   % medium
               0.85, 0.37, 0.01];  % large

group_names = {'Small', 'Medium', 'Large'};

% convert to percent
data_series  = {trap_sm*100, trap_mid*100, trap_lg*100};
up_series    = {trap_smup*100, trap_midup*100, trap_lgup*100};
lo_series    = {trap_smlo*100, trap_midlo*100, trap_lglo*100};

count_series = {count_sm, count_mid, count_lg};

% --- confidence interval patches ---
for g = 1:length(data_series)

    Xpatch = [tfig1(:); flipud(tfig1(:))];
    Ypatch = [up_series{g}(:); flipud(lo_series{g}(:))];

    patch('XData', Xpatch, ...
          'YData', Ypatch, ...
          'FaceColor', colors_trap(g,:), ...
          'FaceAlpha', 0.20, ...
          'EdgeColor', 'none');
end

% --- main lines ---
for g = 1:length(data_series)

    p_trap(g) = plot(tfig1, data_series{g}, ...
                     'Color', colors_trap(g,:), ...
                     'LineWidth', 2, ...
                     'DisplayName', group_names{g});

    for yr = annot_years

        t_idx = find(tfig1 == yr);

        if ~isempty(t_idx) && ~isnan(data_series{g}(t_idx))

            val = count_series{g}(t_idx);

            % label location
            x_p = tfig1(t_idx);
            if g==3
              y_p = data_series{g}(t_idx) - 3;
            else
              y_p = data_series{g}(t_idx) - 5;
            end

            if yr == 1900
                h_al = 'left';
            elseif yr == 2050
                h_al = 'right';
            else
                h_al = 'center';
            end

            text(x_p, y_p, num2str(val), ...
                'FontName', 'Arial', ...
                'FontSize', 8, ...
                'FontWeight', 'bold', ...
                'Color', colors_trap(g,:), ...
                'HorizontalAlignment', h_al, ...
                'VerticalAlignment', 'top');
        end
    end
end

xlabel('Year', ...
       'FontName', 'Arial', ...
       'FontSize', 10, ...
       'FontWeight', 'bold');

ylabel('Trap Efficiency (%)', ...
       'FontName', 'Arial', ...
       'FontSize', 10, ...
       'FontWeight', 'bold');

xlim([1900 2050]);
ylim([15 100]);
yticks(15:10:95);

box on;
grid on;

set(gca, 'FontName', 'Arial', ...
         'FontSize', 8, ...
         'GridLineStyle', ':', ...
         'TickDir', 'out', ...
         'Layer', 'top');

% legend
hConfDummy = patch(nan, nan, [0.7 0.7 0.7], ...
                   'FaceAlpha', 0.2, ...
                   'EdgeColor', 'none');

leg2 = legend([p_trap hConfDummy], ...
              {'Small', 'Medium', 'Large', '95% Conf'}, ...
              'Location', 'southwest', ...
              'FontSize', 8, ...
              'FontName', 'Arial');
set(leg2, 'ItemTokenSize', [18, 8]);

title(leg2, 'Reservoir Size', ...
      'FontWeight', 'bold', ...
      'FontSize', 8);


%% --- Final Title and Margin Check ---
sgtitle('RATTES modeled capacity (sand), sediment volume, and trap efficiency', ...
    'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');
tlo.OuterPosition = [0 0 1 0.95]; 

fignum = fignum + 1;

end
