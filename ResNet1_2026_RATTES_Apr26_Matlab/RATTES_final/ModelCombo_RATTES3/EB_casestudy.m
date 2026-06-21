function [fignum] = Eckland_Figure_EBcasestudy(data5all,calc5all,convert1,convert2,fignum)

%% CASE STUDY OF RESULTS FOR RIO GRANDE SPECIFIC RESERVOIRS

%close all
%set up data

calc=calc5all;
data=data5all;

eb_sid = 17114;
val=find(data.SID==eb_sid);
t1900=find(calc.t==1900);
t2025=find(calc.t==2025);
filename='EB_casestudy/EBR_Capacity_1916-2017.csv';



%% IMPORT ELEPHANT BUTTE SURVEY DATA

opts = detectImportOptions(filename);

opts = setvaropts(opts, 'SurveyDate', ...
    'InputFormat', 'MM/dd/yyyy');

eb_df = readtable(filename, opts);

eb_df.surveycap_m3 = eb_df.Active_capacity_acre_feet * convert1;
eb_df.surveysed_m3 = eb_df.Sediment_volume_acre_feet * convert1;

% extract survey year
eb_df.surveyyear = year(eb_df.SurveyDate);

%% SUBSET ELEPHANT BUTTE RESULTS
%year_cols = string(1900:2050);

cap_vals = calc.cap(val,:);
sed_vals = calc.sed(val,:);
trap_vals = calc.trap(val,:);
trap_brunevals=calc.trapbrune(val,:);
years=calc.t(:);

%% inflow- elephant butte was not in resops
%
Qresops_af=calc.resopsQm3(val,:)/convert1;
Qsuridx=NaN(size(eb_df.surveyyear));

% Qinflow Reclamation from Reclamation
filename2='EB_casestudy/Inflow_Reclamation.csv';

opts2 = detectImportOptions(filename2);

Qusbr_tbl = readtable(filename2, opts2);

% convert table columns to numeric arrays
tusbr = Qusbr_tbl{:,1};   % years
Qusbr_af = Qusbr_tbl{:,2}; % annual inflow acre-ft

%get annual Q for survey years
for i=1:numel(eb_df.surveyyear)
    Qsuridx(i)=find(tusbr==eb_df.surveyyear(i));
end
Qsur=Qusbr_af(Qsuridx);


%% survey trap brown brune
kappa=0.1; %silt

sur_brown=1-1./(1+calc.kappa(val)*(((eb_df.surveycap_m3(:))/convert1)./((data.DA(val))/convert2)));

%constants for brune
a=97;
b=6.42;
V=eb_df.Active_capacity_acre_feet;
tau=V./Qsur; % in af

sur_brune=(a*(1-2*exp(-b*(tau.^0.35))))/100;


cap_valsaf=cap_vals(t1900:t2025)/convert1;
tauusbr=cap_valsaf./(Qusbr_af');
trap_usbrbrune=(a*(1-2*exp(-b*(tauusbr.^0.35))))/100;


%% PLOT STACKED FIGURE

figure(fignum)
clf

% Figure formatting

% White background
set(gcf,'Color','w')

% Figure size in inches
set(gcf,'Units','inches')
set(gcf,'Position',[1 1 6.5 6.3])

% Recenter on screen
movegui(gcf,'center')

% Arial fonts
set(groot,'DefaultAxesFontName','Arial')
set(groot,'DefaultTextFontName','Arial')
set(groot,'DefaultAxesFontSize',10)
set(groot,'DefaultTextFontSize',10)

% Create tiled layout

t = tiledlayout(3,1, ...
    "TileSpacing","compact", ...
    "Padding","loose");

% Colors

cap_color  = "blue";
sed_color  = [0.55 0.27 0.07];   % brown
trap_color = [0 0.39 0];         % dark green
vline_color = [0.4 0.4 0.4];

% ------------------------------------------------------------------------
% Capacity plot

ax1 = nexttile;

h1 = plot(years(t1900:end), ...
    cap_vals(t1900:end)/1e6, ...
    "b-", ...
    "LineWidth",2, ...
    "DisplayName","Model");

hold on

h2 = xline(1973, "--", ...
    "Color",vline_color, ...
    "LineWidth",2, ...
    "DisplayName","Cochiti Dam Closure");

h3 = scatter(eb_df.surveyyear, ...
    eb_df.surveycap_m3/1e6, ...
    25, ...
    cap_color, ...
    "filled", ...
    "DisplayName","Survey");

text(1975, ...
    max(cap_vals/1e6,[],"omitnan"), ...
    "Cochiti Dam closure", ...
    "HorizontalAlignment","left", ...
    "VerticalAlignment","top", ...
    "FontSize",10);

ylabel('Capacity (M m^3)', ...
    'FontWeight','bold');

title('Case Study: Elephant Butte Reservoir, NM, USA', ...
    'FontWeight','bold');

grid on

set(gca,"XTickLabel",[])

hold off

% ------------------------------------------------------------------------
% Sediment plot

ax2 = nexttile;

plot(years(t1900:end), ...
    sed_vals(t1900:end)/1e6, ...
    "-", ...
    "Color",sed_color, ...
    "LineWidth",2);

hold on

xline(1973, "--", ...
    "Color",vline_color, ...
    "LineWidth",2);

scatter(eb_df.surveyyear, ...
    eb_df.surveysed_m3/1e6, ...
    25, ...
    sed_color, ...
    "filled");

ylabel('Sediment volume (M m^3)', ...
    'FontWeight','bold');

ylim([0 950])

grid on

set(gca,"XTickLabel",[])

hold off

% Trap efficiency plot

ax3 = nexttile;

h4 = plot(years(t1900:end), ...
    trap_vals(t1900:end)*100, ...
    "k-", ...
    "LineWidth",2, ...
    "DisplayName","Brown TE + model");

hold on

h5 = plot(tusbr, ...
    trap_usbrbrune*100, ...
    "-", ...
    "Color",'g', ...
    "LineWidth",2, ...
    "DisplayName","Brune TE, Reclamation + model");

burnt_orange = [0.85 0.33 0.10];

h6 = plot(years(t1900:end), ...
    trap_brunevals(t1900:end)*100, ...
    "-", ...
    "Color",burnt_orange, ...
    "LineWidth",2, ...
    "DisplayName","Brune TE, ResOps + model");

xline(1973, "--", ...
    "Color",vline_color, ...
    "LineWidth",2);

h7 = scatter(eb_df.surveyyear, ...
    sur_brune*100, ...
    25, ...
    'go', ...
    'MarkerFaceColor','g', ...
    'DisplayName','Brune TE, Reclamation + Survey');

h8 = scatter(eb_df.surveyyear, ...
    sur_brown*100, ...
    25, ...
    'k', ...
    'filled', ...
    'DisplayName','Brown TE, Survey');

ylabel('Trap efficiency (%)', ...
    'FontWeight','bold');

xlabel('Year', ...
    'FontWeight','bold');

ylim([85 100])

grid on

hold off

% ------------------------------------------------------------------------
% Combined legend below all panels

lgd = legend( ...
    [h1 h3 h2 h4 h5 h6 h7 h8], ...
    { ...
    'Model', ...
    'Survey', ...
    'Cochiti Dam Closure', ...
    'Brown TE + model', ...
    'Brune TE, Reclamation + model', ...
    'Brune TE, ResOps + model', ...
    'Brune TE, Reclamation + Survey', ...
    'Brown TE, Survey'}, ...
    'Orientation','horizontal');

lgd.Layout.Tile = 'south';

% Legend formatting
lgd.NumColumns = 2;
lgd.Box = 'off';
lgd.FontSize = 8;

% Increment figure number

fignum = fignum + 1;
% SAVE FIGURE


end