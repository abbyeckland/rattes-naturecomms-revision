function [calc5all] = Brune(data5all,calc5all,convert1)

%Brune_Brown

%% load data
data=data5all;
calc=calc5all;

storagefile="ResOpsUS/time_series_single_variable_table/DAILY_AV_STORAGE_MCM.csv";
opts = detectImportOptions(storagefile);
opts.VariableNamesLine = 1;     % header row
opts = setvaropts(opts, opts.VariableNames, 'QuoteRule','remove'); % remove '
temp = readtable(storagefile, opts);
temp{:,1}=datetime(temp{:,1}, 'InputFormat','yyyy-MM-dd');
for i = 2:width(temp)
    if iscellstr(temp{:,i}) || isstring(temp{:,i})
        temp.(i) = str2double(temp.(i));
    end
end
datestore=temp.date(2:end);
residstores=temp(1,2:end);
storem3 = temp{2:end, 2:end} * 1e6;
clear i opts storagefile temp

inflowfile="ResOpsUS/time_series_single_variable_table/DAILY_AV_INFLOW_CUMECS.csv";
opts = detectImportOptions(inflowfile);
opts.VariableNamesLine = 1;     % header row
opts = setvaropts(opts, opts.VariableNames, 'QuoteRule','remove'); % remove '
temp = readtable(inflowfile, opts);
temp{:,1}=datetime(temp{:,1}, 'InputFormat','yyyy-MM-dd');
for i = 2:width(temp)
    if iscellstr(temp{:,i}) || isstring(temp{:,i})
        temp.(i) = str2double(temp.(i));
    end
end

dateinflow=temp.date(2:end);
residinflow=temp(1,2:end);
inflowm3 = temp{2:end, 2:end};
clear i opts inflowfile temp




%% Inflow calcs (annual volume with gap filling for years that have at least 300 days worth of data)

yr = year(dateinflow);

years = unique(yr);
ny = numel(years);
nsites = size(inflowm3,2);

annual_inflow = NaN(ny, nsites);   % m³/year

for k = 1:ny
    idx = yr == years(k);   % rows for this year- how many days of data for each year
    if numel(idx)<365
        disp('not enough days in the year')
        years(k)
    end
    
    Q = inflowm3(idx,:);    % m³/s (days × sites)
    
    % Convert to daily volume
    Vday = Q * 86400;       % m³/day
    
    % Count valid days
    nvalid = sum(~isnan(Vday), 1);
    
    % Mean daily volume (ignore NaNs)
    Vmean = mean(Vday, 1, 'omitnan');   % m³/day
    
    % Fill missing days with mean
    Vday_filled = Vday;
    for j = 1:nsites
        if nvalid(j) >= 300
            missing = isnan(Vday(:,j));
            Vday_filled(missing,j) = Vmean(j);
            
            % Sum full year
            annual_inflow(k,j) = sum(Vday_filled(:,j));
        else
            annual_inflow(k,j) = NaN;
        end
    end
end

clear Q Vday Vmean Vday_filled idx k j nsites ny nvalid yr years

%% annual storage

yr = year(datestore);

years = unique(yr);
ny = numel(years);
nsites = size(storem3,2);

annual_storage = NaN(ny, nsites);   % m³/year

for k = 1:ny
    idx = yr == years(k);   % rows for this year- how many days of data for each year
    if numel(idx)<365
        disp('not enough days in the year')
        years(k)
    end
    
    Vday = storem3(idx,:);    % m³ per day (days × sites)
      
    % Count valid days
    nvalid = sum(~isnan(Vday), 1);
    
    % Mean daily storage volume (ignore NaNs)
    Vmean = mean(Vday, 1, 'omitnan');   % m³/day
    
    % Fill missing days with mean
    Vday_filled = Vday;
    for j = 1:nsites
        if nvalid(j) >= 300
            missing = isnan(Vday(:,j));
            Vday_filled(missing,j) = Vmean(j);
            
            % mean annual storage
            annual_storage(k,j) = mean(Vday_filled(:,j));
        else
            annual_storage(k,j) = NaN;
        end
    end
end

clear Q Vday Vmean Vday_filled idx k j nsites ny nvalid yr


%% already ensured that the IDS match in the inflow and in the storage files
%ROID == GID
ids_inflow = string(residinflow{1,:});
ROID = ids_inflow;     % keep one  ROID= ResOps ID
ROID = str2double(ROID);
idxRATTES=NaN(size(ROID));

%ResOps ID is the same as GID (GRanD)
for i=1:numel(ROID)
    tmp=find(data.GID==ROID(i) & data.PermStorag==1 & data.SID<500000);
    if isempty(tmp)~=1
        idxRATTES(i)=tmp;
    end
end

clear ids_inflow

% limit analysis to sites that have annual inflow

check=NaN(size(ROID));
for i=1:numel(ROID)
    check(i)=max(annual_inflow(:,i));
end
a=isnan(check);
delete=find(a==1);

%delete data that doesn't have any inflow
annual_inflow(:,delete)=[];
annual_storage(:,delete)=[];
idxRATTES(delete)=[];
ROID(delete)=[];

%delete data without a RESNET crossref
a=isnan(idxRATTES);
delete=find(a==1);
%delete data 
annual_inflow(:,delete)=[];
annual_storage(:,delete)=[];
idxRATTES(delete)=[];
ROID(delete)=[];

clear a check delete i loadalldata 

%% calc Brune- actually want available storage from model... not live storage.
TEbu=NaN(size(annual_inflow));
tau=NaN(size(annual_inflow));
tauactive=NaN(size(annual_inflow));

%constants for brune
a=97;
b=6.42;

yr_resops=years;

idxtRATTES=NaN(size(yr_resops));
for j=1:numel(yr_resops)
    idxtRATTES(j)=find(calc.t==yr_resops(j));
end

calc.trapbrune=NaN(size(calc.trap));
calc.resopsQm3=NaN(size(calc.trap));
calc.resops_capactivem3=NaN(size(calc.trap));
calc.tauactive=NaN(size(calc.trap));

count0=0;

for i= 1:numel(ROID)
    for j=1:numel(yr_resops)
        if (annual_inflow(j,i))>0
            V=calc.cap(idxRATTES(i),idxtRATTES(j))/convert1; %available capacity for that year in AF
            Q=annual_inflow(j,i)/convert1; %inflow for taht year in AF
            Vactive=annual_storage(j,i)/convert1; %live capacity that year in AF
            tau(j,i)=V/Q;
            tauactive(j,i)=Vactive/Q;
            TEbu(j,i)=(a*(1-2*exp(-b*(tau(j,i)^0.35))))/100;
            if TEbu(j,i)>1
                disp('TE greater than 1')
                return
            elseif TEbu(j,i)<0
                TEbu(j,i)=0;
                count0=count0+1;
            end
            calc.trapbrune(idxRATTES(i),idxtRATTES(j))=TEbu(j,i);
            calc.resopsQm3(idxRATTES(i),idxtRATTES(j))=annual_inflow(j,i);
            calc.resops_capactivem3(idxRATTES(i),idxtRATTES(j))=annual_storage(j,i);
            calc.tauactive(idxRATTES(i),idxtRATTES(j))=tauactive(j,i);
        end
    end
end
disp(['Number of TE that were less than 0: ', num2str(count0)])


calc5all=calc;
save('Brunecalc.mat','calc5all','-v7.3');


disp('*')
disp('Done with Brune Calc')

end



