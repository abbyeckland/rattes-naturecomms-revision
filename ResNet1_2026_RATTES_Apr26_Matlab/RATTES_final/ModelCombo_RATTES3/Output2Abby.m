% export output to Abby

calc=calc5all;
data=data5all;

numdam=numel(data.SID);
numt=numel(calc.t);

cap_timeseries=NaN((numdam+1),(numt+1));
cap_timeseries(1,2:end)=calc.t;
cap_timeseries(2:end,1)=data.SID;
sed_timeseries=cap_timeseries;
trap_timeseries=cap_timeseries;

cap_timeseries(2:end,2:end)=calc.cap(:,:);
sed_timeseries(2:end,2:end)=calc.sed(:,:);
trap_timeseries(2:end,2:end)=calc.trap(:,:);

capfilename='RATTEScap.csv';
trapfilename='RATTEStrap.csv';
sedfilename='RATTESsed.csv';
datafilename='RATTESdatafinal.csv';

writematrix(cap_timeseries,capfilename,'Delimiter',',');
writematrix(trap_timeseries,trapfilename,'Delimiter',',');
writematrix(sed_timeseries,sedfilename,'Delimiter',',');
writetable(data, datafilename,'Delimiter',',');


