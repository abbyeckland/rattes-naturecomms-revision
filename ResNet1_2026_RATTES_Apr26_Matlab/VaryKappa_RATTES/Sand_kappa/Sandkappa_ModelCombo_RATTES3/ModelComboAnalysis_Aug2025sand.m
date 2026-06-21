%Combine MLR and SY models, back project MLR results to year of dam
%completion if yrp>yrc, analysis
%figure and table creation
%figures will save into pre-set 'Figures' subfolders

clear
tic

load('../Sandkappa_SYModel_MLRinput_RATTES_1/SYModelResults_MLRinput_sand.mat')

loadresults2="no"; %load final results if don't need to rerun

if loadresults2=="yes"
    load FinalRATTESsand.mat
    clear loadresults2;
    return
end
clear loadresults2;
%% setup
inputdateMLR='rattes_v1p1_sand'; %MLR model run date

modelcombo="yes";
compare="yes";
riversetup="yes";
riverDA="yes";
termdeltadams="yes";

fignum=1; %figure number for plotting

%% conversions
%conversions,
convert1=1233.482; %converts m3 to ac-ft is convert1*AF=m3, or from m3 is m3/convert1=AF
convert2= 2.59; %converts km2 and mi2 is convert2*mi2= km2 or from km2 is km2/convert2=mi2

%% Calc MLR backwards and combine with site data to match original input dam length
%need to add 95% confidence interval going backwards

if modelcombo == "yes"
    fmt='mlr_trap_decimal_%s.csv';
    filenameMLRtrap=sprintf(fmt,inputdateMLR); clear fmt;
    fmt='mlr_cap_up_95_%s.csv';
    filenameMLRcapup=sprintf(fmt,inputdateMLR); clear fmt;
    fmt='mlr_cap_low_95_%s.csv';
    filenameMLRcaplow=sprintf(fmt,inputdateMLR); clear fmt;
    fmt='mlr_cap_%s.csv';
    filenameMLRcap=sprintf(fmt,inputdateMLR); clear fmt;
    fmt='mlr_sv_up_95_%s.csv';
    filenameMLRsvup=sprintf(fmt,inputdateMLR); clear fmt;
    fmt='mlr_sv_low_95_%s.csv';
    filenameMLRsvlow=sprintf(fmt,inputdateMLR); clear fmt;
    fmt='mlr_sv_%s.csv';
    filenameMLRsv=sprintf(fmt,inputdateMLR); clear fmt;
   [data5all,calc5all] = CombineMLRSY(filenameMLRtrap,filenameMLRcap,filenameMLRcapup,filenameMLRcaplow,filenameMLRsv,filenameMLRsvup,filenameMLRsvlow,data2norivers,calc2_MLRinput,data4SYmodel,calc4_SY2,convert1,convert2);
   clear filenameMLRtrap filenameMLRcap filenameMLRcapup filenameMLRcaplow filenameMLRsv filenameMLRsvup filenameMLRsvlow
else
    load ComboModel.mat
end
clear modelcombo;

%% compare model results
%at dams with perm storage in the US only. 
if compare=="yes"
    [fignum,calc5all,data5all] = CompareSY_MLRAug25(calc5all,data5all,fignum);
else
    load CompareModels.mat;
end
clear compare;

%% Calc sed contributing DA at rivers and all dams
%includes all dams, including lock dams without permanent storage.

if riversetup == "yes"
    [data1wrivers,calc1all_wrivers] = RiverSedDAinput(data1wrivers,calc5all,data5all);%
else
    load riversetup.mat
end
clear riversetup;

%% Calc sed contributing DA at rivers, and upstream sediment
%Sed contributing DA includes all dams, including lock dams without
%permanent storage in US for sediment and yield

if riverDA == "yes"
    [calc1all_wrivers,data1wrivers,data6_delta,calc6_delta,calc5all,fignum] = SDArivers(calc1all_wrivers,data1wrivers,calc5all,data5all,fignum);
    clear doneriver;
else
    load riverSDA.mat
end
clear riverDA

%% Terminal Dams & number of dams above rivers
% this is with permanent storage us dams
if termdeltadams=="yes"
    [calc1all_wrivers,calc6_delta] = DeltaTerm(calc1all_wrivers,calc6_delta,calc5all,data5all,data1wrivers,data6_delta);
else
    load deltaterm.mat
end
clear termdeltadams


%%
clear inputdateMLR fignum
calc1sand=calc1all_wrivers;
calc2sand=calc2_MLRinput;
calc3sand=calc3_SY1;
calc4sand=calc4_SY2;
calc5sand=calc5all;
calc6sand=calc6_delta;
data1sand=data1wrivers;
data2sand=data2norivers;
data3sand=data3SYmodel;
data4sand=data4SYmodel;
data5sand=data5all;
data6sand=data6_delta;
clear calc1all_wrivers calc2_MLRinput calc3_SY1 calc4_SY2 calc5all calc6_delta data1wrivers data2norivers data3SYmodel data4SYmodel data5all data6_delta
save('FinalRATTESsand.mat','-v7.3');
save('RATTESsand5outputs_datacalc.mat','calc5sand','data5sand','-v7.3');


toc