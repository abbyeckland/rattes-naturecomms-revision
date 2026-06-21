function [calc7_OutHUC,data5all] = outletHUC(calc5all,data5all,convert1,convert2)
% sedimentation rate by outlet and HUC

%COUNTRY OUT: where terminal flowline goes: 11 canada to atlantic ocean; 12 canada to pacific ocean; 2 greatlakes; 31 mexico to atlantic ocean; 32 mexico to pacific ocean ; 
% 41 exits US to atlantic ocean; 42 exits US to pacific ocean; 5 exits US to gulf of mexico atlantic ocean; 0 no coast (internally drained)

%% Sediment and yield by outlet
data=data5all;
calc=calc5all;
now=find(calc.t==2025);
numt=numel(calc.t);

outlets=unique(data.countryOut);%listed in numerical order
outlets = sortrows(outlets);
numout=numel(outlets);

DAout=[440317.889179,785127.22263,4653818.86735,353115.19781,717924.41305,810967.695627,1199871.89105]; % in order by outlet #, 0,2,5,11,32,41,42
DAout=DAout';

calc7_OutHUC.sedshedaboveoutlet_km2=zeros(numout,numt);
calc7_OutHUC.seddtaboveoutlet_m3peryr=zeros(numout,numt);
calc7_OutHUC.sedaboveoutlet_km3=zeros(numout,numt);
calc7_OutHUC.sedyieldaboveoutlet_m3perkm2yr=zeros(numout,numt);%last year's sed shed influences this years sedimentation
calc7_OutHUC.percentDAtrappedaboveoutlet=zeros(numout,numt);
calc7_OutHUC.numdamsaboveoutlet_USpermstore=zeros(numout,numt);
calc7_OutHUC.capdesignaboveoutlet=zeros(numout,numt);
calc7_OutHUC.capaboveoutlet=zeros(numout,numt);
calc7_OutHUC.sedshedaboveOUTwCAN_km2=zeros(numout,numt);

DAHUC=[162439.896132,270002.074315,721855.541896,490027.460434,421965.9382,105948.369466,491924.10322,269755.371049,155636.055745,1322148.04891,642212.11842,472994.533693,...
    346736.665278,293568.761839,362987.524568,367048.650828,716517.941665,422292.247268]; %has canada area and mex

%for i=1:numout
for i=1:numout
    idx=find(data.countryOut==outlets(i) & data.PermStorag==1 & data.SID<500000);
    idx2=find(data.countryOut==outlets(i));%has canada b/c needed for percent trapped b/c nhd flowlines have some mex and canada percentage
    for j=2:numt
        calc7_OutHUC.numdamsaboveoutlet_USpermstore(i,j)=sum(calc.builtidx(idx)<=j);
        calc7_OutHUC.sedshedaboveoutlet_km2(i,j)=sum(calc.sedshedALL(idx,j));
        calc7_OutHUC.seddtaboveoutlet_m3peryr(i,j)=sum(calc.seddt(idx,j));
        calc7_OutHUC.sedaboveoutlet_km3(i,j)=sum(calc.sed(idx,j))/1e9;
        calc7_OutHUC.sedyieldaboveoutlet_m3perkm2yr(i,j)=calc7_OutHUC.seddtaboveoutlet_m3peryr(i,j)/calc7_OutHUC.sedshedaboveoutlet_km2(i,j-1);%last year's sed shed influences this years sedimentation
        calc7_OutHUC.capdesignaboveoutlet(i,j)=sum(calc.capdesign(idx,j));
        calc7_OutHUC.capaboveoutlet(i,j)=sum(calc.cap(idx,j));

        calc7_OutHUC.sedshedaboveOUTwCAN_km2(i,j)=sum(calc.sedshedALL(idx2,j));
        calc7_OutHUC.percentDAtrappedaboveoutlet(i,j)=calc7_OutHUC.sedshedaboveOUTwCAN_km2(i,j)/(DAout(i))*100;
    end
end

calc7_OutHUC.caplossaboveoutlet=(calc7_OutHUC.capdesignaboveoutlet-calc7_OutHUC.capaboveoutlet)./(calc7_OutHUC.capdesignaboveoutlet).*100;

calc7_OutHUC.SedTotalAboveOutlet_2025km3=calc7_OutHUC.sedaboveoutlet_km3(:,now);
calc7_OutHUC.SedTotalAboveOutlet_2050km3=calc7_OutHUC.sedaboveoutlet_km3(:,end);
calc7_OutHUC.YieldAboveOutlet_2025_m3perkm2yr=calc7_OutHUC.sedyieldaboveoutlet_m3perkm2yr(:,now);
calc7_OutHUC.YieldAboveOutlet_2025_AFpermi2yr=calc7_OutHUC.YieldAboveOutlet_2025_m3perkm2yr/convert1*convert2;
calc7_OutHUC.CapLossAboveOutlet_2025_percent=calc7_OutHUC.caplossaboveoutlet(:,now);
calc7_OutHUC.OutletID=outlets;
calc7_OutHUC.SeddtAboveOutlet_2025_Mm3peryer=calc7_OutHUC.seddtaboveoutlet_m3peryr(:,now)/1e6;

clear numout idx i a j

%% sediment and yield by HUC
filename='import_hydroNHD_082625.csv';
hydrodata=readtable(filename,'Format','auto');

hydrodata=renamevars(hydrodata,["Field1","Field2"],...
                           ["SID","HydroSeq_dam"]);
hydrodata.Field3=[];
hydrodata = sortrows(hydrodata,"SID");

%check that data fields match
check=data.SID-hydrodata.SID;
a=find(check~=0);
if ~isempty(a)
    disp('mismatch betweeen SIDs')
    return
end

DAHUC=[198830.78,276482.47,739945.75,842544.51,421965.94,105948.37,491924.1,276037.48,259212.63,1349417.31,642212.12,474542.49,597883.22,293568.76,424355.65,367048.65,836516.04,436625.08]; %in order by HUC number 1 to 18, clipped to conus boundary
DAHUC=DAHUC'; %has canada areas
% get the HUC from the first two values in the reach code.
% Extract first two characters
hydrodata.HUC = floor(hydrodata.REACHCODE ./ 1e12);  % adjust divisor to match digits

hucs=unique(hydrodata.HUC);%listed in numerical order
a=~isnan(hucs);
b=find(a==1);
hucs=hucs(b);
numhuc=numel(hucs);

calc7_OutHUC.sedshedabovehuc_km2=zeros(numhuc,numt);
calc7_OutHUC.seddtabovehuc_m3peryr=zeros(numhuc,numt);
calc7_OutHUC.sedabovehuc_km3=zeros(numhuc,numt);
calc7_OutHUC.sedyieldabovehuc_m3perkm2yr=zeros(numhuc,numt);%last year's sed shed influences this years sedimentation
calc7_OutHUC.percentDAtrappedabovehuc=zeros(numhuc,numt);
calc7_OutHUC.numdamsabovehuc_USpermstore=zeros(numhuc,numt);
calc7_OutHUC.capdesignabovehuc=zeros(numhuc,numt);
calc7_OutHUC.capabovehuc=zeros(numhuc,numt);
calc7_OutHUC.sedshedaboveHUCwCAN_km2=zeros(numhuc,numt);


for i=1:numhuc
    idx=find(hydrodata.HUC==hucs(i) & data.PermStorag==1 & data.SID<500000);
    idx2=find(hydrodata.HUC==hucs(i));
    for j=2:numt
        calc7_OutHUC.numdamsabovehuc_USpermstore(i,j)=sum(calc.builtidx(idx)<=j);
        calc7_OutHUC.sedshedabovehuc_km2(i,j)=sum(calc.sedshedALL(idx,j));
        calc7_OutHUC.seddtabovehuc_m3peryr(i,j)=sum(calc.seddt(idx,j));
        calc7_OutHUC.sedabovehuc_km3(i,j)=sum(calc.sed(idx,j))/1e9;
        calc7_OutHUC.sedyieldabovehuc_m3perkm2yr(i,j)=calc7_OutHUC.seddtabovehuc_m3peryr(i,j)/calc7_OutHUC.sedshedabovehuc_km2(i,j-1);%last year's sed shed influences this years sedimentation
        calc7_OutHUC.capdesignabovehuc(i,j)=sum(calc.capdesign(idx,j));
        calc7_OutHUC.capabovehuc(i,j)=sum(calc.cap(idx,j));

        calc7_OutHUC.sedshedaboveHUCwCAN_km2(i,j)=sum(calc.sedshedALL(idx2,j));
        calc7_OutHUC.percentDAtrappedabovehuc(i,j)=calc7_OutHUC.sedshedaboveHUCwCAN_km2(i,j)/(DAHUC(i))*100;
    end
end

calc7_OutHUC.caplossabovehuc=(calc7_OutHUC.capdesignabovehuc-calc7_OutHUC.capabovehuc)./(calc7_OutHUC.capdesignabovehuc).*100;

calc7_OutHUC.SedTotalAboveHUC_2025km3=calc7_OutHUC.sedabovehuc_km3(:,now);
calc7_OutHUC.SedTotalAboveHUC_2050km3=calc7_OutHUC.sedabovehuc_km3(:,now);
calc7_OutHUC.YieldAboveHUC_2025_m3perkm2yr=calc7_OutHUC.sedyieldabovehuc_m3perkm2yr(:,now);
calc7_OutHUC.YieldAboveHUC_2025_AFpermi2yr=calc7_OutHUC.YieldAboveHUC_2025_m3perkm2yr/convert1*convert2;
calc7_OutHUC.CapLossAboveHUC_2025_percent=calc7_OutHUC.caplossabovehuc(:,now);
calc7_OutHUC.HUCid=hucs;
calc7_OutHUC.SeddtAboveHUC_2025_Mm3peryer=calc7_OutHUC.seddtabovehuc_m3peryr(:,now)/1e6;

calc7_OutHUC.DAout=DAout;
calc7_OutHUC.DAHUC=DAHUC;

%%
data5all.HUC=hydrodata.HUC;
save('OutHUC.mat','calc7_OutHUC','data5all','calc5all','-v7.3');

end