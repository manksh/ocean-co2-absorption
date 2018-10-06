%====================== add path to some tools ============================

addpath /Users/peter/Documents/MATLAB/m_map   % mapping tool for plots
addpath /Users/peter/Documents/MATLAB/mexcdf/mexnc % to read file info and variables in matlab
addpath /Users/peter/Documents/MATLAB/mexcdf/snctools % idem dito
addpath /Users/peter/Documents/MATLAB/nansuite % mean ignoring nan
addpath /Users/peter/Documents/MATLAB/functions
addpath /Users/peter/Documents/MPI/MPI-SOM-FFN/v2017/functions
addpath /Users/peter/Documents/MATLAB
%==========================================================================

chl=nc_varget('training_data/Chl_2D_mon_CESM002_1x1_198201-201701.nc', 'Chl');
mld=nc_varget('training_data/MLD_2D_mon_CESM002_1x1_198201-201701.nc','MLD');
pco2=nc_varget('training_data/pCO2_2D_mon_CESM002_1x1_198201-201701.nc','pCO2_socat');
pco2(pco2==0)=NaN;
sss=nc_varget('training_data/SSS_2D_mon_CESM002_1x1_198201-201701.nc','SSS');
sst=nc_varget('training_data/SST_2D_mon_CESM002_1x1_198201-201701.nc','SST');
xco22=nc_varget('training_data/XCO2_1D_mon_CESM002_native_198201-201701.nc','XCO2');
load('training_data/biomes.mat');

lon=nc_varget('training_data/MLD_2D_mon_CESM002_1x1_198201-201701.nc','xlon');
lat=nc_varget('training_data/MLD_2D_mon_CESM002_1x1_198201-201701.nc','ylat');

for tt=1:421
    
    atm_co2(tt,1:180,1:360)=xco22(tt);
end

%--------------------------------------------------------------------------
% 1.5) change dimensions
%--------------------------------------------------------------------------

lon_tmp=lon;
lon_tmp(1:end/2)=lon(end/2+1:end)-360;
lon_tmp(end/2+1:end)=lon(1:end/2);
clear lon
lon=lon_tmp;
clear lon_tmp

pco2_tmp=pco2;
pco2_tmp(:,:,1:end/2)=pco2(:,:,end/2+1:end);
pco2_tmp(:,:,end/2+1:end)=pco2(:,:,1:end/2);
clear pco2
pco2=pco2_tmp;
clear pco2_tmp

mld_tmp=mld;
mld_tmp(:,:,1:end/2)=mld(:,:,end/2+1:end);
mld_tmp(:,:,end/2+1:end)=mld(:,:,1:end/2);
clear mld
mld=log(mld_tmp);
clear mld_tmp

chl_tmp=chl;
chl_tmp(:,:,1:end/2)=chl(:,:,end/2+1:end);
chl_tmp(:,:,end/2+1:end)=chl(:,:,1:end/2);
clear chl
chl=log(chl_tmp);
chl=real(chl);
clear chl_tmp

sss_tmp=sss;
sss_tmp(:,:,1:end/2)=sss(:,:,end/2+1:end);
sss_tmp(:,:,end/2+1:end)=sss(:,:,1:end/2);
clear sss
sss=sss_tmp;
clear sss_tmp

sst_tmp=sst;
sst_tmp(:,:,1:end/2)=sst(:,:,end/2+1:end);
sst_tmp(:,:,end/2+1:end)=sst(:,:,1:end/2);
clear sst
sst=sst_tmp;
clear sst_tmp

%--------------------------------------------------------------------------
% 1.7) anomalies
%--------------------------------------------------------------------------

for uu=1:12
   sst_mean(uu,:,:)=nanmean(sst(uu:12:end-1,:,:)); 
   mld_mean(uu,:,:)=nanmean(mld(uu:12:end-1,:,:));
   chl_mean(uu,:,:)=nanmean(chl(uu:12:end-1,:,:));
   sss_mean(uu,:,:)=nanmean(sss(uu:12:end-1,:,:));
   atm_mean(uu,:,:)=nanmean(atm_co2(uu:12:end-1,:,:));
end

for uu=1:12:420
    sst_anom(uu:uu+11,:,:)=sst(uu:uu+11,:,:)-sst_mean;
    sss_anom(uu:uu+11,:,:)=sss(uu:uu+11,:,:)-sss_mean;
    chl_anom(uu:uu+11,:,:)=chl(uu:uu+11,:,:)-chl_mean;
    mld_anom(uu:uu+11,:,:)=mld(uu:uu+11,:,:)-mld_mean;
    atm_anom(uu:uu+11,:,:)=atm_co2(uu:uu+11,:,:)-atm_mean;
end

clear sst_mean sss_mean chl_mean mld_mean atm_mean

%--------------------------------------------------------------------------
% 2) lat-lon
%--------------------------------------------------------------------------


latt=zeros(420,180,360);
lonn=zeros(420,180,360);
months=zeros(420,180,360);

for j=1:180
    latt(1:420,j,1:360)=lat(j);
end

for k=1:360
    lonn(1:420,1:180,k)=lon(k);
end

for uu=1:12
    months(uu:12:end,1:end,1:end)=uu;
end

years=1982:1/12:2017;
years=floor(years);

% -------------------------------------------------------------------------
% 2. load observation files
% -------------------------------------------------------------------------

for year2go=1982:2016

    year2go
    
y_ind=find(years==year2go);
co21=pco2(y_ind,:,:);
lat1=latt(y_ind,:,:);
lon1=lonn(y_ind,:,:);
months1=months(y_ind,:,:);

%=============== sst is 3 dim array - reshape to 1D =======================

co2_new = reshape(co21, prod(size(co21)), 1);
lat_new = reshape(lat1, prod(size(lat1)), 1);
lon_new = reshape(lon1, prod(size(lon1)), 1);
months_new = reshape(months1, prod(size(months1)), 1);
    
year2go_new=zeros(length(co2_new),1);
year2go_new(:)=year2go;

% columns: 1.year/2.month/3.lat/4.lon/5.fco2

Ldata_binned=[year2go_new months_new lat_new lon_new co2_new];
%Ldata_binned1=double(Ldata_binned1);
%Ldata_binned1(:,3)=datenum(Ldata_binned1(:,1),Ldata_binned1(:,2),1);
%Ldata_binned=unique(Ldata_binned1, 'rows');
%clear temp Ldata_binned1;
% column 6 in Ldata_All is fco2
%ind=find(Ldata_binned(:,5)<50);
%Ldata_binned(ind,:)=[];

Ldata_binned(:,6)=NaN;
Ldata_coincided = Ldata_binned(:,[1 2 3 4 6 6 6 6 6 6 6 6 6 6 6 5 6]);

% -------------------------------------------------------------------------
% 3. Coinciding other training parameter
% -------------------------------------------------------------------------

mld1=mld(y_ind,:,:);
atm_co21=atm_co2(y_ind,:,:);
chl1=chl(y_ind,:,:);
sss1=sss(y_ind,:,:);
sst1=sst(y_ind,:,:);
%data_taka1=data_taka(y_ind,:,:);
sst_anom1=sst_anom(y_ind,:,:);
mld_anom1=mld_anom(y_ind,:,:);
chl_anom1=chl_anom(y_ind,:,:);
sss_anom1=sss_anom(y_ind,:,:);
atm_anom1=atm_anom(y_ind,:,:);
biomes1=biomes(y_ind,:,:);

mld_new = reshape(mld1, prod(size(mld1)), 1);
atm_co2_new = reshape(atm_co21, prod(size(atm_co21)), 1);
chl_new = reshape(chl1, prod(size(chl1)), 1);
sss_new = reshape(sss1, prod(size(sss1)), 1);
sst_new = reshape(sst1, prod(size(sst1)), 1);
%data_taka_new = reshape(data_taka1, prod(size(data_taka1)), 1);
chl_anom_new = reshape(chl_anom1, prod(size(chl_anom1)), 1);
mld_anom_new = reshape(mld_anom1, prod(size(chl_anom1)), 1);
sst_anom_new = reshape(sst_anom1, prod(size(sst_anom1)), 1);
sss_anom_new = reshape(sss_anom1, prod(size(sss_anom1)), 1);
atm_anom_new = reshape(atm_anom1, prod(size(atm_anom1)), 1);
biomes_new = reshape(biomes1, prod(size(biomes1)), 1);

% 1=year 2=month 3=lat 4=lon 5=SST 6=mld 7=chl 8=sss 9=atm_co2 10=data_taka 11=sst anom
% 12=mld anom 13 =chl anom 14 =sss anom 15=atm_co2 anom 16=obs 17=biomes

Ldata_coincided(:,5)=sst_new;
Ldata_coincided(:,6)=mld_new;
Ldata_coincided(:,7)=chl_new;
Ldata_coincided(:,8)=sss_new;
Ldata_coincided(:,9)=atm_co2_new;
%Ldata_coincided(:,10)=data_taka_new;
Ldata_coincided(:,11)=sst_anom_new;
Ldata_coincided(:,13)=chl_anom_new;
Ldata_coincided(:,12)=mld_anom_new;
Ldata_coincided(:,14)=sss_anom_new;
Ldata_coincided(:,15)=atm_anom_new;
Ldata_coincided(:,17)=biomes_new;

%======================= remove co2 nans ==================================

ind2= isnan(Ldata_coincided(:,16));
Ldata_coincided(ind2,:)=[];

Ldata_coincided = sortrows(Ldata_coincided,[1 2 3 4]);
 
clear ind2;

% -------------------------------------------------------------------------
% 4. save output
% -------------------------------------------------------------------------
eval(['Ldata_' num2str(year2go) ' =  Ldata_coincided;']);

saveloc=['Ldata/Ldata_' num2str(year2go)];

savedat=['Ldata_' num2str(year2go)];

save(saveloc, savedat);

clear saveloc savedat

end
% -------------------------------------------------------------------------
% 5. graphical check if ok
% -------------------------------------------------------------------------

lat =-89.5:1:89.5;
lon =-179.5:1:179.5;

xx=size(lat);
yy=size(lon);

month2plot=4;

indxx=find(Ldata_coincided(:,2)==month2plot & Ldata_coincided(:,1)==year2go);

LON=Ldata_coincided(indxx,4);
LAT=Ldata_coincided(indxx,3);

[field]=vec_to_array1(Ldata_coincided(indxx,5:17),lon,lat,LON,LAT);


map_plot_new(lon, lat, squeeze(field(1,:,:)), 'equidistant', -180, 180, -90, 90, 1, 0, 30, 'SST',1,2)
map_plot_new(lon, lat, squeeze(field(2,:,:)), 'equidistant', -180, 180, -90, 90, 3, 0, 8, 'MLD',1,2)
map_plot_new(lon, lat, squeeze(field(12,:,:)), 'equidistant', -180, 180, -90, 90, 19, 280, 440, 'biomes',1,2)