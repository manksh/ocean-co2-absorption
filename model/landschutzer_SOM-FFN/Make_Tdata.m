
%====================== add path to some tools ============================

addpath /Users/peter/Documents/MATLAB/m_map   % mapping tool for plots
addpath /Users/peter/Documents/MATLAB/mexcdf/mexnc % to read file info and variables in matlab
addpath /Users/peter/Documents/MATLAB/mexcdf/snctools % idem dito
addpath /Users/peter/Documents/MATLAB/nansuite % mean ignoring nan
addpath /Users/peter/Documents/MATLAB/functions
addpath /Users/peter/Documents/MPI/MPI-SOM-FFN/v2017/functions
addpath /Users/peter/Documents/MATLAB
% -------------------------------------------------------------------------
% 1. load the data which are in (time,lat,lon) format
% -------------------------------------------------------------------------

% V2 used: ECCO, Raynolds, Globalviev, Taka, Globecolor

chl=nc_varget('training_data/Chl_2D_mon_CESM002_1x1_198201-201701.nc', 'Chl');
mld=nc_varget('training_data/MLD_2D_mon_CESM002_1x1_198201-201701.nc','MLD');
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

years=1982:1/12:2017;

latt=zeros(420,180,360);
lonn=zeros(420,180,360);
months=zeros(420,180,360);
time=zeros(420,180,360);

for j=1:180
    latt(1:420,j,1:360)=lat(j);
end

for k=1:360
    lonn(1:420,1:180,k)=lon(k);
end

for uu=1:12
    months(uu:12:end,1:end,1:end)=uu;
end

for oo=1:420
    
    time(oo,:,:)=years(oo);
    
end

years=floor(years);


for year2go=1982:2016

    year2go
    
% -------------------------------------------------------------------------
% 3. The first parameter is SST
% -------------------------------------------------------------------------

y_ind=find(years==year2go);
sst1=sst(y_ind,:,:);
time1=time(y_ind,:,:);
lat1=latt(y_ind,:,:);
lon1=lonn(y_ind,:,:);
months1=months(y_ind,:,:);

%=============== sst is 3 dim array - reshape to 1D =======================

sst_new = reshape(sst1, prod(size(sst1)), 1);
time_new = reshape(time1, prod(size(time1)), 1);
lat_new = reshape(lat1, prod(size(lat1)), 1);
lon_new = reshape(lon1, prod(size(lon1)), 1);
months_new = reshape(months1, prod(size(months1)), 1);

%============== create an array with the same length for year =============

    year2go_new=zeros(length(sst_new),1);
    year2go_new(:)=year2go;
    
    %year2go_new=year2go_new';

%=================== put data in the same array ===========================
% 1=year; 2=months; 3=matlab time 4=lat; 5=lon; 6 =sst

sst_NA_ALL=[year2go_new months_new time_new lat_new lon_new sst_new];


clear time1 i ;

%==== Insert columns at the end of the existing data for parameter: =======
% 7=mld 8=chl 9=sss 10=atm_co2 11=data_taka 12=sst anom
% 13=mld anom 14 =chl anom 15 =sss anom 16=atm_co2 anom 17=NaN 18=biomes

sst_NA_ALL = sst_NA_ALL(:,[1 2 3 4 5 6 6 6 6 6 6 6 6 6 6 6 6 6]); 

sst_NA_ALL(:,7:18)=NaN;

%--------------------------------------------------------------------------
% 4. Coinciding training data
%--------------------------------------------------------------------------

mld1=mld(y_ind,:,:);
atm_co21=atm_co2(y_ind,:,:);
chl1=chl(y_ind,:,:);
sss1=sss(y_ind,:,:);
%data_taka1=data_taka(y_ind,:,:);
sst_anom1=sst_anom(y_ind,:,:);
mld_anom1=mld_anom(y_ind,:,:);
chl_anom1=chl_anom(y_ind,:,:);
sss_anom1=sss_anom(y_ind,:,:);
atm_anom1=atm_anom(y_ind,:,:);
biomes1=biomes(y_ind,:,:);
% 

mld_new = reshape(mld1, prod(size(mld1)), 1);
atm_co2_new = reshape(atm_co21, prod(size(atm_co21)), 1);
chl_new = reshape(chl1, prod(size(chl1)), 1);
sss_new = reshape(sss1, prod(size(sss1)), 1);
%data_taka_new = reshape(data_taka1, prod(size(data_taka1)), 1);
mld_anom_new = reshape(mld_anom1, prod(size(mld_anom1)), 1);
chl_anom_new = reshape(chl_anom1, prod(size(chl_anom1)), 1);
sst_anom_new = reshape(sst_anom1, prod(size(sst_anom1)), 1);
sss_anom_new = reshape(sss_anom1, prod(size(sss_anom1)), 1);
atm_anom_new = reshape(atm_anom1, prod(size(atm_anom1)), 1);
biomes_new = reshape(biomes1, prod(size(biomes1)), 1);


sst_NA_ALL(:,7)=mld_new;
sst_NA_ALL(:,8)=chl_new;
sst_NA_ALL(:,9)=sss_new;
sst_NA_ALL(:,10)=atm_co2_new;
%sst_NA_ALL(:,11)=data_taka_new;
sst_NA_ALL(:,12)=sst_anom_new;
sst_NA_ALL(:,13)=mld_anom_new;
sst_NA_ALL(:,14)=chl_anom_new;
sst_NA_ALL(:,15)=sss_anom_new;
sst_NA_ALL(:,16)=atm_anom_new;
sst_NA_ALL(:,18)=biomes_new;

clear a b c d i j k lat_new lon_new time_new sst_new year2go_new
clear bathy_new mld_new chl_new data_taka_new sss_new spec_rate_new atm_co2_new ssh_new curr_speed_new wind_new

%======== order the data (as L), first in time, then in lon and lat =======

Tdata = sortrows(sst_NA_ALL,[3 4 5]);

% 1=year 2=month 3=lat 4=lon 5=SST 6=mld 7=chl 8=sss 9=atm_co2 10=data_taka 11=sst anom
% 12=mld anom 13 =chl anom 14 =sss anom 15=atm_co2 anom 16=NaN 17=biomes

Tdata=Tdata(:,[1 2 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18]);

% -------------------------------------------------------------------------
% 6. save output
% -------------------------------------------------------------------------

%eval(['Tdata_' num2str(year2go) '_500m = Tdata_500m;']);
eval(['Tdata_' num2str(year2go) ' = Tdata;']);

saveloc=['Tdata/Tdata_' num2str(year2go) '.mat'];
savedat=['Tdata_' num2str(year2go)];

save(saveloc, savedat);

clear saveloc savedat Tdata_500m

end
% -------------------------------------------------------------------------
% 7. graphical check if ok
% -------------------------------------------------------------------------

xx=size(lat);
yy=size(lon);

s = ['Tdata_' num2str(year2go)];
plot_array=eval(s);
clear s
indxx=find(plot_array(:,2)==2 & plot_array(:,1)==year2go);

LON=plot_array(indxx,4);
LAT=plot_array(indxx,3);

[field]=vec_to_array1(plot_array(indxx,5:17),lon,lat,LON,LAT);

map_plot_2016(lon, lat, squeeze(field(1,:,:)), 'equidistant', -180, 180, -90, 90, 1, 0, 30,2, 'SST',1,2)
map_plot_2016(lon, lat, squeeze(field(2,:,:)), 'equidistant', -180, 180, -90, 90, 3, 0, 8,1, 'MLD',1,2)
map_plot_2016(lon, lat, squeeze(field(3,:,:)), 'equidistant', -180, 180, -90, 90, 4, -4, 4,1, 'CHL',1,2)