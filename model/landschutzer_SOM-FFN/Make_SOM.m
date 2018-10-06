%====================== add path to some tools ============================

addpath /Users/peter/Documents/MATLAB/m_map   % mapping tool for plots
addpath /Users/peter/Documents/MATLAB/mexcdf/mexnc % to read file info and variables in matlab
addpath /Users/peter/Documents/MATLAB/mexcdf/snctools % idem dito
addpath /Users/peter/Documents/MATLAB/nansuite % mean ignoring nan
addpath /Users/peter/Documents/MATLAB/functions
addpath /Users/peter/Documents/MPI/MPI-SOM-FFN/v2017/functions
addpath /Users/peter/Documents/MATLAB
%--------------------------------------------------------------------------
% 1) load cluster data 
%--------------------------------------------------------------------------

mld=nc_varget('training_data/MLD_2D_mon_CESM002_1x1_198201-201701.nc','MLD');
pco2=nc_varget('training_data/pCO2_2D_mon_CESM002_1x1_198201-201701.nc','pCO2');
sss=nc_varget('training_data/SSS_2D_mon_CESM002_1x1_198201-201701.nc','SSS');
sst=nc_varget('training_data/SST_2D_mon_CESM002_1x1_198201-201701.nc','SST');
chl=nc_varget('training_data/Chl_2D_mon_CESM002_1x1_198201-201701.nc', 'Chl');

lon=nc_varget('training_data/MLD_2D_mon_CESM002_1x1_198201-201701.nc','xlon');
lat=nc_varget('training_data/MLD_2D_mon_CESM002_1x1_198201-201701.nc','ylat');
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
mld=mld_tmp;
clear mld_tmp

pco2_tmp=pco2;
pco2_tmp(:,:,1:end/2)=pco2(:,:,end/2+1:end);
pco2_tmp(:,:,end/2+1:end)=pco2(:,:,1:end/2);
clear pco2
pco2=pco2_tmp;
clear pco2_tmp

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
% 2) lat-lon
%--------------------------------------------------------------------------


latt=zeros(420,180,360);
lonn=zeros(420,180,360);

for j=1:180
    latt(1:420,j,1:360)=lat(j);
end

for k=1:360
    lonn(1:420,1:180,k)=lon(k);
end

%--------------------------------------------------------------------------
% 2.5) take average
%--------------------------------------------------------------------------

timevec=1982:1/12:2017-1/12;
ind=find(timevec>=1982 &timevec<2017);

sst=sst(ind,:,:);
mld=mld(ind,:,:);
sss=sss(ind,:,:);
data_taka=pco2(ind,:,:);


for uu=1:12
   sst_an(uu,:,:)=nanmean(sst(uu:12:end,:,:));
   mld_an(uu,:,:)=nanmean(mld(uu:12:end,:,:));
   sss_an(uu,:,:)=nanmean(sss(uu:12:end,:,:));
   taka_an(uu,:,:)=nanmean(data_taka(uu:12:end,:,:));
end

lonsst_an=lonn(1:12,:,:);
latsst_an=latt(1:12,:,:);
monthx=zeros(size(latsst_an));

for ww=1:12
    monthx(ww,:,:)=ww;
end

%--------------------------------------------------------------------------
% 3) reshape and rearrange for SOM
%--------------------------------------------------------------------------


sst_new = reshape(sst_an, prod(size(sst_an)), 1);
sss_new = reshape(sss_an, prod(size(sss_an)), 1);
mld_new = reshape(mld_an, prod(size(mld_an)), 1);
taka_new = reshape(taka_an, prod(size(taka_an)), 1);
lat_new = reshape(latsst_an, prod(size(latsst_an)), 1);
lon_new = reshape(lonsst_an, prod(size(lonsst_an)), 1);
month_new = reshape(monthx, prod(size(monthx)), 1);

clear lon lat month

temp2=[sst_new mld_new sss_new taka_new lat_new lon_new month_new];


[index2,Train_SOM]=nanremoveV3_SOM(temp2);

LAT2=lat_new(index2);
LON2=lon_new(index2);
month=month_new(index2);
%Label_SOM=fCO2;

%Train_SOM1=mapminmax(Train_SOM');
%clear Train_SOM
%Train_SOM=Train_SOM1';

clear index1 index2 temp1 temp2 sum_index sum_index1

clear idx1 idx2
%--------------------------------------------------------------------------
% 5) SOM part to identify biomes
%--------------------------------------------------------------------------

    
net = selforgmap([4 4],100,3,'hextop','dist');
net.trainParam.epochs = 200;
[net,tr] = train(net,Train_SOM');
    
y = net(Train_SOM');
classes = vec2ind(y);    

for month2plot=1:12
 
            index4=find(month==month2plot);

            classesYYY=classes(index4);

                lat =-89.5:1:89.5;
                lon =-179.5:1:179.5;
                LATY = LAT2(index4);
                LONY = LON2(index4);
                [classesY]=vec_to_array2(classesYYY,lon,lat,LONY,LATY);

                biomeY=classesY;
                biomeY2(month2plot,:,:)=biomeY(:,:);
end

lon=-179.5:1:179.5;
lat=-89.5:1:89.5;

%--------------------------------------------------------------------------
% 5.5) arctic 1 biome
%--------------------------------------------------------------------------

[lonn latt]=meshgrid(lon,lat);
%crop Arctic everywhere else

ind2=find(lonn>30 & latt>65);
ind3=find(lonn<-90 & latt>65);
ind4=find(latt>80);

indtot2=[ind2; ind3; ind4];
indtot=unique(indtot2);

biomeset=biomeY2(:,indtot);
biomeset(isnan(biomeset))=[];
mode_b=mode(biomeset);


biomeY2(:,indtot)=mode_b;

%--------------------------------------------------------------------------
% no smoothing
%--------------------------------------------------------------------------

    new3=biomeY2;

for ii=1:12:420
    biomes(ii:ii+11,:,:)=new3(:,:,:); 
end

%--------------------------------------------------------------------------
% 7) save and plot 3-D biomes
%--------------------------------------------------------------------------


%array_test=squeeze(biomes(2,:,:));
array_test=squeeze(mode(biomes));

saveloc=['training_data/biomes.mat'];
save(saveloc,'biomes');


figure(6);clf; 
m_proj('equidistant','lon',[-180 180],'lat',[-90 90]);
h=m_pcolor(lon,lat,array_test);
set(h,'edgecolor','none');
colorbar('hori');caxis([1 16]);
%colormap(redblue);
m_coast('patch','w');
m_grid('xtick',[ -135 -90 -45 0 45 90 135], ...
    'ytick',[-90 -60 -30 0 30 60 90],'fontsize',10);



