%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% SECTIONS OF PI VALUE AS COLOURED PLOT AND TOD as contourlines %%%%%%%%%%%%%%%%%%%%%


% Author: Leonardo Bertini
% Modified: 29th April 2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CODE
depth=load('/Volumes/LaCie_Leonardo/NorESM/all_ramps/Depth_Levels.mat'); %get depth array from one of the files in the folder
depth=depth.depths;
depth=depth';

depth_interp=min(depth):200:max(depth);

%parea=ncread('/Volumes/LaCie_Leonardo/NorESM/all_ramps/grid_gx1v6.nc','parea');
plon=ncread('/Volumes/LaCie_Leonardo/NorESM/all_ramps/grid_gx1v6.nc','plon') ;
plat=ncread('/Volumes/LaCie_Leonardo/NorESM/all_ramps/grid_gx1v6.nc','plat') ;

plat=plat';
plon=plon';
plon=plon;

addpath(genpath('/Volumes/LaCie_Leonardo/NorESM'))
addpath /Volumes/LaCie_Leonardo/NorESM/scripts_jerry
addpath /Volumes/LaCie_Leonardo/NorESM


domain_mask=load ('/Volumes/LaCie_Leonardo/NorESM/Mask_regions/Mask_corrected.mat');
domain_mask=domain_mask.domain_mask;

atl_ind=load('/Volumes/LaCie_Leonardo/NorESM/scripts_jerry/noresm_atl_ind_fixed.asc');

varlist={'ph';'o2';'omegac';'detoc';'templvl';'AOU'}  ;
for varIDX=1:6

if strcmp(varlist{varIDX},'ph')==1
    varname='pH';
    varunit='';
    axmin=7.3;
    axmax=8.3;
    colorscheme=parula(20);
    contourlevels=[20 50 70 100];
    smoothnumbsurface=3;
    smoothnumbsect=1;
    smoothnumbcontour=7;
    
elseif strcmp(varlist{varIDX},'o2')==1
    varname='DO';
    varunit='(\mumol O_2 kg^-^1)';
    axmin=0;
    axmax=350;
    colorscheme=parula(14);
    contourlevels=0:40:280;
    smoothnumbsurface=3;
    smoothnumbsect=1;
    smoothnumbcontour=7;
    
elseif strcmp(varlist{varIDX},'omegac')==1
    varname='\Omega_C';
    varunit='';
    axmin=0.5;
    axmax=1.5;
    colorscheme=flip(redblue,2);
    contourlevels=0:20:280;
    smoothnumbsurface=3;
    smoothnumbsect=1;
    smoothnumbcontour=7;

elseif strcmp(varlist{varIDX},'detoc')==1
    varname='POC flux';
    varunit='(mgC m^-^2 day^-^1)';
    colorscheme=parula(20);
    axmin=0;
    axmax=0.5;
    smoothnumbsurface=3;
    smoothnumbsect=1;
    smoothnumbcontour=7;
    contourlevels=[50 100 150 200 250];

    
elseif strcmp(varlist{varIDX},'templvl')==1
    varname='T';
    varunit='(\circC)';
    colorscheme=parula(25);
    axmin=0;
    axmax=25;
    contourlevels=0:20:280;
    smoothnumbsurface=3;
    smoothnumbsect=1;
    smoothnumbcontour=10;

elseif strcmp(varlist{varIDX},'AOU')==1
    varname='AOU';
    varunit='(\mumol O2 kg^-^1)';
    colorscheme=parula(10);
    axmin=0;
    axmax=200;
   contourlevels=[20 60 100 140 200 280];
    smoothnumbsurface=3;
    smoothnumbsect=1;
    smoothnumbcontour=7;
 
end


%% PI_mean

    PI_mean=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/Pre_industrial/30yrMeans/PI_%s_first_5yr_mean.mat',varlist{varIDX}));
   
    %if O2 convert to umol O2 kg-1
    if varIDX==2
    PI_mean=((PI_mean.PI_mean)/1024)*1000000;
    else
    PI_mean=PI_mean.PI_mean;
    end
    
    PI_mean=permute(PI_mean,[2,1,3]);
    
    %unflipping domain mask
    bit1=domain_mask(:,122:end); %bit1 (200:end,:,:); %reshaping -- position 200 is where Pacific is
    bit2=domain_mask(:,1:121);   %bit2 (1:199,:,:);
    domain_mask_1=cat(2,bit1,bit2);
    
    for l=1:384
    for c=1:320
        for k=1:70
            if domain_mask_1(l,c)~=2 
                PI_mean(l,c,k)=NaN;
            end
        end
    end
    end

%pcolor(PI_mean(:,:,1));shading flat
PI_mean_1=permute(PI_mean,[3,1,2]);
    
    
%% TOD
timevarname='TOD';
fieldTOD=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/all_ramps/filtered/Results/New_time_var_results/%s_results_updated_v1D_poly1_new.mat',varlist{varIDX}),timevarname);
fieldTOD=fieldTOD.TOD;

fieldTOD1=permute(fieldTOD,[3,1,2]); 
%unflipping
bit1=fieldTOD1(:,:,122:end); %bit1 (200:end,:,:); %reshaping -- position 200 is where Pacific is
bit2=fieldTOD1(:,:,1:121);   %bit2 (1:199,:,:);
fieldTOD2=cat(3,bit1,bit2);

%pcolor(squeeze(fieldTOD2(1,:,:)));shading flat

for k=1:70
    for l=1:384
        for c=1:320
            if domain_mask_1(l,c)==2 && ~isnan(PI_mean_1(k,l,c)) && isnan(fieldTOD2(k,l,c))
                fieldTOD2(k,l,c)=-10; %flag of no departure
            end
        end
    end
end



%%
% figure
% pcolor(fieldZERO(:,:,55))
% shading flat
% colorbar
% colormap([[1 1 1]; parula(14) ; [177/255 47/255 245/255]]) %
% caxis([-10 150])
% hold on
% scatter(c,l,'r','filled')
% title(sprintf('%s at line %d column %d depth %d m',varlist{varIDX}, l, c, depths(k)))   


k=[find(depth==500) find(depth==1000) find(depth==2000)];  

hFig = figure(1); clf;
set(hFig, 'Units','centimeters','Position',[0,0,31,23]); %
set(hFig,'visible','off')
m1=0.01; p1=.002; s1=0.005; pb=.02; 

%% surface 1
subaxis(3,2,1,'Spacing',s1,'Padding',p1,'Margin',m1,'PaddingBottom',2*pb,'PaddingRight',0,'PaddingTop',15*p1,'PaddingLeft',0);
m_proj('Lambert','lon',[-100 15], 'lat', [-10 65]) ;
get(gca,'Position');
set(gca,'Position',[-0.08    0.7067-0.01    0.4875    0.2533]);
field3=squeeze(PI_mean_1(k(1),:,:));

for smoothnumber=1:smoothnumbsurface %smoothing 7 times at the specific depth level
    field3=smooth2d(field3);
end


if varIDX==4 %if POC then calculate data in mgC m-2 day-1
m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(field3*12.0107*5*1000)); shading interp
else
m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(field3)); shading interp
end

colormap(colorscheme) % [177/255 47/255 245/255]
caxis([axmin axmax])

m_coast('patch',[.77 .77 .77],'edgecolor','k');

hold on

%Tcontour
fieldTOD3=squeeze(fieldTOD2(k(1),:,:));
for smoothnumber=1:smoothnumbcontour %smoothing 7 times at the specific depth level
    fieldTOD3=smooth2d(round(fieldTOD3));
end
[M,c]=m_contour(shift_atl(plon),shift_atl(plat),shift_atl(fieldTOD3));

c.LevelList=contourlevels; %levels
c.LineColor=[.9 .9 .9];%'grey';
c.LineStyle='-';
clabel(M,c,'FontSize',8) 

%Grid box
%m_grid('backgroundcolor','k','box','fancy','tickdir','out'); %if brown [160/255 82/255 45/255]
m_grid('xtick',-90:30:0,'ytick',0:30:90,'backgroundcolor','k');
set(gca,'layer','top','FontName','Helvetica','fontsize',8,'TickDir','in','TickLength',[.015 .02],'XMinorTick','on','YMinorTick','on');
hold on
m_line(plon(atl_ind),plat(atl_ind),'Color',[102/255 0 0 ],'LineWidth',2);
ttl=title(['a)' sprintf(' %s at %i m',varname,round(depth(k(1))))],'fontsize',14);
set(ttl,'Position',[1.42467418893584e-06,0.90,-4.50359962737050e+15])
if varIDX==3
set(ttl,'Position',[1.42467418893584e-06,0.83,-4.50359962737050e+15])
end


%%  surface 2
subaxis(3,2,3,'Spacing',s1,'Padding',p1,'Margin',m1,'PaddingBottom',2*pb,'PaddingRight',0,'PaddingTop',15*p1,'PaddingLeft',0);
m_proj('Lambert','lon',[-100 15], 'lat', [-10 65]);
set(gca,'Position',[-0.08     0.3783-0.01    0.4875    0.2533]);

field3=squeeze(PI_mean_1(k(2),:,:));

for smoothnumber=1:smoothnumbsurface %smoothing 7 times at the specific depth level
    field3=smooth2d(field3);
end


if varIDX==4 %if POC then calculate data in mgC m-2 day-1
m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(field3*12.0107*5*1000)); shading interp
else
m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(field3)); shading interp
end

colormap(colorscheme) % [177/255 47/255 245/255]
caxis([axmin axmax])

m_coast('patch',[.77 .77 .77],'edgecolor','k');
hold on

%Tcontour
fieldTOD3=squeeze(fieldTOD2(k(2),:,:));
for smoothnumber=1:smoothnumbcontour %smoothing 7 times at the specific depth level
    fieldTOD3=smooth2d(round(fieldTOD3));
end
[M,c]=m_contour(shift_atl(plon),shift_atl(plat),shift_atl(fieldTOD3));

c.LevelList=contourlevels; %levels
c.LineColor=[.9 .9 .9];%'grey';
c.LineStyle='-';
clabel(M,c,'FontSize',8) 

%Grid box
%m_grid('backgroundcolor',[.2 .2 .2],'box','fancy','tickdir','out'); %if brown [160/255 82/255 45/255]
m_grid('xtick',-90:30:0,'ytick',0:30:90,'backgroundcolor','k');
set(gca,'layer','top','FontName','Helvetica','fontsize',8,'TickDir','in','TickLength',[.015 .02],'XMinorTick','on','YMinorTick','on');
hold on
m_line(plon(atl_ind),plat(atl_ind),'Color',[102/255 0 0 ],'LineWidth',2);
ttl=title(['b)' sprintf(' %s at %i m',varname,round(depth(k(2))))],'fontsize',14);
set(ttl,'Position',[1.42467418893584e-06,0.90,-4.50359962737050e+15])
if varIDX==3
set(ttl,'Position',[1.42467418893584e-06,0.83,-4.50359962737050e+15])
end

%% surface 3
subaxis(3,2,5,'Spacing',s1,'Padding',p1,'Margin',m1,'PaddingBottom',2*pb,'PaddingRight',0,'PaddingTop',15*p1,'PaddingLeft',0);
m_proj('Lambert','lon',[-100 15], 'lat', [-10 65]) ;
set(gca,'Position',[-0.08    0.0500-0.01    0.4875    0.2533]);

field3=squeeze(PI_mean_1(k(3),:,:));

for smoothnumber=1:smoothnumbsurface %smoothing 7 times at the specific depth level
    field3=smooth2d(field3);
end


if varIDX==4 %if POC then calculate data in mgC m-2 day-1
m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(field3*12.0107*5*1000)); shading interp
else
m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(field3)); shading interp
end

colormap(colorscheme) % [177/255 47/255 245/255]
caxis([axmin axmax])

m_coast('patch',[.77 .77 .77],'edgecolor','k');
hold on

%Tcontour
fieldTOD3=squeeze(fieldTOD2(k(3),:,:));
for smoothnumber=1:smoothnumbcontour %smoothing 7 times at the specific depth level
    fieldTOD3=smooth2d(round(fieldTOD3));
end
[M,c]=m_contour(shift_atl(plon),shift_atl(plat),shift_atl(fieldTOD3));

c.LevelList=contourlevels; %levels
c.LineColor=[.9 .9 .9];%'grey';
c.LineStyle='-';
clabel(M,c,'FontSize',8) 

%Grid box
m_grid('xtick',-90:30:0,'ytick',0:30:90,'backgroundcolor','k');
set(gca,'layer','top','FontName','Helvetica','fontsize',8,'TickDir','in','TickLength',[.015 .02],'XMinorTick','on','YMinorTick','on');
hold on
m_line(plon(atl_ind),plat(atl_ind),'Color',[102/255 0 0 ],'LineWidth',2)
ttl=title(['c)' sprintf(' %s at %i m',varname,round(depth(k(3))))],'fontsize',14);
set(ttl,'Position',[1.42467418893584e-06,0.90,-4.50359962737050e+15])
if varIDX==3
set(ttl,'Position',[1.42467418893584e-06,0.83,-4.50359962737050e+15])
end


%% section
subplot(3,2,[2,4,6])
get(gca,'Position');
% % %option2
% % subplot(3,2,[2,4,6])
% % get(gca,'Position')
% % set(gca,'Position',[0.5 0.09760 0.4073 0.7914])

bottomflag=-1111;
PI_mean_1(PI_mean_1==bottomflag)=NaN; %this makes the points representing the bottom 'transparent'

%changing all the negative flags to avoid the concentric rings when interpolating
PI_mean_1(PI_mean_1<0)=600;

%getting section data
[sect_fld_original,sect_lat,sect_dep]=makesect_atl_noresm_fixed_z(PI_mean_1,depth,plat,atl_ind); %the function asks for (field, depth, lat or long, file with section indexes)

for smoothnumber=1:smoothnumbsect %smoothing 3 times
    sect_fld_original=smooth2d(sect_fld_original);
end

%section plot
if varIDX==4 %if POC then calculate data in mgC m-2 day-1
pcolor(sect_lat,-sect_dep,sect_fld_original*12.0107*5*1000);
else
pcolor(sect_lat,-sect_dep,sect_fld_original);
end
shading interp
hold on

%TOD
[sect_ctr_original,sect_lat,sect_dep]=makesect_atl_noresm_fixed_z(fieldTOD2,depth,plat,atl_ind); %the function asks for (field, depth, lat or long, file with section indexes)

for smoothnumber=1:smoothnumbcontour %smoothing 7 times
    sect_ctr_original=smooth2d(sect_ctr_original);
end

[M,c]=contour(sect_lat,-sect_dep,sect_ctr_original);

c.LevelList=contourlevels; %levels
c.LineColor=[.9 .9 .9];%'grey';
c.LineStyle='-';
c.LineWidth=2;
clabel(M,c,'FontSize',10) 
hold on
title(['d) North Atlantic section of ' sprintf('%s', varname)],'fontsize',14)

%contourlines
% [sect_fld_original,C,h]=makesect_atl_noresm_fixed_z_COUNTOUR_Tmin(field2,depth,plat,atl_ind);
% h.LevelList=[250 350 450 480];
% h.Color=[.7 .7 .7];
% clabel(C,h,'Color','white','FontSize',8);
% clabel(C,h,'LabelSpacing',300)
% h.LineStyle='-';


colormap(colorscheme) % [177/255 47/255 245/255]
caxis([axmin axmax])
 
h=colorbar;
caxis([axmin axmax])
h.Label.String=sprintf('Pre-industrial mean %s %s',varname,varunit);
h.Location='southoutside';
set(h,'fontsize',12);

xlabel('Latitude','fontsize',12)
ylabel('Depth (m)','fontsize',12)

if varIDX==4 %if detoc then change axis to end at 2000m deep
    axis([-10 65 -1800 0])
else
    axis([-10 65 -6000 0])
end
set(gca, 'color','k','layer','top','FontName','Helvetica','fontsize',12,'TickDir','in','TickLength',[.005 .01],'XMinorTick','on','YMinorTick','on','xtick',-10:10:60,'xticklabel',{'10\circS' '0\circ' '10\circN' '20\circN' '30\circN' '40\circN' '50\circN' '60\circN'  });
set(gca,'Position',[0.40    0.18    0.55    0.77]);

%saving image

set(gcf,'color','w')
hFig.InvertHardcopy='off';
supersizeme(1.2)
print(sprintf('/Users/leonardobertini/OneDrive/IMBRSea/Modules/Master Thesis/Figures_Paper_LEO_NORESM/new_figs/TOD_PImean_%s.png',varlist{varIDX}),'-dpng','-r300')

close all
end









% % %%%%%%%%%%%%ORIGINAL LAYOUT WITHOUT SUBAXIS
% % 
% % %%
% % fig = figure('Units','centimeters','Position',[0,0,21,29.7]); %that's a4 size
% % 
% % %% surface 1
% % subplot(2,3,1)
% % m_proj('Lambert','lon',[-100 20], 'lat', [-10 65]) 
% % 
% % field3=squeeze(PI_mean_1(k(1),:,:));
% % 
% % for smoothnumber=1:smoothnumbsurface %smoothing 7 times at the specific depth level
% %     field3=smooth2d(field3);
% % end
% % 
% % 
% % if varIDX==4 %if POC then calculate data in mgC m-2 day-1
% % m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(field3*12.0107*5*1000)); shading interp
% % else
% % m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(field3)); shading interp
% % end
% % 
% % colormap(colorscheme) % [177/255 47/255 245/255]
% % caxis([axmin axmax])
% % 
% % m_coast('patch',[.7 .7 .7],'edgecolor','k');
% % hold on
% % 
% % %Tcontour
% % fieldTOD3=squeeze(fieldTOD2(k(1),:,:));
% % for smoothnumber=1:smoothnumbcontour %smoothing 7 times at the specific depth level
% %     fieldTOD3=smooth2d(round(fieldTOD3));
% % end
% % [M,c]=m_contour(shift_atl(plon),shift_atl(plat),shift_atl(fieldTOD3));
% % 
% % c.LevelList=contourlevels; %levels
% % c.LineColor=[.7 .7 .7];%'grey';
% % c.LineStyle='-';
% % clabel(M,c,'FontSize',10) 
% % 
% % %Grid box
% % m_grid('backgroundcolor','k','box','fancy','tickdir','out'); %if brown [160/255 82/255 45/255]
% % hold on
% % m_line(plon(atl_ind),plat(atl_ind),'Color',[102/255 0 0 ],'LineWidth',2)
% % title(sprintf('Depth= %.2f m',depth(k(1))))
% % 
% % %%  surface 2
% % subplot(2,3,2)
% % m_proj('Lambert','lon',[-100 20], 'lat', [-10 65]) 
% % 
% % field3=squeeze(PI_mean_1(k(2),:,:));
% % 
% % for smoothnumber=1:smoothnumbsurface %smoothing 7 times at the specific depth level
% %     field3=smooth2d(field3);
% % end
% % 
% % 
% % if varIDX==4 %if POC then calculate data in mgC m-2 day-1
% % m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(field3*12.0107*5*1000)); shading interp
% % else
% % m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(field3)); shading interp
% % end
% % 
% % colormap(colorscheme) % [177/255 47/255 245/255]
% % caxis([axmin axmax])
% % 
% % m_coast('patch',[.7 .7 .7],'edgecolor','k');
% % hold on
% % 
% % %Tcontour
% % fieldTOD3=squeeze(fieldTOD2(k(2),:,:));
% % for smoothnumber=1:smoothnumbcontour %smoothing 7 times at the specific depth level
% %     fieldTOD3=smooth2d(round(fieldTOD3));
% % end
% % [M,c]=m_contour(shift_atl(plon),shift_atl(plat),shift_atl(fieldTOD3));
% % 
% % c.LevelList=contourlevels; %levels
% % c.LineColor=[.7 .7 .7];%'grey';
% % c.LineStyle='-';
% % clabel(M,c,'FontSize',10) 
% % 
% % %Grid box
% % m_grid('backgroundcolor','k','box','fancy','tickdir','out'); %if brown [160/255 82/255 45/255]
% % hold on
% % m_line(plon(atl_ind),plat(atl_ind),'Color',[102/255 0 0 ],'LineWidth',2)
% % title(sprintf('Depth= %.2f m',depth(k(2))))
% % 
% % 
% % %% surface 3
% % subplot(2,3,3)
% % m_proj('Lambert','lon',[-100 20], 'lat', [-10 65]) 
% % 
% % field3=squeeze(PI_mean_1(k(3),:,:));
% % 
% % for smoothnumber=1:smoothnumbsurface %smoothing 7 times at the specific depth level
% %     field3=smooth2d(field3);
% % end
% % 
% % 
% % if varIDX==4 %if POC then calculate data in mgC m-2 day-1
% % m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(field3*12.0107*5*1000)); shading interp
% % else
% % m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(field3)); shading interp
% % end
% % 
% % colormap(colorscheme) % [177/255 47/255 245/255]
% % caxis([axmin axmax])
% % 
% % m_coast('patch',[.7 .7 .7],'edgecolor','k');
% % hold on
% % 
% % %Tcontour
% % fieldTOD3=squeeze(fieldTOD2(k(3),:,:));
% % for smoothnumber=1:smoothnumbcontour %smoothing 7 times at the specific depth level
% %     fieldTOD3=smooth2d(round(fieldTOD3));
% % end
% % [M,c]=m_contour(shift_atl(plon),shift_atl(plat),shift_atl(fieldTOD3));
% % 
% % c.LevelList=contourlevels; %levels
% % c.LineColor=[.7 .7 .7];%'grey';
% % c.LineStyle='-';
% % clabel(M,c,'FontSize',10) 
% % 
% % %Grid box
% % m_grid('backgroundcolor','k','box','fancy','tickdir','out'); %if brown [160/255 82/255 45/255]
% % hold on
% % m_line(plon(atl_ind),plat(atl_ind),'Color',[102/255 0 0 ],'LineWidth',2)
% % title(sprintf('Depth= %.2f m',depth(k(3))))
% % 
% % %% section
% % subplot(2,3,4:6)
% % get(gca,'Position')
% % set(gca,'Position',[0.1300    0.08    0.7750    0.5])
% % 
% % % % %option2
% % % % subplot(3,2,[2,4,6])
% % % % get(gca,'Position')
% % % % set(gca,'Position',[0.5 0.09760 0.4073 0.7914])
% % 
% % % % % % bottomflag=-1111;
% % % % % % PI_mean_1(PI_mean_1==bottomflag)=NaN; %this makes the points representing the bottom 'transparent'
% % % % % 
% % % % % %changing all the negative flags to avoid the concentric rings when interpolating
% % % % % %PI_mean_1(PI_mean_1<0)=600;
% % 
% % %getting section data
% % [sect_fld_original,sect_lat,sect_dep]=makesect_atl_noresm_fixed_z(PI_mean_1,depth,plat,atl_ind); %the function asks for (field, depth, lat or long, file with section indexes)
% % 
% % for smoothnumber=1:smoothnumbsect %smoothing 3 times
% %     sect_fld_original=smooth2d(sect_fld_original);
% % end
% % 
% % %section plot
% % if varIDX==4 %if POC then calculate data in mgC m-2 day-1
% % pcolor(sect_lat,-sect_dep,sect_fld_original*12.0107*5*1000);
% % else
% % pcolor(sect_lat,-sect_dep,sect_fld_original);
% % end
% % shading interp
% % hold on
% % 
% % %TOD contour lines
% % [sect_ctr_original,sect_lat,sect_dep]=makesect_atl_noresm_fixed_z(fieldTOD2,depth,plat,atl_ind); %the function asks for (field, depth, lat or long, file with section indexes)
% % 
% % for smoothnumber=1:smoothnumbcontour %smoothing 7 times
% %     sect_ctr_original=smooth2d(sect_ctr_original);
% % end
% % 
% % [M,c]=contour(sect_lat,-sect_dep,sect_ctr_original)
% % 
% % c.LevelList=contourlevels; %levels
% % c.LineColor=[.7 .7 .7];%'grey';
% % c.LineStyle='-';
% % clabel(M,c,'FontSize',10) 
% % hold on
% % 
% % colormap(colorscheme) % [177/255 47/255 245/255]
% % caxis([axmin axmax])
% %  
% % h=colorbar;
% % caxis([axmin axmax])
% % h.Location='northoutside';
% % h.Label.String=sprintf('Pre-industrial mean %s %s',varname,varunit);
% % h.Location='northoutside';
% % 
% % xlabel('Latitude (\circ)')
% % ylabel('Depth (m)')
% % if varIDX==4 %if detoc then change axis to end at 2000m deep
% %     axis([-10 65 -2000 0])
% % else
% %     axis([-10 65 -6000 0])
% % end
% % set(gca, 'color','k' );
% % 
% % 
% % % %%DRAFT TO MAKE GIF FILES%%%%%%%%
% % % 
% % % 
% % % for k = 1 : (length(depth)-1)
% % % %     [A,map] = rgb2ind(im{k},256);
% % % %     if k == 1
% % % %         imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
% % % %     else
% % % %         imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
% % % %     end
% % % %  end
% % % % 
% % % %  close all
% % % %  
% % %  
% % %  
% % % %% section pcolor
% % % 
% % % fig=figure('Units','centimeters','Position',[0,0,29.7,21]); %that's a4 size
% % % [sect_fld_original]=makesect_atl_noresm_fixed_z(field2,depth,plat,atl_ind); %the function asks for (field, depth, lat or long, file with section indexes)
% % % hold on
% % % 
% % % %% Tmin line contour
% % % fieldTmin=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/all_ramps/filtered/%s_results_updated_v1D_polyfit.mat',varlist{varIDX}),'TpeakMAX');
% % % fieldTmin=fieldTmin.TpeakMAX;
% % % 
% % % fieldTmin=permute(fieldTmin,[3,1,2]);
% % % 
% % % %unflipping
% % % bit1=fieldTmin(:,:,120:end); %bit1 (200:end,:,:); %reshaping -- position 200 is where Pacific is
% % % bit2=fieldTmin(:,:,1:119);   %bit2 (1:199,:,:);
% % % fieldTmin=cat(3,bit1,bit2);
% % % 
% % % [sect_fld_Tmin,C,h]=makesect_atl_noresm_fixed_z_COUNTOUR_Tmin(fieldTmin,depth,plat,atl_ind); %the function asks for (field, depth, lat or long, file with section indexes)
% % % h.LineColor='white';
% % % h.LevelList=100:50:500;
% % % clabel(C,h,'color','white')
% % % %%
% % % 
% % % 
% % % shading interp
% % % 
% % % if varIDX~=2
% % % colordummy=parula(10);
% % % colordummy=colordummy(1:9,:); %this removes bright yellow
% % % colormap(colordummy)
% % % caxis([0 90])
% % % h=colorbar;
% % % h.Location='northoutside';
% % % h.Ticks=0:10:90;
% % % h.Label.String=sprintf('%s ToD (yr)',varname);
% % % %h.TickLabels(1)={'No departure'};
% % % h.TickLabels(end)={'\geq90 yr'};
% % % 
% % % else
% % % colordummy=parula(15);
% % % colordummy=colordummy(1:13,:); %this removes bright yellow
% % % colormap(colordummy)
% % % caxis([0 130])
% % % h=colorbar;
% % % h.Location='northoutside';
% % % h.Ticks=0:10:130;
% % % h.Label.String=sprintf('%s ToD (yr)',varname);
% % % %h.TickLabels(1)={'No departure'};
% % % h.TickLabels(end)={'\geq130 yr'};
% % % end
% % % 
% % % 
% % % xlabel('Latitude (\circ)')
% % % ylabel('Depth (m)')
% % % axis([-10 75 -6000 0])
% % % set(gca, 'color', 'k');
% % % hold on
% % %yline(-depth(k),'--','LineWidth',1)
% % %print(sprintf('/Volumes/LaCie_Leonardo/NorESM/Figures/Section_North_Atlantic_%s_%s.png',timevarname,varlist{varIDX}),'-dpng','-r300')

    
