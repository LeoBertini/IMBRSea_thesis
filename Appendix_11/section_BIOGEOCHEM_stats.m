%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SECTIONS STATISTICAL ANALYSES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%PRODUCES SURFACE AND SECTION PLOTS FOR THE RESULTS OF THE STATISTICAL
%ANALYSES

% Author: Leonardo Bertini
% Modified: 29th April 2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CODE

%% Loading Initial Variables

depth=load('/Volumes/LaCie_Leonardo/NorESM/all_ramps/Depth_Levels.mat'); %get depth array from one of the files in the folder
depth=depth.depths;
depth=depth';

depth_interp=min(depth):100:max(depth); %this is the depth array corresponding to the evenly spaced depths

%parea=ncread('/Volumes/LaCie_Leonardo/NorESM/all_ramps/grid_gx1v6.nc','parea');
plon=ncread('/Volumes/LaCie_Leonardo/NorESM/all_ramps/grid_gx1v6.nc','plon') ;
plat=ncread('/Volumes/LaCie_Leonardo/NorESM/all_ramps/grid_gx1v6.nc','plat') ;

plat=plat';
plon=plon';
plon=plon;

addpath(genpath('/Volumes/LaCie_Leonardo/NorESM/'))


atl_ind=load('/Volumes/LaCie_Leonardo/NorESM/scripts_jerry/noresm_atl_ind_fixed.asc');

pairs={'pair1';'pair2';'pair3'};
deltaphmode='off';
deltamode='off';
surfacemode='off';
%PAIR 1 - End of mitigation
%PAIR 2 - Middle of extension
%PAIR 3 - End of extension
pairIDX=3;


if pairIDX==1
    period='Pre-Industrial x End of Mitigation phase';
elseif pairIDX==2
    period='Pre-Industrial x Middle of Extension phase';
elseif pairIDX==3
    period='Pre-Industrial x End of Extension phase';
end

varlist={'ph';'o2';'omegac';'detoc';'templvl';'AOU';'o2sat'}  ;
for varIDX=[1]

%% Specifying parameters for each variable

if strcmp(varlist{varIDX},'ph')==1
    varname='[H^+]';%varname='pH';
    varunit='';
    axmin=-25;
    axmax=+25;
    colorscheme=parula(20);
    contourlevels=280:20:480;
    smoothnumbsurface=3;
    smoothnumbsect=3;
    smoothnumbcontour=7;
    
elseif strcmp(varlist{varIDX},'o2')==1
    varname='DO';
    varunit='(mol O2 m^-^3)';
    axmin=-15;
    axmax=15;
    colorscheme=parula(20);
    contourlevels=280:20:480;
    smoothnumbsurface=3;
    smoothnumbsect=3;
    smoothnumbcontour=7;
    
elseif strcmp(varlist{varIDX},'omegac')==1
    varname='\Omega_C';
    varunit='';
    axmin=-25;
    axmax=25;
    colorscheme=flip(redblue,2);
    contourlevels=280:20:480;
    smoothnumbsurface=3;
    smoothnumbsect=3;
    smoothnumbcontour=7;

elseif strcmp(varlist{varIDX},'detoc')==1
    varname='POC flux';
    varunit='(mgC m^-^2 day^-^1)';
    colorscheme=parula(20);
    axmin=0;
    axmax=0.5;
     smoothnumbsurface=3;
    smoothnumbsect=3;
    smoothnumbcontour=7;
    contourlevels=280:20:480;

    
elseif strcmp(varlist{varIDX},'templvl')==1
    varname='T';
    varunit='(\circC)';
    colorscheme=parula(25);
    axmin=-30;
    axmax=30;
    contourlevels=280:20:480;
    smoothnumbsurface=3;
    smoothnumbsect=1;
    smoothnumbcontour=10;
    
    if strcmp(surfacemode,'on') 
    axmin=-15;
    axmax=15;

    end

elseif strcmp(varlist{varIDX},'AOU')==1
    varname='AOU';
    varunit='(\mumol O2 kg^-^1)';
    colorscheme=parula(20);
    axmin=-25;
    axmax=25;
    contourlevels=280:20:480;
    smoothnumbsurface=3;
    smoothnumbsect=1;
    smoothnumbcontour=7;
    
    elseif strcmp(varlist{varIDX},'o2sat')==1
    varname='O_2 saturation';
    varunit='(\mumol O2 kg^-^1)';
    colorscheme=parula(20);
    axmin=-25;
    axmax=25;
    contourlevels=280:20:480;
    smoothnumbsurface=3;
    smoothnumbsect=1;
    smoothnumbcontour=7;
   
end

%%

%loading everything
%load /Volumes/LaCie_Leonardo/NorESM/StatsResults/ph_stats_updated_v1D_vertinterp.mat

%% Loading Percentage change data for a specific period

if strcmp(varlist{varIDX},'ph')==1 %If it's pH, data is shown in terms of hydrogen ion percentage change
    %Percentage_change_hidrogen_ion
    field=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/StatsResults/%s_stats_updated_v1D_vertinterp.mat',varlist{varIDX})...
        ,sprintf('Percentage_change_hidrogen_ion_%s',pairs{pairIDX}));
    field=eval(sprintf('field.Percentage_change_hidrogen_ion_%s',pairs{pairIDX}));
    
else
    %Data is shown in terms of Percentage_change_mean
    field=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/StatsResults/%s_stats_updated_v1D_vertinterp.mat',varlist{varIDX})...
        ,sprintf('Percentage_change_mean_%s',pairs{pairIDX}));
    field=eval(sprintf('field.Percentage_change_mean_%s',pairs{pairIDX}));
end

fieldZERO=field;

%=======DELTA PH ============
if strcmp(varlist{varIDX},'ph') && strcmp(deltaphmode,'on') %then transform percentage change into Delta pH
    
    for l=1:size(fieldZERO,1)
        for c=1:size(fieldZERO,2)
            for k=1:size(fieldZERO,3)
                
                if fieldZERO(l,c,k)<0
                    delta(l,c,k)=2-log10(fieldZERO(l,c,k)+100);
                elseif fieldZERO(l,c,k)>0
                    delta(l,c,k)=-2+log10(fieldZERO(l,c,k)-100);
                elseif isnan(fieldZERO(l,c,k))
                    delta(l,c,k)=NaN;
                end
                
            end
        end
    end
    delta=real(delta);
    fieldZERO=delta;
    varname='\DeltapH';
    varunit='';
    
    if pairIDX==1
    axmin=-1;
    axmax=+1;
    elseif pairIDX==2
    axmin=-.25;
    axmax=+.25;
    elseif pairIDX==3
    axmin=-.25;
    axmax=+.25;
    end
    
    colorscheme=parula(12);
    contourlevels=280:20:480;
    smoothnumbsurface=3;
    smoothnumbsect=3;
    smoothnumbcontour=7;
end

%=======DELTA other variables  ============

if ~strcmp(varlist{varIDX},'ph') && strcmp(deltamode,'on')
 field=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/StatsResults/%s_stats_updated_v1D_vertinterp.mat',varlist{varIDX}));
    field=eval(sprintf('field.delta_%s',pairs{pairIDX}));
end
 
fieldZERO=field; %overwrtie field which had previously values of percentage change with delta values. 

if strcmp(varlist{varIDX},'o2') && strcmp(deltamode,'on') 
    fieldZERO=(fieldZERO/1024)*1000000; %if it is oxygen, convert to umol O2 kg using mean sw rho of 1024kg m3
    varunit='\mumol O_2 kg^-^1';
    varname='\DeltaDO';
    axmin=-50;
    axmax=+50;
end

if strcmp(varlist{varIDX},'AOU') && strcmp(deltamode,'on') 
    varunit='\mumol O_2 kg^-^1';
    varname='\DeltaAOU';
    axmin=-50;
    axmax=+50;
end


if strcmp(varlist{varIDX},'templvl') && strcmp(deltamode,'on') 
    varunit='\circC';
    varname='\DeltaT';
    axmin=-2.5;
    axmax=+2.5;
end


field1=permute(fieldZERO,[3,1,2]);
% figure
% pcolor(squeeze(field1(1,:,:)))
% shading flat
% hold on
% scatter(110,246,'filled')

%unflipping
bit1=field1(:,:,122:end); %bit1 (200:end,:,:); %reshaping -- position 200 is where Pacific is
bit2=field1(:,:,1:121);   %bit2 (1:199,:,:);
field2=cat(3,bit1,bit2);
% figure
% pcolor(squeeze(field2(1,:,:))); shading flat

%% Loading stippling mask for surface plots

Pflags=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/StatsResults/%s_stats_updated_v1D_vertinterp.mat',varlist{varIDX})...
    ,sprintf('Pflags_%s',pairs{pairIDX}));
Pflags1=eval(sprintf('Pflags.Pflags_%s',pairs{pairIDX}));

% figure
% pcolor(Pflags1(:,:,1));shading flat
%
Ptest_ind=permute(Pflags1,[3,1,2]);
bit1=Ptest_ind(:,:,122:end); %bit1 (200:end,:,:); %reshaping -- position 200 is where Pacific is
bit2=Ptest_ind(:,:,1:121);   %bit2 (1:199,:,:);
Ptest_ind=cat(3,bit1,bit2);

%% loading Trecovery data
timevarname='Trecovery';
fieldTrec=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/all_ramps/filtered/Results/New_time_var_results/%s_results_updated_v1D_poly1_NEW_NEW.mat',varlist{varIDX}),timevarname);
fieldTrec=fieldTrec.Trecovery{2};

fieldTrec1=permute(fieldTrec,[3,1,2]); 
%unflipping
bit1=fieldTrec1(:,:,122:end); %bit1 (200:end,:,:); %reshaping -- position 200 is where Pacific is
bit2=fieldTrec1(:,:,1:121);   %bit2 (1:199,:,:);
fieldTrec2=cat(3,bit1,bit2);


for k=1:70
    for l=1:384
        for c=1:320
            if fieldTrec2(k,l,c)>480 || fieldTrec2(k,l,c)<0
                fieldTrec2(k,l,c)=NaN; %flag of no departure
            end
        end
    end
end


%% adjusting colormap
redblue(20)

redblueleo=zeros(10,3);
%shades of red
redblueleo(end,:)=[1,0.0514285714285714,0.0514285714285714];
redblueleo(end-1,:)=[1,0.214285714285714,0.214285714285714];
redblueleo(end-2,:)=[1,0.357142857142857,0.357142857142857];
redblueleo(end-3,:)=[1,0.500000000000000,0.500000000000000];
redblueleo(end-4,:)=[1,0.642857142857143,0.642857142857143];

%shades of blue
redblueleo(1,:)=[0.0714285714285714,0.0714285714285714,1];
redblueleo(2,:)=[0.214285714285714,0.214285714285714,1];
redblueleo(3,:)=[0.357142857142857,0.357142857142857,1];
redblueleo(4,:)=[0.500000000000000,0.500000000000000,1];
redblueleo(5,:)=[0.642857142857143,0.642857142857143,1];


%%%IF SURFACE MODE
if strcmp(surfacemode,'on') 
k=[find(depth_interp==100) find(depth_interp==300) find(depth_interp==500)];
elseif ~strcmp(surfacemode,'on')  
%%%IF MESOPELAGIC MODE
k=[find(depth_interp==500) find(depth_interp==1000) find(depth_interp==2000)];
%pcolor(squeeze(Ptest_ind(k,:,:)));shading flat
end



%% figure

hFig = figure(1); clf;
set(hFig, 'Units','centimeters','Position',[0,0,31,23]); %
set(hFig,'visible','on')
m1=0.01; p1=.002; s1=0.005; pb=.02; 



% surface1
subaxis(3,2,1,'Spacing',s1,'Padding',p1,'Margin',m1,'PaddingBottom',2*pb,'PaddingRight',0,'PaddingTop',15*p1,'PaddingLeft',0);
m_proj('Lambert','lon',[-100 15], 'lat', [-10 65]) ;
get(gca,'Position');
set(gca,'Position',[-0.08    0.7067-0.01    0.4875    0.2533]);

field3=squeeze(field2(k(1),:,:));

for smoothnumber=1:smoothnumbsurface %smoothing 7 times at the specific depth level
    field3=smooth2d(field3);
end

m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(field3)); shading interp

hold on
%drawing stippling
mask1=Ptest_ind(k(1),:,:)==100;%% these are the stippling masks for each depth level k
mask1=squeeze(mask1);
m_scatter(plon(mask1),plat(mask1),2,[0 0 0],'.')

colormap(redblue(20))
%colormap(redblueleo)
caxis([axmin axmax])
m_coast('patch',[.7 .7 .7],'edgecolor','k');

% %Tcontour
% fieldTrec3=squeeze(fieldTrec2(k(1),:,:));
% for smoothnumber=1:smoothnumbcontour %smoothing 7 times at the specific depth level
%     fieldTrec3=smooth2d(round(fieldTrec3));
% end
% [M,c]=m_contour(shift_atl(plon),shift_atl(plat),shift_atl(fieldTrec3));
% 
% c.LevelList=contourlevels; %levels
% c.LineColor=[.7 .7 .7];%'grey';
% c.LineStyle='-';
% clabel(M,c,'FontSize',10) 


%Grid box
m_grid('xtick',-90:30:0,'ytick',0:30:90,'backgroundcolor','k');
set(gca,'layer','top','FontName','Helvetica','fontsize',8,'TickDir','in','TickLength',[.015 .02],'XMinorTick','on','YMinorTick','on');
hold on
m_line(plon(atl_ind),plat(atl_ind),'Color',[102/255 0 0 ],'LineWidth',2);
ttl=title(['a) ' sprintf('%s at %i m',varname,round(depth_interp(k(1))))],'fontsize',14);
set(ttl,'Position',[1.42467418893584e-06,0.83,-4.50359962737050e+15])


%%
subaxis(3,2,3,'Spacing',s1,'Padding',p1,'Margin',m1,'PaddingBottom',2*pb,'PaddingRight',0,'PaddingTop',15*p1,'PaddingLeft',0);
m_proj('Lambert','lon',[-100 15], 'lat', [-10 65]);
set(gca,'Position',[-0.08     0.3783-0.01    0.4875    0.2533]);

field3=squeeze(field2(k(2),:,:));

for smoothnumber=1:7 %smoothing 7 times at the specific depth level
    field3=smooth2d(field3);
end


m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(field3)); shading interp

hold on
%drawing stippling
mask1=Ptest_ind(k(2),:,:)==100;%% these are the stippling masks for each depth level k
mask1=squeeze(mask1);
m_scatter(plon(mask1),plat(mask1),2,[0 0 0], '.')

colormap(redblueleo)
caxis([axmin axmax])
m_coast('patch',[.7 .7 .7],'edgecolor','k');

% %Tcontour
% fieldTrec3=squeeze(fieldTrec2(k(2),:,:));
% for smoothnumber=1:smoothnumbcontour %smoothing 7 times at the specific depth level
%     fieldTrec3=smooth2d(round(fieldTrec3));
% end
% [M,c]=m_contour(shift_atl(plon),shift_atl(plat),shift_atl(fieldTrec3));
% 
% c.LevelList=contourlevels; %levels
% c.LineColor=[.7 .7 .7];%'grey';
% c.LineStyle='-';
% clabel(M,c,'FontSize',10) 


%Grid box
%m_grid('backgroundcolor',[.2 .2 .2],'box','fancy','tickdir','out'); %if brown [160/255 82/255 45/255]
m_grid('xtick',-90:30:0,'ytick',0:30:90,'backgroundcolor','k');
set(gca,'layer','top','FontName','Helvetica','fontsize',8,'TickDir','in','TickLength',[.015 .02],'XMinorTick','on','YMinorTick','on');
hold on
m_line(plon(atl_ind),plat(atl_ind),'Color',[102/255 0 0 ],'LineWidth',2);
ttl=title(['b) ' sprintf('%s at %i m',varname,round(depth_interp(k(2))))],'fontsize',14);
set(ttl,'Position',[1.42467418893584e-06,0.83,-4.50359962737050e+15])


%%
%% surface 3
subaxis(3,2,5,'Spacing',s1,'Padding',p1,'Margin',m1,'PaddingBottom',2*pb,'PaddingRight',0,'PaddingTop',15*p1,'PaddingLeft',0);
m_proj('Lambert','lon',[-100 15], 'lat', [-10 65]) ;
set(gca,'Position',[-0.08    0.0500-0.01    0.4875    0.2533]);


field3=squeeze(field2(k(3),:,:));

for smoothnumber=1:7 %smoothing 7 times at the specific depth level
    field3=smooth2d(field3);
end


m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(field3)); shading interp

hold on
%drawing stippling
mask1=Ptest_ind(k(3),:,:)==100;%% these are the stippling masks for each depth level k
mask1=squeeze(mask1);
m_scatter(plon(mask1),plat(mask1),2,[0 0 0], '.')


colormap(redblueleo)
caxis([axmin axmax])
m_coast('patch',[.7 .7 .7],'edgecolor','k');

% %Tcontour
% fieldTrec3=squeeze(fieldTrec2(k(3),:,:));
% for smoothnumber=1:smoothnumbcontour %smoothing 7 times at the specific depth level
%     fieldTrec3=smooth2d(round(fieldTrec3));
% end
% [M,c]=m_contour(shift_atl(plon),shift_atl(plat),shift_atl(fieldTrec3));
% 
% c.LevelList=contourlevels; %levels
% c.LineColor=[.7 .7 .7];%'grey';
% c.LineStyle='-';
% clabel(M,c,'FontSize',10) 


%Grid box
m_grid('xtick',-90:30:0,'ytick',0:30:90,'backgroundcolor','k');
set(gca,'layer','top','FontName','Helvetica','fontsize',8,'TickDir','in','TickLength',[.015 .02],'XMinorTick','on','YMinorTick','on');
hold on
m_line(plon(atl_ind),plat(atl_ind),'Color',[102/255 0 0 ],'LineWidth',2)
ttl=title(['c) ' sprintf('%s at %i m',varname,round(depth_interp(k(3))))],'fontsize',14);
set(ttl,'Position',[1.42467418893584e-06,0.83,-4.50359962737050e+15])



%% section
subplot(3,2,[2,4,6])
get(gca,'Position');


%Obtaining section fields
[sectfield,sect_lat,sect_dep]=makesect_atl_noresm_fixed_z(field2,depth_interp,plat,atl_ind); %this gets the orginial section field
dummyfield=sectfield;


% [depthQ, platQ]=meshgrid(depth_interp, (min(min(sect_lat2)):1:max(max(sect_lat2)))); %this makes the new grid
%
% dummyfield=griddata(sect_depth,sect_lat2,sectfield,depthQ,platQ) %%this calculates interpolated values


for smoothnumber=1:smoothnumbsect %smoothing 7 times at the specific depth level
    dummyfield=smooth2d(dummyfield);
end
pcolor(sect_lat,-sect_dep,sectfield)
%makesect_atl_noresm_fixed_z(field2,depth_interp,plat,atl_ind); %the function asks for (field, depth, lat or long, file with section indexes)
shading interp
hold on


%% stippling mask for section
%%1st make section field
[sect_lat2,sect_depth,section_stpl_field]=makesect_atl_noresm_fixed_z_stippl(Ptest_ind,depth_interp,plat,atl_ind);
mask2=section_stpl_field==100; %this is the stippling mask where significant differences were labelled with flag=100

stipple(sect_lat,-sect_dep,mask2,'color','k','density',1000,'markersize',3);
hold on


% %Tcontour
% [sect_ctr_original,sect_lat,sect_dep]=makesect_atl_noresm_fixed_z(fieldTrec2,depth,plat,atl_ind); %the function asks for (field, depth, lat or long, file with section indexes)
% 
% for smoothnumber=1:smoothnumbcontour %smoothing 7 times
%     sect_ctr_original=smooth2d(sect_ctr_original);
% end
% 
% [M,c]=contour(sect_lat,-sect_dep,sect_ctr_original)
% %c.LevelList=contourlevels; %levels
% c.LineColor=[.7 .7 .7];%'grey';
% clabel(M,c,'Color','w','FontSize',10);
% clabel(M,c,'LabelSpacing',300)
% c.LineStyle='-';
hold on

%%colormapping
colormap(redblueleo)
caxis([axmin axmax])
h=colorbar;
h.Location='southoutside';

unit='%';

if varIDX==1 && strcmp(deltaphmode,'off') && (strcmp(surfacemode,'on') || ~strcmp(surfacemode,'on'))
    h.Label.String=sprintf('Mean percentage change of [ H^+] (%s)',unit);
elseif varIDX==1 && strcmp(deltaphmode,'on') && (strcmp(surfacemode,'on') || ~strcmp(surfacemode,'on'))
        h.Label.String=sprintf('Mean %s ',varname);
elseif strcmp(deltamode,'on') && (strcmp(surfacemode,'on') || ~strcmp(surfacemode,'on'))
        h.Label.String=sprintf('Mean %s (%s)',varname,varunit);
elseif strcmp(deltaphmode,'off') && (strcmp(surfacemode,'on') || ~strcmp(surfacemode,'on'))
    h.Label.String=sprintf('Mean percentage change of %s (%s)',varname,unit);
end

xlabel('Latitude ( \circ )')
ylabel('Depth ( m )')
title(['d) ' sprintf('%s', period)],'fontsize',14)

if varIDX==4 %if detoc then change axis to end at 2000m deep
    axis([-10 65 -2000 0])
else
    axis([-10 65 -6000 0])
end

set(gca, 'color','k','layer','top','FontName','Helvetica','fontsize',12,'TickDir','in','TickLength',[.005 .01],'XMinorTick','on','YMinorTick','on','xtick',-10:10:60,'xticklabel',{'10\circS' '0\circ' '10\circN' '20\circN' '30\circN' '40\circN' '50\circN' '60\circN'  });
set(gca,'Position',[0.40    0.18    0.55    0.77]);

h.TickLabels(1)={['\leq ' num2str(axmin)]};
h.TickLabels(end)={['\geq ' num2str(axmax)]};

colormap(redblue(50))

%saving image

set(gcf,'color','w')
hFig.InvertHardcopy='off';
supersizeme(1.2)
if strcmp(deltaphmode,'on') || strcmp(deltamode,'on')
print(sprintf('/Users/leonardobertini/OneDrive/IMBRSea/Modules/Master Thesis/Figures_Paper_LEO_NORESM/new_figs/Change_%s_period_%d.png',varname(2:end),pairIDX),'-dpng','-r300')
elseif ~strcmp(surfacemode,'on') 
print(sprintf('/Users/leonardobertini/OneDrive/IMBRSea/Modules/Master Thesis/Figures_Paper_LEO_NORESM/new_figs/Percentage_Change_%s_period_%d.png',varlist{varIDX},pairIDX),'-dpng','-r300')
elseif strcmp(surfacemode,'on') 
print(sprintf('/Users/leonardobertini/OneDrive/IMBRSea/Modules/Master Thesis/Figures_Paper_LEO_NORESM/new_figs/Percentage_Change_%s_period_%d_surface_layers.png',varlist{varIDX},pairIDX),'-dpng','-r300')
end



clear redblueleo h
close all
end



