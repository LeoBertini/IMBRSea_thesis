%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SECTIONS T-recovery %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Produces stressor fields plots (SECTIONS AND SURFACES) for the variables outisde envelope at end of
%experiment) interpolatating in latsetp deg of latitude, lonstep deg of longitude
%and at every depthstep m of depth

% Author: Leonardo Bertini
% Modified: 29th April 2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CODE
%% latitudinal section

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

%getting zonal section
bit1=domain_mask(:,122:end); %bit1 (200:end,:,:); %reshaping -- position 200 is where Pacific is
bit2=domain_mask(:,1:121);   %bit2 (1:199,:,:);
dm_mask1=cat(2,bit1,bit2);

atl_ind=load('/Volumes/LaCie_Leonardo/NorESM/scripts_jerry/noresm_atl_ind_fixed.asc');

varlist={'ph';'o2';'omegac';'detoc';'templvl';'AOU'}  ;

for varIDX=1:6

 if strcmp(varlist{varIDX},'ph')==1
    varname='pH';
    smoothnumbsurface=3;
    smoothnumbsect=3;
    smoothnumbcontour=7;
    isolinessct=[200 400 480 650 750 1000];
    isolinessfc=[200 300 500 700];

elseif strcmp(varlist{varIDX},'o2')==1
    varname='DO';
    smoothnumbsurface=3;
    smoothnumbsect=3;
    smoothnumbcontour=7;
    
elseif strcmp(varlist{varIDX},'omegac')==1
    varname='\Omega_C';
    smoothnumbsurface=3;
    smoothnumbsect=3;
    smoothnumbcontour=7;
    isolinessct=[160 200 320 380 400 440 480 500 750 1000];
    isolinessfc=[200 300 500 700];

elseif strcmp(varlist{varIDX},'detoc')==1
    varname='POC flux';
    smoothnumbsurface=3;
    smoothnumbsect=3;
    smoothnumbcontour=7;
  isolinessct=[160 200 320 380 400 440 480 500 750 1000];
    isolinessfc=[200 300 500 700];
    
elseif strcmp(varlist{varIDX},'templvl')==1
    varname='T';
    smoothnumbsurface=3;
    smoothnumbsect=3;
    smoothnumbcontour=10;
     isolinessct=[160 200 320 380 400 440 480 500 650 750 1000];
    isolinessfc=[200 300 500 700];

elseif strcmp(varlist{varIDX},'AOU')==1
    varname='AOU';
    smoothnumbsurface=3;
    smoothnumbsect=5;
    smoothnumbcontour=7;
    isolinessct=[160 200 380 480 750 ];
    isolinessfc=[140 480 700];

 
 end

 %COMMON PARAMETERS FOR ALL FIGURES
    axmin=120;
    axmax=720;
    varunit='Time (yr)';
    
    %colorscheme=parula(10);
    %%%option 2
    pinky=pink(20);
    use=pinky(5:10,:);use=flip(use);

    % const=255;
    % pinky=[251/const 114/const 101/const ;...
    %         235/const 87/const 77/const;...
    %         220/const 60/const 54/const;...
    %         204/const 33/const 30/const;...
    %         201/const 28/const 25/const;...
    %         189/const 7/const 7/const];
    dummycolor=parula(13);
    colorscheme=[dummycolor(2:10,:);use];
    
 
 
fields=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/all_ramps/filtered/Results/New_time_var_results/%s_results_updated_v1D_poly1_new.mat',varlist{varIDX}));

%Possible flags included in the Trecovery dataset
badflagsnegslope=-7777;
badflagtooearly=-4444;
bottomflag=-1111;


%% Importing Trevocery dataset
fieldZERO=fields.Trecovery;

fieldZERO=real(fieldZERO);
fieldZERO(fieldZERO==badflagsnegslope)=120; %this makes all different flags the same value %Definite no recovery
fieldZERO(fieldZERO==badflagtooearly)=120; %this makes all different flags the same value %Gets

field1=permute(fieldZERO,[3,1,2]);

% figure
% pcolor(squeeze(field1(55,:,:)))
% shading flat
% colorbar
% caxis([160 560])


%unflipping
bit1=field1(:,:,122:end); %bit1 (200:end,:,:); %reshaping -- position 200 is where Pacific is
bit2=field1(:,:,1:121);   %bit2 (1:199,:,:);
field2=cat(3,bit1,bit2);
field2(field2<0 & field2~=bottomflag)=600; %this makes all negative flags which are not bottom positive to indicate very late or no recovery

% figure
% pcolor(squeeze(field2(55,:,:)))
% shading flat
% colorbar
% caxis([160 560])


k = [find(depth==500) find(depth==1000) find(depth==2000)]; %depth levels for surface plots

hFig = figure(1); clf;
set(hFig, 'Units','centimeters','Position',[0,0,35,23]); %
set(hFig,'visible','off')
m1=0.01; p1=.002; s1=0.005; pb=.02; 

%% surface 1
subaxis(3,2,1,'Spacing',s1,'Padding',p1,'Margin',m1,'PaddingBottom',2*pb,'PaddingRight',0,'PaddingTop',15*p1,'PaddingLeft',0);
m_proj('Lambert','lon',[-100 15], 'lat', [-10 65]) ;
get(gca,'Position');
set(gca,'Position',[-0.08    0.7067-0.01    0.4875    0.2533]);
field3=squeeze(field2(k(1),:,:));
field3(field3==-1111)=NaN; %this makes all bottom NaN and removes it before the smoothing

for smoothnumber=1:smoothnumbsurface %smoothing 3 times at the specific depth level
    field3=smooth2d(field3);
end

m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(field3)); shading interp
colormap(colorscheme)
caxis([axmin axmax])
m_coast('patch',[.7 .7 .7],'edgecolor','k');
hold on

%Tcontour
[M,c]=m_contour(shift_atl(plon),shift_atl(plat),shift_atl(field3));
c.LevelList=isolinessfc; %levels
c.LineColor=[.9 .9 .9];%'grey';
c.LineStyle='-';
clabel(M,c,'FontSize',8) 

%Grid box
m_grid('xtick',-90:30:0,'ytick',0:30:90,'backgroundcolor','k');
set(gca,'layer','top','FontName','Helvetica','fontsize',8,'TickDir','in','TickLength',[.015 .02],'XMinorTick','on','YMinorTick','on');
hold on
m_line(plon(atl_ind),plat(atl_ind),'Color',[102/255 0 0 ],'LineWidth',2)
ttl=title(['a)' sprintf(' %s Trec at %i m',varname,round(depth(k(1))))],'fontsize',14)
set(ttl,'Position',[1.42467418893584e-06,0.90,-4.50359962737050e+15])
if varIDX==3
set(ttl,'Position',[1.42467418893584e-06,0.83,-4.50359962737050e+15])
end

%% surface2
subaxis(3,2,3,'Spacing',s1,'Padding',p1,'Margin',m1,'PaddingBottom',2*pb,'PaddingRight',0,'PaddingTop',15*p1,'PaddingLeft',0);
m_proj('Lambert','lon',[-100 15], 'lat', [-10 65]);
set(gca,'Position',[-0.08     0.3783-0.01    0.4875    0.2533]);

field3=squeeze(field2(k(2),:,:));
field3(field3==-1111)=NaN; %this makes all bottom NaN and removes it before the smoothing

for smoothnumber=1:smoothnumbsurface %smoothing 3 times at the specific depth level
    field3=smooth2d(field3);
end

m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(field3)); shading interp
colormap(colorscheme)
caxis([axmin axmax])
m_coast('patch',[.7 .7 .7],'edgecolor','k');

hold on

%Tcontour
[M,c]=m_contour(shift_atl(plon),shift_atl(plat),shift_atl(field3));
c.LevelList=isolinessfc; %levels
c.LineColor=[.9 .9 .9];%'grey';
c.LineStyle='-';
clabel(M,c,'FontSize',8) 

%Grid box
m_grid('xtick',-90:30:0,'ytick',0:30:90,'backgroundcolor','k');
set(gca,'layer','top','FontName','Helvetica','fontsize',8,'TickDir','in','TickLength',[.015 .02],'XMinorTick','on','YMinorTick','on');
hold on
m_line(plon(atl_ind),plat(atl_ind),'Color',[102/255 0 0 ],'LineWidth',2)
ttl=title(['b)' sprintf(' %s Trec at %i m',varname,round(depth(k(2))))],'fontsize',14)
set(ttl,'Position',[1.42467418893584e-06,0.90,-4.50359962737050e+15])
if varIDX==3
set(ttl,'Position',[1.42467418893584e-06,0.83,-4.50359962737050e+15])
end

%% surface3
subaxis(3,2,5,'Spacing',s1,'Padding',p1,'Margin',m1,'PaddingBottom',2*pb,'PaddingRight',0,'PaddingTop',15*p1,'PaddingLeft',0);
m_proj('Lambert','lon',[-100 15], 'lat', [-10 65]) ;
set(gca,'Position',[-0.08    0.0500-0.01    0.4875    0.2533]);

field3=squeeze(field2(k(3),:,:));
field3(field3==-1111)=NaN; %this makes all bottom NaN and removes it before the smoothing

for smoothnumber=1:smoothnumbsurface %smoothing 3 times at the specific depth level
    field3=smooth2d(field3);
end

m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(field3)); shading interp

colormap(colorscheme)
caxis([axmin axmax])
m_coast('patch',[.7 .7 .7],'edgecolor','k');

hold on
%Tcontour
[M,c]=m_contour(shift_atl(plon),shift_atl(plat),shift_atl(field3));

c.LevelList=isolinessfc; %levels
c.LineColor=[.9 .9 .9];%'grey';
c.LineStyle='-';
clabel(M,c,'FontSize',8) 

%Grid box
m_grid('xtick',-90:30:0,'ytick',0:30:90,'backgroundcolor','k');
set(gca,'layer','top','FontName','Helvetica','fontsize',8,'TickDir','in','TickLength',[.015 .02],'XMinorTick','on','YMinorTick','on');
hold on
m_line(plon(atl_ind),plat(atl_ind),'Color',[102/255 0 0 ],'LineWidth',2)
ttl=title(['c)' sprintf(' %s Trec at %i m',varname,round(depth(k(3))))],'fontsize',14)
set(ttl,'Position',[1.42467418893584e-06,0.90,-4.50359962737050e+15])
if varIDX==3
set(ttl,'Position',[1.42467418893584e-06,0.83,-4.50359962737050e+15])
end

%% section
subplot(3,2,[2,4,6])
get(gca,'Position');

field2(field2==bottomflag)=NaN; %this makes the points representing the bottom 'transparent'

%changing all the negative flags to avoid the concentric rings when interpolating
field2(field2<0)=600;

%getting section data
[sect_fld_original,sect_lat,sect_dep]=makesect_atl_noresm_fixed_z(field2,depth,plat,atl_ind); %the function asks for (field, depth, lat or long, file with section indexes)

for smoothnumber=1:smoothnumbsect %smoothing 7 times
    sect_fld_original=smooth2d(sect_fld_original);
end

%section plot
pcolor(sect_lat,-sect_dep,sect_fld_original);
shading interp
hold on

%contourlines
[sect_ctr_original,sect_lat,sect_dep]=makesect_atl_noresm_fixed_z(field2,depth,plat,atl_ind); %the function asks for (field, depth, lat or long, file with section indexes)

for smoothnumber=1:smoothnumbcontour %smoothing 7 times
    sect_ctr_original=smooth2d(sect_ctr_original);
end

[M,c]=contour(sect_lat,-sect_dep,sect_ctr_original);

c.LevelList=isolinessct; %levels
c.LineColor=[.9 .9 .9];%'grey';
c.LineStyle='-';
c.LineWidth=1.5;
clabel(M,c,'FontSize',10) 

colormap(colorscheme)
caxis([axmin axmax])

h=colorbar;
h.Location='southoutside';
h.Label.String=sprintf('Time of Recovery (yr)');
h.Ticks=[120 160:40:axmax];

%h.TickLabels(1)={'Early Rcvr'};
h.TickLabels(end)={['\geq' num2str(axmax)]};
h.FontSize
%h.TickLabels(end)={'\geq Simulation window'};

xlabel('Latitude (\circ)')
ylabel('Depth (m)')
title(['d) ' sprintf('%s Trec', varname) ': North Atlantic section' ],'fontsize',14)

if varIDX==4 %if detoc then change axis to end at 2000m deep
    axis([-10 65 -2000 0]) 
else
    axis([-10 65 -6000 0])
end

set(gca, 'color','k','layer','top','FontName','Helvetica','fontsize',12,'TickDir','in','TickLength',[.005 .01],'XMinorTick','on','YMinorTick','on','xtick',-10:10:60,'xticklabel',{'10\circS' '0\circ' '10\circN' '20\circN' '30\circN' '40\circN' '50\circN' '60\circN'  });
set(gca,'Position',[0.40    0.18    0.55    0.77]);

%saving image
set(gcf,'color','w')
hFig.InvertHardcopy='off';
supersizeme(1.2)

print(sprintf('/Users/leonardobertini/OneDrive/IMBRSea/Modules/Master Thesis/Figures_Paper_LEO_NORESM/new_figs/Trec_surfc_%s.png',varlist{varIDX}),'-dpng','-r300')

close all
end
%hold on
%yline(-depth(k),'--','LineWidth',1) %this line is not drawing
%print(sprintf('/Volumes/LaCie_Leonardo/NorESM/Figures/Section_North_Atlantic_%s_%s.png',timevarname,varlist{varIDX}),'-dpng','-r300')


%####### CODE FOR GIF GENERATION
%     filenameind =sprintf('/Volumes/LaCie_Leonardo/NorESM/Figures/Section_North_Atlantic_%s_%s_depth_%.2f.png',timevarname,varlist{varIDX},depth(k));
%     print(filenameind,'-dpng','-r300')

%     drawnow
%         frame = getframe(fig);
%         im{k} = frame2im(frame);
%         hold all
%
%     pause( 0.5 )


% filename =sprintf('Section_North_Atlantic_Trecovery_%s.gif',varlist{varIDX}); % Specify the output file name
%
%  for k = 1:3%: (length(depth)-1)
%     [A,map] = rgb2ind(im{k},256);
%     if k == 1
%         imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
%     else
%         imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
%     end
%  end
%  close all
