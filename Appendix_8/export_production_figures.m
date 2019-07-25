%% CODE

parea=ncread('/Volumes/LaCie_Leonardo/NorESM/all_ramps/grid_gx1v6.nc','parea');
plon=ncread('/Volumes/LaCie_Leonardo/NorESM/all_ramps/grid_gx1v6.nc','plon') ;
plat=ncread('/Volumes/LaCie_Leonardo/NorESM/all_ramps/grid_gx1v6.nc','plat') ;
parea=parea';
plat=plat';
plon=plon';
plon=plon;

addpath(genpath('/Volumes/LaCie_Leonardo/NorESM/'))



domain_mask=load ('/Volumes/LaCie_Leonardo/NorESM/Mask_regions/Mask_corrected.mat');
domain_mask=domain_mask.domain_mask;

atl_ind=load('/Volumes/LaCie_Leonardo/NorESM/scripts_jerry/noresm_atl_ind_fixed.asc');

varlist={'epc100'}  ;
varIDX=1
    
    varname='POC export';
    varunit='mgC m^-^2 day^-^1';
    axmin=0;
    axmax=250;
    colorscheme=parula(10);
    contourlevels=[50 100 150 200 250 300 400];
    smoothnumbsurface=3;
    smoothnumbsect=1;
    smoothnumbcontour=7;
    
    
 %% PI_mean

    PI_mean=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/Pre_industrial/30yrMeans/PI_%s_first_5yr_mean.mat',varlist{varIDX}));
    PI_mean=PI_mean.PI_mean;
    PI_mean=permute(PI_mean,[2,1]);
    
    %unflipping domain mask
    bit1=domain_mask(:,122:end); %bit1 (200:end,:,:); %reshaping -- position 200 is where Pacific is
    bit2=domain_mask(:,1:121);   %bit2 (1:199,:,:);
    domain_mask_1=cat(2,bit1,bit2);
    
    for l=1:384
        for c=1:320
            if domain_mask_1(l,c)~=2
                PI_mean(l,c)=NaN;
            end
        end
    end


    PI_sd=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/Pre_industrial/30yrMeans/PI_%s_30yr_sd.mat',varlist{varIDX}));
    PI_sd=PI_sd.PI_sd;   
    PI_sd=permute(PI_sd,[2,1]);
    
    
    
%% initial figure - time variables spread depending on the chosen SD

 %latitudinal mean
fields=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/all_ramps/filtered/Results/New_time_var_results/%s_results_updated_v1D_poly1_NEW_NEW.mat',varlist{varIDX}));

TOD=fields.TOD;


Trec=fields.Trecovery;
Trec{1,1}(Trec{1,1}<0)=NaN;
Trec{1,2}(Trec{1,2}<0)=NaN;
Trec{1,3}(Trec{1,3}<0)=NaN;

Trec{1,1}(Trec{1,1}>1000)=600;
Trec{1,2}(Trec{1,2}>1000)=600;
Trec{1,3}(Trec{1,3}>1000)=600;

%calculating area weighted mean export production per latitude

% % numerator=nan(384/2);
% % denominator=nan(384/2);
% 
% for idx=1:3
%     
%     for l=1:384
%         for c=1:320
%             numerator1(l,c)=TOD{1,idx}(l,c)*parea(l,c);
%             numerator2(l,c)=Trec{1,idx}(l,c)*parea(l,c);
%         end
%     end   
%     
%   for l=384/2:384
%         dummy1(l-384/2+1)=nansum(numerator1(l-384/2+1,:))/nansum(parea(l-384/2+1,:));
%         dummy2(l-384/2+1)=nansum(numerator2(l-384/2+1,:))/nansum(parea(l-384/2+1,:));
%         dummylat(l-384/2+1)=nanmean(plat(l-384/2+1,:));
%   end
%     
%     
%     
%     TODlatitudinal{idx}=dummy1;
%     Treclatitudinal{idx}=dummy2;
%     meanlat=dummylat;
%     clear dummy*
% end




for idx=1:3
for l=384/2:384
    
  
    dumb(l-384/2+1)=nanmean(TOD{1,idx}(l,:));
    dumb2(l-384/2+1)=nanmean(Trec{1,idx}(l,:));
    meanlat(l-384/2+1)=nanmean(plat(l,:));

end
    TODlatitudinal{idx}=dumb(~isnan(dumb));
    Treclatitudinal{idx}=dumb2(~isnan(dumb2));

end


% 
% figure
% plot(meanlat(1:381),TODlatitudinal{1,1})
% hold on
% plot(meanlat,dummyresult{1,2})
% hold on
% plot(meanlat,dummyresult{1,3})


fig=figure('Position',[ 66   255   742   450])

axmaxred=75;
axmaxblue=75;

%filled area
    y1=smooth(TODlatitudinal{1,1}(1:190));
    y3=smooth(TODlatitudinal{1,3}(1:190));
    Y = [y1(:).', fliplr(y3(:).')];
    X = [meanlat(1:190), fliplr(meanlat(1:190))];
    shaded1=fill(X, Y,[255/255 204/255 204/255])
    shaded1.FaceAlpha=0.6;
    shaded1.LineStyle='none';      

y2=smooth(TODlatitudinal{1,2});
    p2=line(meanlat(1:190),y2(1:190),'linewidth',2,'color','r','LineStyle','-'); %TOD mean       
      
 grid on
    ax1 = gca; % current axes
    ax1.Position=[    0.1307    0.1100    0.7743    0.6583]
    ax1.XColor = 'r';
    ax1.YColor = 'r';
    ax1.LineWidth=2;
    ax1.XLim=[0 axmaxred];
    ax1.YLim=[0 200];
    ax1.GridColor=[.5 .5 .5];
    ax1.YMinorTick='on';
    ax1.MinorGridColor=[.5 .5 .5];
    ax1.GridAlpha=0.4;
    ax1.FontName='Helvetica Neue';
    ax1.FontSize=16;
    ax1.YLabel.String='ToD (yr)';
    ax1.XLabel.String='Latitude';
    ax1.XTick=linspace(0,axmaxred,6);
    ax1.XTickLabel=strcat(ax1.XTickLabel,'\circN');
    ax1.YTick=linspace(0,200,6);

    grid minor   
    
    
     ax1_pos = ax1.Position; % position of first axes
    ax2 = axes('Position',ax1_pos,...
        'XAxisLocation','top',...
        'YAxisLocation','right',...
        'Color','none');

    hold on
    
    y10=smooth(Treclatitudinal{1,1}(1:190));
    y30=smooth(Treclatitudinal{1,3}(1:190));
    Y = [y10(:).', fliplr(y30(:).')];
    X = [meanlat(1:190), fliplr(meanlat(1:190))];
    shaded1=fill(X, Y,[102/255 178/255 225/255])
    shaded1.FaceAlpha=0.4;
    shaded1.LineStyle='none'; 
   
    hold on
    y20=smooth(Treclatitudinal{1,2}(1:190));
    p20=line(meanlat(1:190),y20,'linewidth',2,'color','b','LineStyle','-'); %TREC
   
    ax2.Position=[    0.1307    0.1100    0.7743    0.6583]
    ax2.FontName='Helvetica Neue';
    ax2.FontSize=16;
    ax2.LineWidth=2;
    
    ax2.YLabel.String='Trec (yr)';
    ax2.XLabel.String='Latitude';
    ax2.XColor = 'b';
    ax2.YColor = 'b';
    ax2.XLim=[0 axmaxblue];
    ax2.YLim=[0 1000];
    ax2.GridColor=[.5 .5 .5];
    ax2.YMinorTick='on';
    ax2.MinorGridColor=[.5 .5 .5];
    ax2.GridAlpha=0.4;
    ax2.XTick=linspace(0,axmaxblue,6);
    ax2.XTickLabel=strcat(ax2.XTickLabel,'\circN');
    ax2.YTick=linspace(0,1000,6);
    grid minor
    
    lgd=legend([p2 p20],'ToD','Trec')
    lgd.Location='southeast';
    lgd.Color=[.8 .8 .8];
    
    varname='POC_e_x_p_o_r_t'; 
    title(sprintf('a) NAtl Departure and Recovery time scales of %s : meridional distribution',varname),'FontSize',14)
    uistack(lgd,'top')
    
    set(gcf,'color','w')
    fig.InvertHardcopy = 'off';
    print('/Volumes/LaCie_Leonardo/NorESM/export_production/meridional_profile_ep100_smooth.png','-dpng','-r300')

%%%%%%======================%%%%% %%%%%%======================%%%%%        
       

%%
hFig = figure(1); clf;
set(hFig, 'Units','centimeters','Position',[0 0 37.0769   20.7433]); %
set(hFig,'color','w','visible','on')
m1=0.01; p1=.002; s1=0.005; pb=.02; 

%% surface 1 --> PImean
subaxis(2,2,1,'Spacing',s1,'Padding',p1,'Margin',m1,'PaddingBottom',2*pb,'PaddingRight',0,'PaddingTop',15*p1,'PaddingLeft',0);
m_proj('Lambert','lon',[-100 15], 'lat', [-10 65]) ;
get(gca,'Position');
set(gca,'Position',[0.01    0.6067-0.01    0.4875    0.35]);

field3=PI_mean(:,:);

for smoothnumber=1:smoothnumbsurface %smoothing 7 times at the specific depth level
    field3=smooth2d(field3);
end

%from molC m-2 s-1 to mgC m-2 day-1
m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(field3*12.0107*1000*24*3600)); shading interp

colormap(gca,colorscheme) % [177/255 47/255 245/255]
caxis([axmin axmax])

m_coast('patch',[.7 .7 .7],'edgecolor','k');
hold on

c=colorbar
c.Location='westoutside';
c.Label.String=[{sprintf('Pre-industrial mean %s',varname); sprintf('(%s)',varunit)}];
c.TickLabels(end)={['\geq ' num2str(axmax)]};


%Grid box
m_grid('backgroundcolor','k','tickdir','out'); %if brown [160/255 82/255 45/255]
hold on
m_line(plon(atl_ind),plat(atl_ind),'Color',[102/255 0 0 ],'LineWidth',2)
title({sprintf('Depth = %.2f m',100)})
hold on

%% surface 2 --> TOD

subaxis(2,2,2,'Spacing',s1,'Padding',p1,'Margin',m1,'PaddingBottom',2*pb,'PaddingRight',0,'PaddingTop',15*p1,'PaddingLeft',0);
m_proj('Lambert','lon',[-100 15], 'lat', [-10 65]) ;
get(gca,'Position');
set(gca,'Position',[0.48   0.6067-0.01    0.4875    0.35]);
m_proj('Lambert','lon',[-100 20], 'lat', [-10 65]) 

timevarname='TOD';
varunit2='yr';
fieldTOD=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/all_ramps/filtered/Results/New_time_var_results/%s_results_updated_v1D_poly1_NEW_NEW.mat',varlist{varIDX}),timevarname);
fieldTOD=fieldTOD.TOD{1,2};

axmin2=-250;
axmax2=250;

%unflipping
bit1=fieldTOD(:,122:end); %bit1 (200:end,:,:); %reshaping -- position 200 is where Pacific is
bit2=fieldTOD(:,1:121);   %bit2 (1:199,:,:);
fieldTOD2=cat(2,bit1,bit2);

% figure
% pcolor(fieldTOD2);shading flat
 
    for l=1:384
        for c=1:320
            if domain_mask_1(l,c)==2 && ~isnan(PI_mean(l,c)) && isnan(fieldTOD2(l,c))
                fieldTOD2(l,c)=-200; %flag of no departure
            end
        end
    end
    
    


field3=fieldTOD2;

for smoothnumber=1:4 %smoothing 7 times at the specific depth level
    field3=smooth2d(field3);
end

%surface
m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(field3)); shading flat


colorscheme2=parula(12);
colornodept=[204/255 229/255 255/255];
colormap(gca,[colornodept;colornodept;colornodept;colornodept;colornodept;colorscheme2(8:12,:)])
caxis([axmin2 axmax2])

h=colorbar;
h.Location='westoutside';
h.Label.Position=[-1.137500047683716,0.000238418579102,0];
h.Label.String=[{sprintf('ToD %s (%s) \n',varname,varunit2)}];
h.Ticks=[-250:50:250];
h.TickLabels(end)={['\geq 250']};
h.TickLabels(1:5)={''};

m_coast('patch',[.7 .7 .7],'edgecolor','k');


hold on
%Grid box
m_grid('backgroundcolor','k','tickdir','out'); %if brown [160/255 82/255 45/255]
hold on
m_line(plon(atl_ind),plat(atl_ind),'Color',[102/255 0 0 ],'LineWidth',2)
title({sprintf('Depth= %.2f m',100)})
hold off

%% surface 3 --> Peak value surface
 
subaxis(2,2,3,'Spacing',s1,'Padding',p1,'Margin',m1,'PaddingBottom',2*pb,'PaddingRight',0,'PaddingTop',15*p1,'PaddingLeft',0);
m_proj('Lambert','lon',[-100 15], 'lat', [-10 65]) ;
set(gca,'Position',[0.01    0.10-0.01    0.4875    0.35]);
m_proj('Lambert','lon',[-100 20], 'lat', [-10 65]) 

timevarname='PeakVALUE';
varunit3=varunit;
fieldPeak=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/all_ramps/filtered/Results/New_time_var_results/%s_results_updated_v1D_poly1_NEW_NEW.mat',varlist{varIDX}),timevarname);
fieldPeak=fieldPeak.PeakVALUE;

axmin3=0;
axmax3=250;

%unflipping
bit1=fieldPeak(:,122:end); %bit1 (200:end,:,:); %reshaping -- position 200 is where Pacific is
bit2=fieldPeak(:,1:121);   %bit2 (1:199,:,:);
fieldPeak=cat(2,bit1,bit2);

field3=fieldPeak;

for smoothnumber=1:smoothnumbsurface %smoothing 7 times at the specific depth level
    field3=smooth2d(field3);
end

%surface
m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(field3*12.0107*1000*24*3600)); shading interp

colorscheme3=parula(10);
colormap(gca,colorscheme3)
caxis([axmin3 axmax3])

h=colorbar;
h.Location='westoutside';
h.Label.String=[{sprintf('Max_c_h_a_n_g_e %s',varname); sprintf('(%s)',varunit)}];
h.TickLabels(end)={['\geq ' num2str(axmax3)]};


m_coast('patch',[.7 .7 .7],'edgecolor','k');


hold on
%Grid box
m_grid('backgroundcolor','k','tickdir','out'); %if brown [160/255 82/255 45/255]
hold on
m_line(plon(atl_ind),plat(atl_ind),'Color',[102/255 0 0 ],'LineWidth',2)
title({sprintf('Depth= %.2f m',100)})
hold off



%% surface 4 - Time of recovery 
%surface
subaxis(2,2,4,'Spacing',s1,'Padding',p1,'Margin',m1,'PaddingBottom',2*pb,'PaddingRight',0,'PaddingTop',15*p1,'PaddingLeft',0);
m_proj('Lambert','lon',[-100 15], 'lat', [-10 65]) ;
set(gca,'Position',[0.48    0.10-0.01    0.4875    0.35]);
m_proj('Lambert','lon',[-100 20], 'lat', [-10 65]) 


axmin5=120;
axmax5=520;
timevarname='Trecovery';
fieldTrec=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/all_ramps/filtered/Results/New_time_var_results/%s_results_updated_v1D_poly1_NEW_NEW.mat',varlist{varIDX}),timevarname);

%Possible flags included in the Trecovery dataset
badflagsnegslope=-7777;
badflagtooearly=-4444;
bottomflag=-1111;


%Importing Trevocery dataset
%fieldZERO=Trecovery;
fieldZERO=fieldTrec.Trecovery{1,2};
fieldZERO=real(fieldZERO);
%fieldZERO(fieldZERO==badflagtooearly)=NaN; %this makes all different flags the same value %Gets

field1=permute(fieldZERO,[1,2]);

% figure
% pcolor(squeeze(field1(55,:,:)))
% shading flat
% colorbar
% caxis([160 560])


%unflipping
bit1=field1(:,122:end); %bit1 (200:end,:,:); %reshaping -- position 200 is where Pacific is
bit2=field1(:,1:121);   %bit2 (1:199,:,:);
field2=cat(2,bit1,bit2);

for l=1:384
    for c=1:320
    if fieldTOD2(l,c)<0;
    field2(l,c)=80; %this makes all different flags the same value %Definite no recovery
    end
    end
end


field2(field2<0 & field2~=bottomflag)=600; %this makes all negative flags which are not bottom positive to indicate very late or no recovery


field3=field2;
field3(field3==-1111)=NaN; %this makes all bottom NaN and removes it before the smoothing

for smoothnumber=1:smoothnumbsurface-1 %smoothing 3 times at the specific depth level
    field3=smooth2d(field3);
end

m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(field3)); shading interp
dummycolor=parula(10);

colormap(gca,[colornodept; colornodept; colornodept;dummycolor(1:8,:); repmat([255/255 140/255 0],1,1); repmat([255/255 99/255 71/255],1,1)]) % [177/255 47/255 245/255]
%colormap(gca,[[204/255 229/255 255/255];dummycolor(1:7,:); repmat([255/255 140/255 0],1,1); repmat([255/255 99/255 71/255],1,1)]) % [177/255 47/255 245/255]


caxis([0 axmax5])
m_coast('patch',[.7 .7 .7],'edgecolor','k');

%Grid box
m_grid('backgroundcolor','k','tickdir','out'); %if brown [160/255 82/255 45/255]
hold on
m_line(plon(atl_ind),plat(atl_ind),'Color',[102/255 0 0 ],'LineWidth',2)
title(sprintf('Depth= %.2f m',100))

h=colorbar;
h.Location='westoutside';
h.Label.String=sprintf('Trec %s (yr)',varname);
h.Ticks=160:40:axmax5;
%h.TickLabels(1)={'\leq 140'};
h.TickLabels(end)={'\geq 480'};

hold off
hFig.InvertHardcopy = 'off';
supersizeme(1.4)
print('/Volumes/LaCie_Leonardo/NorESM/export_production/surface_maps_ep100.png','-dpng','-r300')
print('/Volumes/LaCie_Leonardo/NorESM/export_production/surface_maps_ep100.pdf','-dpdf','-r300')


%%%%%%======================%%%%% %%%%%%======================%%%%%        

%% FIGURE

%% surface Percentage change + stipplings
% Loading Percentage change data for a specific period
hFig = figure(1); clf;
set(hFig, 'Units','centimeters','Position',[0 0  20 35  ]); %
set(hFig,'color','w','visible','on')
m1=0.01; p1=.002; s1=0.005; pb=.02; 



pairs={'pair1';'pair2';'pair3'};
letter={'a';'b';'c'};
axmin4=-25;
axmax4=25;
%PAIR 1 - End of mitigation
%PAIR 2 - Middle of extension
%PAIR 3 - End of extension
for pairIDX=1:3
    
if pairIDX==1
    period='Pre-Industrial x End of Mitigation phase';
elseif pairIDX==2
    period='Pre-Industrial x Middle of Extension phase';
elseif pairIDX==3
    period='Pre-Industrial x End of Extension phase';
end

%Data is shown in terms of Percentage_change_mean
    pchange=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/StatsResults/%s_stats_updated_v1D.mat',varlist{varIDX})...
        ,sprintf('Percentage_change_mean_%s',pairs{pairIDX}));
    pchange=eval(sprintf('pchange.Percentage_change_mean_%s',pairs{pairIDX}));

varunit4='%';

%unflipping
bit1=pchange(:,122:end); %bit1 (200:end,:,:); %reshaping -- position 200 is where Pacific is
bit2=pchange(:,1:121);   %bit2 (1:199,:,:);
pchange=cat(2,bit1,bit2);

%Loading stippling mask for surface plots

Pflags=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/StatsResults/%s_stats_updated_v1D.mat',varlist{varIDX})...
    ,sprintf('Pflags_%s',pairs{pairIDX}));
Pflags1=eval(sprintf('Pflags.Pflags_%s',pairs{pairIDX}));

% figure
% pcolor(Pflags1(:,:,1));shading flat
%
Ptest_ind=permute(Pflags1,[1,2]);
bit1=Ptest_ind(:,122:end); %bit1 (200:end,:,:); %reshaping -- position 200 is where Pacific is
bit2=Ptest_ind(:,1:121);   %bit2 (1:199,:,:);
Ptest_ind=cat(2,bit1,bit2);

% figure
% pcolor(Ptest_ind);shading flat
 
subaxis(3,1,pairIDX,'Spacing',s1,'Padding',p1,'Margin',m1,'PaddingBottom',2*pb,'PaddingRight',0,'PaddingTop',15*p1,'PaddingLeft',0);
m_proj('Lambert','lon',[-100 20], 'lat', [-10 65]) 

if pairIDX==1
a=get(gca,'Position');
set(gca,'Position',[a(1) a(2) 0.9800    0.2026]);
end
if pairIDX==2
set(gca,'Position',[a(1) a(2)-0.30 0.9800    0.2026]);
end
if pairIDX==3
set(gca,'Position',[a(1) a(2)-0.60 0.9800    0.2026]);
end

field3=pchange;

for smoothnumber=1:smoothnumbsurface %smoothing 7 times at the specific depth level
    field3=smooth2d(field3);
end

%surface
m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(field3)); shading interp
hold on

%drawing stippling
mask1=Ptest_ind(:,:)==100;%% these are the stippling masks for each depth level k
m_scatter(plon(mask1),plat(mask1),2,[0 0 0], '.')

%adjusting colormap
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


colormap(gca,redblueleo)
caxis([axmin4 axmax4])
m_coast('patch',[.7 .7 .7],'edgecolor','k');

if pairIDX==3
h=colorbar;
h.Location='southoutside';
h.Label.String=sprintf('%s Percentage change (%s) \n',varname,varunit4);
h.Label.FontSize=14;
h.Label.FontWeight='normal';
h.Ticks=[-25 -15 0 15 25]
h.TickLabels(end)={['\geq ' num2str(axmax4)]};
h.TickLabels(1)={['\leq ' num2str(axmin4)]};
h.FontSize=12;

end   

hold on
%Grid box
m_grid('backgroundcolor','k','tickdir','out','FontSize',14,'xtick',[-90 -60 -30 0], 'ytick',[60 30 0]); %if brown [160/255 82/255 45/255]
hold on
m_line(plon(atl_ind),plat(atl_ind),'Color',[102/255 0 0 ],'LineWidth',2)
title(sprintf('%s) %s \n Depth = %.2f m',letter{pairIDX},period,100),'Fontsize',14)
if pairIDX==3
set(gca,'Position',[a(1) a(2)-0.6 0.9800    0.2026]);
end

hold off
colormap(redblue(25))
end
print('/Volumes/LaCie_Leonardo/NorESM/export_production/surface_maps_ep100_percentage_change.png','-dpng','-r300')
print('/Volumes/LaCie_Leonardo/NorESM/export_production/surface_maps_ep100_percentage_change.pdf','-dpdf','-r300')


  



%% surface DELTA change + stipplings
% Loading Percentage change data for a specific period
hFig = figure(1); clf;
set(hFig, 'Units','centimeters','Position',[0 0  20.7433 37.0769  ]); %
set(hFig,'color','w','visible','on')
m1=0.01; p1=.002; s1=0.005; pb=.02; 



pairs={'pair1';'pair2';'pair3'};
letter={'a';'b';'c'};
axmin4=-25;
axmax4=25;
%PAIR 1 - End of mitigation
%PAIR 2 - Middle of extension
%PAIR 3 - End of extension
for pairIDX=1:3
    
if pairIDX==1
    period='Pre-Industrial x End of Mitigation phase';
elseif pairIDX==2
    period='Pre-Industrial x Middle of Extension phase';
elseif pairIDX==3
    period='Pre-Industrial x End of Extension phase';
end

%Data is shown in terms of Percentage_change_mean
    pchange=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/StatsResults/%s_stats_updated_v1D.mat',varlist{varIDX})...
        ,sprintf('delta_%s',pairs{pairIDX}));
    pchange=eval(sprintf('pchange.delta_%s',pairs{pairIDX}));
        
varunit5='mgC m^-^2 day^-^1';

%unflipping
bit1=pchange(:,122:end); %bit1 (200:end,:,:); %reshaping -- position 200 is where Pacific is
bit2=pchange(:,1:121);   %bit2 (1:199,:,:);
pchange=cat(2,bit1,bit2);

%Loading stippling mask for surface plots

Pflags=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/StatsResults/%s_stats_updated_v1D.mat',varlist{varIDX})...
    ,sprintf('Pflags_%s',pairs{pairIDX}));
Pflags1=eval(sprintf('Pflags.Pflags_%s',pairs{pairIDX}));

% figure
% pcolor(Pflags1(:,:,1));shading flat
%
Ptest_ind=permute(Pflags1,[1,2]);
bit1=Ptest_ind(:,122:end); %bit1 (200:end,:,:); %reshaping -- position 200 is where Pacific is
bit2=Ptest_ind(:,1:121);   %bit2 (1:199,:,:);
Ptest_ind=cat(2,bit1,bit2);

% figure
% pcolor(Ptest_ind);shading flat
 
subaxis(3,1,pairIDX,'Spacing',s1,'Padding',p1,'Margin',m1,'PaddingBottom',2*pb,'PaddingRight',0,'PaddingTop',15*p1,'PaddingLeft',0);
m_proj('Lambert','lon',[-100 15], 'lat', [15 65]) ;

if pairIDX==1
a=get(gca,'Position');
set(gca,'Position',[a(1) a(2) 0.9800    0.2026]);
end
if pairIDX==2
set(gca,'Position',[a(1) a(2)-0.30 0.9800    0.2026]);
end
if pairIDX==3
set(gca,'Position',[a(1) a(2)-0.60 0.9800    0.2026]);
end

field3=pchange;

for smoothnumber=1:smoothnumbsurface %smoothing 7 times at the specific depth level
    field3=smooth2d(field3);
end

%surface
m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(field3*12.0107*1000*24*3600)); shading interp
hold on

%drawing stippling
mask1=Ptest_ind(:,:)==100;%% these are the stippling masks for each depth level k
m_scatter(plon(mask1),plat(mask1),10,[0 0 0], '.')

%adjusting colormap
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


colormap(gca,redblueleo)
caxis([axmin4 axmax4])
m_coast('patch',[.7 .7 .7],'edgecolor','k');

if pairIDX==3
h=colorbar;
h.Location='southoutside';
h.Label.String=sprintf('\\Delta %s (%s)',varname,varunit5);
h.Label.FontSize=14;
h.Label.FontWeight='normal';
%h.TickLabels(end)={['\geq ' num2str(axmax4)]};
%h.TickLabels(1)={['\leq ' num2str(axmin4)]};
h.FontSize=12;

end   

hold on
%Grid box
m_grid('backgroundcolor','k','tickdir','out','FontSize',14); %if brown [160/255 82/255 45/255]
hold on
m_line(plon(atl_ind),plat(atl_ind),'Color',[102/255 0 0 ],'LineWidth',2)
title(sprintf('%s) %s \n Depth = %.2f m',letter{pairIDX},period,100),'Fontsize',14)
if pairIDX==3
set(gca,'Position',[a(1) a(2)-0.6 0.9800    0.2026]);
end

hold off
end
print('/Volumes/LaCie_Leonardo/NorESM/export_production/surface_maps_ep100_delta.png','-dpng','-r300')



%% time series methods poster

k=100
varlist={'epc100'};

 %% Importing Pre-Industrial 5-yr mean

    PI_mean=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/Pre_industrial/30yrMeans/PI_%s_first_5yr_mean.mat',varlist{varIDX}));
    PI_mean=PI_mean.PI_mean*12.0107*1000*24*3600; 

    
      PI_mean_bit1=PI_mean(200:end,:,:); %reshaping -- position 200 is where Pacific is
      PI_mean_bit2=PI_mean(1:199,:,:);
    
    PI_mean=[PI_mean_bit1 ; PI_mean_bit2];
    PI_mean=permute(PI_mean,[2,1,3]);      %tranposing variables from (320X384 to 384X320) to match domain_mask


    %% Importing Pre-Industrial 30-yr Standard Deviation

    
    PI_sd=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/Pre_industrial/30yrMeans/PI_%s_30yr_sd.mat',varlist{varIDX}));
    PI_sd=PI_sd.PI_sd*12.0107*1000*24*3600;
    
    PI_sd_bit1=PI_sd(200:end,:,:); %reshaping -- position 200 is where Pacific is
    PI_sd_bit2=PI_sd(1:199,:,:);
    PI_sd=[PI_sd_bit1 ; PI_sd_bit2];
    
    PI_sd=permute(PI_sd,[2,1,3]);      %tranposing variables from (320X384 to 384X320) to match domain_mask


   
   clear *bit*
   
   
varIDX=1;
load(sprintf('/Volumes/LaCie_Leonardo/NorESM/all_ramps/filtered/new_series_%s_k_%d.mat',varlist{varIDX},k))
l=309;c=112
series=new_series{l,c}*12.0107*1000*24*3600;


figure
plot(series)

series_complete=series;  
%checking time of departure 



%% TESTING SINGLE POINT


 % SINGLE POINT FIGURE

 pre_ind_gridcell=PI_mean(l,c);  
 pre_ind_sd=PI_sd(l,c);
 PI_lb=PI_mean(l,c)-2*pre_ind_sd;
 PI_ub=PI_mean(l,c)+2*pre_ind_sd;
 
 
f1=figure
box on
%drawing PI horizontal line
h1=yline(PI_mean(l,c))
hold on
%drawing PI horizontal standar deviation line
yline(PI_lb,'--')
hold on
yline(PI_ub,'--')
hold on

%drawing hatched areas
%upper
xx=1:700;
X=fliplr(xx);
y1=zeros(1,700)+PI_mean(l,c);
y2=zeros(1,700)+PI_ub;
h2=fill([xx X], [y1 fliplr(y2)], [.9 .9 .9], 'linestyle', 'none')
hold on
%lower
xx=1:700;
X=fliplr(xx);
y1=zeros(1,700)+PI_mean(l,c);
y2=zeros(1,700)+PI_lb;
h2=fill([xx X], [y1 fliplr(y2)], [.9 .9 .9], 'linestyle', 'none')
hold on

%drawing ramps
h3=plot(1:140,series_complete(1:140))
hold on
h4=plot(140:280,series_complete(140:280))
hold on
h5=plot(280:length(series_complete),series_complete(280:end))
hold on

grid on

%drawing vertical line at end of ramp-down
xline(140,':','Mitigation start','LabelHorizontalAlignment','left','LabelVerticalAlignment','middle','FontSize',9)
hold on

%finding and drawing Time when series leaves hatched area
[colTOD]=(find(series_complete<PI_lb | series_complete>PI_ub));
gap=9;
dummy=diff(colTOD)==1;
M = movsum(dummy,[0 gap]);
noreturn=find(M==10,1,'first');
if ~isempty(noreturn)
xline(colTOD(noreturn),':','ToD')
end
%finding and drawing Time of gloabl absolute maximum
% [valmin, col]=min(series_complete); 
% xline(col,':','T_l_o_w_e_s_t')
hold on
[valmax , colmax]=max(abs(series_complete-PI_mean(l,c)));
xline(colmax,':','T_p_e_a_k')
hold on

%finding and drawing Time when series gets back to hatched area
              
               
%              
%                 if isempty(point_envelope_single) %if no point enters the envelope AND last part of series is below envelope
%                
%                     %calculate regression based on post ramp points
%                     curve=fit(time_regress_single(end-30:end),pHseries_regress_single(end-50:end),'poly1');
%                     p1 = curve.p1;
%                     p2 = curve.p2;
%                     trecovery_dummy=round((PI_lb-p2)/p1);
%                
%                 %if regression time is higher than 3000 yrs or slope of regression is negative it means no recovery and grid cell receives bad flag 
%                     if  trecovery_dummy>3000 || p1<0
%                     trecovery(l,c)=badflag;
% 
%                     else %if there regression crosses envelope at some point ealrier than 3000y
%                     trecovery(l,c)=trecovery_dummy; 
%                     end
%                 
%                
%                 else %in case it enters the envelope take the index difference to find the first 10 consecutive points
%                      dummy2=diff(point_envelope_single)==1;
%                      M2 = movsum(dummy2,[0 gap]);
%                      noreturn2=find(M2==10,1,'first');
%                      
%                     trecovery_dummy=point_envelope_single(noreturn2)+mitigationstart-1;
%                     trecovery(l,c)=trecovery_dummy;  % if there is a point inside the envelope
%                     curve=[];
%                 end
%                      
                
   %% aaa
       mitigationstart=140;
       regwindow=100;

                time_regress_single=load('/Volumes/LaCie_Leonardo/NorESM/all_ramps/time_array.mat');
                time_regress_single=time_regress_single.time_lowest_regress;
                time_regress_single=time_regress_single+1; 
                
                %finding points after Tpeak to determine when series enters envelope            
                
                if series_complete(colmax)<PI_lb
                    
                    %PI_mean(l,c)-valmax<PI_lb
                
                point_envelope_single=find(series_complete(colmax:end)>PI_lb);

                 noreturn2=[];
                    trecovery_dummy=[];
                    curve=[]; p1=NaN;
                    
                    if ~isempty(point_envelope_single)  %in case there are points that enter the envelope take the index difference to find the first 10 consecutive points after the mitigation start
                        dummy2=diff(point_envelope_single)==1;
                        M2 = movsum(dummy2,[0 gap]);
                        noreturn2=find(M2==10,1,'first');
                        
                        if ~isempty(noreturn2) %if there are 10 consecutive points
                            trecovery_dummy=point_envelope_single(noreturn2)+mitigationstart-1;
                            trecovery(l,c)=trecovery_dummy;
                            
                            
                        else %if there are no 10 consecutive points
                            %calculate regression based on extension ramp points
                            curve=fit(time_regress_single(end-regwindow:end),series_complete(end-regwindow:end)','poly1');
                            p1 = curve.p1;
                            p2 = curve.p2;
                            trecovery_dummy=round((PI_lb-p2)/p1);
                            trecovery(l,c)=trecovery_dummy; % there regression crosses envelope at some point ealrier than 1000y
                        end
                        
                        
                    elseif isempty(point_envelope_single) || isempty(noreturn2) %if there are no points to enter envelope do a regression using last 50 yr of the series anyway
                        %calculate regression based on post ramp points
                        curve=fit(time_regress_single(end-regwindow:end),series_complete(end-regwindow:end)','poly1');
                        p1 = curve.p1;
                        p2 = curve.p2;
                        trecovery_dummy=round((PI_lb-p2)/p1);
                        trecovery(l,c)=trecovery_dummy; % there regression crosses envelope at some point ealrier than 1000y
                        
                    end
                end
                    
              %%      
                if series_complete(colmax)>PI_ub
                    
                    %PI_mean(l,c)+valmax>PI_ub
                    
                 point_envelope_single=find(series_complete(colmax:end)<PI_ub);

                    noreturn2=[];
                    trecovery_dummy=[];
                    curve=[]; p1=NaN;
                    
                    if ~isempty(point_envelope_single)  %in case there are points that enter the envelope take the index difference to find the first 10 consecutive points after the mitigation start
                        dummy2=diff(point_envelope_single)==1;
                        M2 = movsum(dummy2,[0 gap]);
                        noreturn2=find(M2==10,1,'first');
                        
                        if ~isempty(noreturn2) %if there are 10 consecutive points
                            trecovery_dummy=point_envelope_single(noreturn2)+mitigationstart-1;
                            trecovery(l,c)=trecovery_dummy;
                            
                            
                        else %if there are no 10 consecutive points
                            %calculate regression based on extension ramp points
                            curve=fit(time_regress_single(end-regwindow:end),series_complete(end-regwindow:end)','poly1');
                            p1 = curve.p1;
                            p2 = curve.p2;
                            trecovery_dummy=round((PI_ub-p2)/p1);
                            trecovery(l,c)=trecovery_dummy; % there regression crosses envelope at some point ealrier than 1000y
                        end
                        
                        
                    elseif isempty(point_envelope_single) || isempty(noreturn2) %if there are no points to enter envelope do a regression using last 50 yr of the series anyway
                        %calculate regression based on post ramp points
                        curve=fit(time_regress_single(end-regwindow:end),series_complete(end-regwindow:end)','poly1');
                        p1 = curve.p1;
                        p2 = curve.p2;
                        trecovery_dummy=round((PI_ub-p2)/p1);
                        trecovery(l,c)=trecovery_dummy; % there regression crosses envelope at some point ealrier than 1000y
                        
                    end
                 end
                    
                    
                    
                    
                    %% flagging
                    
                         
%         badflagslow=-9999;
%         badflagsnegslope=-7777;
%         badflagtooearly=-4444;
%         bottomflag=-1111;
%                     
%                     %if trecovery_dummy happens before mitigation start
%                     %then flag early
%                     if trecovery_dummy<mitigationstart+gap
%                         trecovery(l,c)=badflagtooearly;
%                     end
%                     
%                     %if regression time turns out to be higher
%                     %than 1000 yrs AND slope of regression is
%                     %positive it means very slow recovery
%                     if  trecovery_dummy>1000 && p1>0
%                         trecovery(l,c)=badflagslow;
%                     end
%                     
%                     %if regression curve has negative slope,
%                     %means it will never recover.
%                     if  trecovery_dummy<mitigationstart || p1<0
%                         trecovery(l,c)=badflagsnegslope;
%                     end
                    
         %% 
         
         
                
%continuing with the plot     
if  trecovery(l,c)>0
xline(trecovery_dummy,':','T_r_e_c')
end
hold on
if ~isempty(curve)
    %plot(
    plot(curve,[time_regress_single(end-regwindow:end) ; trecovery_dummy ], [series_complete(end-regwindow:end)' ; PI_lb] )
end   
hold on
ytickformat('%.2f')
legend('PI_m_e_a_n','', '', '','-2*STDev(PI)', 'RAMP-UP', 'RAMP-DOWN', 'EXTENSION','','','','','','Regression fit', 'Location', 'southeast') 

set(gcf,'color','w')
ylabel('POC export at 100m (mgC m^-^2 day^-^1)')
xlabel('Time (yr)')
xlim([0 500])
supersizeme(1.4)

%% AVI 
  workingDir = '/Volumes/LaCie_Leonardo/NorESM/export_production/';
  outputVideo = VideoWriter(['/Volumes/LaCie_Leonardo/NorESM/export_production\' 'Series_export_production.avi']); %name of the output file
  outputVideo.FrameRate = 1; %number of frames per seconds 480/15 approx 30 sec
  open(outputVideo)

  for ii = 1:3
   img = imread([workingDir sprintf('series_poc_export_%d.png',ii)]);
   writeVideo(outputVideo,img)
  end
  close(outputVideo)
  
  %%GIF
  filename = [workingDir 'POC_envelope.gif'];
   for ii = 1:3
   img = imread([workingDir sprintf('series_poc_export_%d.png',ii)]);
  % Capture the plot as an image 
      %frame = getframe(img); 
      %im = frame2im(frame); 
      [imind,cm] = rgb2ind(img,256); 
      % Write to the GIF File 
      if ii == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',1.5); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',1.5); 
      end 
  end
