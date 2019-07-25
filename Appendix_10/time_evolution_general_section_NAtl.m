depth=load('/Volumes/LaCie_Leonardo/NorESM/all_ramps/Depth_Levels.mat'); %get depth array from one of the files in the folder
depth=depth.depths;
depth=depth';


%% First plot
%% First plot is CO2 atm
co2start=284.7;
co2annual=zeros(480,1);

co2annual(1)=co2start;

for i=2:length(co2annual)
    if i<=140
    co2annual(i)=co2annual(i-1)*1.01;
    end
    
    if i>=141 && i<=280 && co2annual(i-1)>=co2start
    co2annual(i)=co2annual(i-1)*0.99;
    elseif i>=141 && i<=280 && co2annual(i-1)<=co2start
      co2annual(i)=co2start; 
    end

    if i>=281
    co2annual(i)=co2start;
    end
end

hFig = figure(1); clf;
set(hFig, 'Units','centimeters','Position',[0,0,40,15]); %
m1=0.02; p1=.002; s1=0.03; pb=.05; 

%% Panel 1
% %subaxis(2,3,1,'Spacing',s1,'Padding',p1,'Margin',m1,'PaddingBottom',2*pb,'PaddingRight',0,'PaddingTop',15*p1,'PaddingLeft',0);rampup_sst=plot(1:140,co2annual(1:140));
% subplot(2,3,1)
% rampup=plot(1:140,co2annual(1:140),'r','linewidth',8);
% hold on
% rampdown=plot(140:280,co2annual(140:280),'b','linewidth',8);
% hold on
% extension=plot(280:480,co2annual(280:480),'c','linewidth',8);
% hold on 
% xlabel('Time (yr)')
% ylabel('Atmospheric CO_2 (ppm)')
% grid on
% title('a) Anthropogenic forcing','FontSize',20)
% yline(co2start,'--', 'Pre-industrial baseline 284.7 ppm','LabelVerticalAlignment','bottom','FontWeight','bold','FontSize',16,'LineWidth',2)
% yline(max(co2annual),'--',['CO_2 peak ' sprintf('%.2f ppm',max(co2annual))],'LabelVerticalAlignment','bottom','FontWeight','bold','FontSize',16,'LineWidth',2)
% lgd=legend([rampup rampdown extension],'Ramp-up','Ramp-down','Extension')
% set(lgd,'Position',[0.195,0.80,0.1133,0.0576])
% set(gca,'FontName','Helvetica','fontsize',20,'TickDir','in','TickLength',[.005 .01],'XMinorTick','on','YMinorTick','on',...
%     'xtick',0:100:500,'LineWidth',2)
% set(gca,'Position',[0.13,0.5838,0.18,0.3412]);
% grid on



%SST + MEAN OCEAN TEMPERATURE
temp=load('/Volumes/LaCie_Leonardo/NorESM/Initial_figs/area_and_vol_NAtl_templvl.mat');
subplot(1,3,1)
[AX,H1,H2]=plotyy(1:481,temp.seriesSURF,1:481,temp.ocn_volseries);
set(H1,'linewidth',2)
set(H2,'linewidth',2)
title('b)Temperature','FontSize',20)
set(AX(:),'xlim',([0 480]))
xlabel('Time (yr)')
ylabel(AX(1),'SST (\circC)')
ylabel(AX(2),'NAtl T (\circC)')
set(AX(:),'FontName','Helvetica','fontsize',20,'TickDir','in','TickLength',[.005 .01],'XMinorTick','on','YMinorTick','on',...
    'xtick',0:100:500,'LineWidth',2)
title('a) Temperature','FontSize',20)
set(gca,'Position',[0.05,0.20,0.20,0.70]);
grid on


%SURFACE PH + MEAN OCEAN PH
ph=load('/Volumes/LaCie_Leonardo/NorESM/Initial_figs/area_and_vol_NAtl_ph.mat');
subplot(1,3,2)
[AX,H1,H2]=plotyy(1:481,ph.seriesSURF,1:481,ph.ocn_volseries);
set(AX(:),'xlim',([0 480]),'ylim',([7.6 8.2]))
set(AX(:),'YTick',[7.6 7.8 8.0 8.2])
set(H1,'linewidth',2)
set(H2,'linewidth',2)
xlabel('Time (yr)')
ylabel(AX(1),'Surface pH')
ylabel(AX(2),'NAtl pH')
set(AX(:),'FontName','Helvetica','fontsize',20,'TickDir','in','TickLength',[.005 .01],'XMinorTick','on','YMinorTick','on',...
    'xtick',0:100:500,'LineWidth',2)
set(gca,'Position',[0.1+0.28,0.2,0.20,0.70]);
grid on
title('b) pH','FontSize',20)



%OXYGEN+AOU
%subaxis(2,3,2,'Spacing',s1,'Padding',p1,'Margin',m1,'PaddingBottom',2*pb,'PaddingRight',0,'PaddingTop',15*p1,'PaddingLeft',0);rampup_sst=plot(1:140,co2annual(1:140));
subplot(1,3,3)
Oxy=load('/Volumes/LaCie_Leonardo/NorESM/Initial_figs/area_and_vol_NAtl_o2.mat');
Aou=load('/Volumes/LaCie_Leonardo/NorESM/Initial_figs/area_and_vol_NAtl_AOU.mat');
[AX,H1,H2]=plotyy(1:481,((Oxy.ocn_volseries/1024)*1000000),1:481,Aou.ocn_volseries);
set(AX(:),'xlim',([0 480]))
set(H1,'linewidth',2)
set(H2,'linewidth',2)
xlabel(AX(1),'Time (yr)')
ylabel(AX(1),'DO (\mumol O_2 kg^-^1)')
ylabel(AX(2),'AOU (\mumol O_2 kg^-^1)')
set(AX(:),'FontName','Helvetica','fontsize',20,'TickDir','in','TickLength',[.005 .01],'XMinorTick','on','YMinorTick','on',...
    'xtick',0:100:500,'LineWidth',2)
title('c) DO vs AOU','FontSize',20)
set(gca,'Position',[0.1+0.62 0.20 0.20 0.70]);
grid on

set(gcf,'color','w')
print('/Volumes/LaCie_Leonardo/NorESM/Initial_figs/fig1_NAtl.png','-dpng','-r300')

%=====================================================================================================================================================================
%=====================================================================================================================================================================
%=====================================================================================================================================================================

% %OMEGA
% %subaxis(2,3,3,'Spacing',s1,'Padding',p1,'Margin',m1,'PaddingBottom',2*pb,'PaddingRight',0,'PaddingTop',15*p1,'PaddingLeft',0);rampup_sst=plot(1:140,co2annual(1:140));
% subplot(2,3,4:6)
% omg=load('/Volumes/LaCie_Leonardo/NorESM/Initial_figs/area_and_vol_NAtl_omegac.mat');
% addpath /Volumes/LaCie_Leonardo/NorESM/all_ramps
% panel3=pcolor(1:481,-depth,omg.ocn_volseries_depth);shading interp
% 
% caxis([0 2])
% colorscheme=redblue;
% colorscheme2=[colorscheme(1:50,:);colorscheme(71:end,:)];
% colorscheme=flip(colorscheme2,2);
% colormap(colorscheme)
% 
% h=colorbar;
% h.Label.String='\Omega_c';
% %h.TickLabels(end)={['\geq' num2str(max(caxis))]};
% h.Location='eastoutside';
% h.FontSize=18;
% 
% xlim([0 480])
% ylim([-3000 0])
% xlabel('Time (yr)')
% ylabel('Depth (m)')
% title('e) Evolution of NAtl mean \Omega_c','FontSize',14)
% set(gca,'FontName','Helvetica','fontsize',18,'TickDir','in','TickLength',[.008 .01],'XMinorTick','on','YMinorTick','on',...
%     'xtick',0:100:500,'XGrid','on','Ygrid','on','layer','top','LineWidth',2)
% %set(gca,'Position',[0.45,0.11,0.42,0.3412]);
% set(gcf,'color','w')
% 
% print('/Volumes/LaCie_Leonardo/NorESM/Initial_figs/fig1_NAtl.png','-dpng','-r300')


%%
%=====================================================================================================================================================================
%=====================================================================================================================================================================
%=====================================================================================================================================================================

fig=figure('color', 'white','units','centimeters', 'Position', [0 0 13 32])

%Temperature section
%subaxis(2,3,3,'Spacing',s1,'Padding',p1,'Margin',m1,'PaddingBottom',2*pb,'PaddingRight',0,'PaddingTop',15*p1,'PaddingLeft',0);rampup_sst=plot(1:140,co2annual(1:140));
subplot(5,1,1)
omg=load('/Volumes/LaCie_Leonardo/NorESM/Initial_figs/area_and_vol_templvl.mat');
addpath /Volumes/LaCie_Leonardo/NorESM/all_ramps
panel1=pcolor(1:481,-depth,omg.ocn_volseries_depth);shading interp

caxis([0 10])
h=colorbar;
h.Label.String='T (\circC)';
h.FontSize=12;
h.Ticks=[0 5 10];
h.Location='eastoutside';
colormap(gca,parula(10))
xlim([0 480])
ylim([-5000 0])
%xlabel('Time (yr)')
ylabel('Depth (m)')
title('a) Evolution of NAtl mean T','FontSize',20)
set(gca,'FontName','Helvetica','fontsize',12,'TickDir','in','TickLength',[.008 .01],'XMinorTick','on','YMinorTick','on',...
    'xtick',0:100:500,'XGrid','on','Ygrid','on','layer','top','LineWidth',2,'ytick',-5000:1000:0)% ...
    %'Position',[    0.1501    0.85    0.6364    0.1137]);
hold on
l1=xline(140,':','CO_2 peak','LabelHorizontalAlignment','left','FontWeight','bold');
l2=xline(280,':','End of mitigation','LabelHorizontalAlignment','left','FontWeight','bold');

%%
%pH section
%subaxis(2,3,3,'Spacing',s1,'Padding',p1,'Margin',m1,'PaddingBottom',2*pb,'PaddingRight',0,'PaddingTop',15*p1,'PaddingLeft',0);rampup_sst=plot(1:140,co2annual(1:140));
subplot(5,1,2)
omg=load('/Volumes/LaCie_Leonardo/NorESM/Initial_figs/area_and_vol_ph.mat');
addpath /Volumes/LaCie_Leonardo/NorESM/all_ramps
panel2=pcolor(1:481,-depth,omg.ocn_volseries_depth);shading interp

caxis([7.2 8.2])
h=colorbar;
h.Label.String='pH';
h.FontSize=12;
h.Ticks=[7.2 7.7 8.2];

%h.TickLabels(end)={['\geq' num2str(max(caxis))]};
h.Location='eastoutside';
colormap(gca,parula(10));
xlim([0 480])
ylim([-5000 0])

%xlabel('Time (yr)')
ylabel('Depth (m)')
title('b) Evolution of NAtl mean pH','FontSize',20)
set(gca,'FontName','Helvetica','fontsize',12,'TickDir','in','TickLength',[.008 .01],'XMinorTick','on','YMinorTick','on',...
    'xtick',0:100:500,'XGrid','on','Ygrid','on','layer','top','LineWidth',2,'ytick',-5000:1000:0, 'XTickLabel',0:100:500)%,...
   % 'Position',[    0.1501    0.68    0.6364    0.1137])
hold on
l1=xline(140,':','CO_2 peak','LabelHorizontalAlignment','left','FontWeight','bold');
l2=xline(280,':','End of mitigation','LabelHorizontalAlignment','left','FontWeight','bold');

%%
%O2 section
%subaxis(2,3,3,'Spacing',s1,'Padding',p1,'Margin',m1,'PaddingBottom',2*pb,'PaddingRight',0,'PaddingTop',15*p1,'PaddingLeft',0);rampup_sst=plot(1:140,co2annual(1:140));
subplot(5,1,3)
omg=load('/Volumes/LaCie_Leonardo/NorESM/Initial_figs/area_and_vol_o2.mat');
addpath /Volumes/LaCie_Leonardo/NorESM/all_ramps
panel3=pcolor(1:481,-depth,(omg.ocn_volseries_depth/1024)*1000000);shading interp

caxis([100 300])
h=colorbar;
h.Label.String='DO (\mumol O_2 kg^-^1)';
h.FontSize=12;

%h.TickLabels(end)={['\geq' num2str(max(caxis))]};
h.Location='eastoutside';
h.Ticks=100:40:300;
colormap(gca,parula(10));
xlim([0 480])
ylim([-5000 0])

%xlabel('Time (yr)')
ylabel('Depth (m)')
title('c) Evolution of NAtl mean DO','FontSize',20)
set(gca,'FontName','Helvetica','fontsize',12,'TickDir','in','TickLength',[.008 .01],'XMinorTick','on','YMinorTick','on',...
    'xtick',0:100:500,'XGrid','on','Ygrid','on','layer','top','LineWidth',2,'ytick',-5000:1000:0, 'XTickLabel',0:100:500)%,...
    %'Position',[    0.1501    0.50    0.6364    0.1137])
hold on
l1=xline(140,':','CO_2 peak','LabelHorizontalAlignment','left','FontWeight','bold');
l2=xline(280,':','End of mitigation','LabelHorizontalAlignment','left','FontWeight','bold');

%%
%AOU section
%subaxis(2,3,3,'Spacing',s1,'Padding',p1,'Margin',m1,'PaddingBottom',2*pb,'PaddingRight',0,'PaddingTop',15*p1,'PaddingLeft',0);rampup_sst=plot(1:140,co2annual(1:140));
subplot(5,1,4)
omg=load('/Volumes/LaCie_Leonardo/NorESM/Initial_figs/area_and_vol_AOU.mat');
addpath /Volumes/LaCie_Leonardo/NorESM/all_ramps
panel4=pcolor(1:481,-depth,omg.ocn_volseries_depth);shading interp

caxis([0 150])
h=colorbar;
h.Label.String='AOU (\mumol O_2 kg^-^1)';
h.FontSize=12;
%h.TickLabels(end)={['\geq' num2str(max(caxis))]};
h.Location='eastoutside';
colormap(gca,parula(15))
xlim([0 480])
ylim([-5000 0])

%xlabel('Time (yr)')
ylabel('Depth (m)')
title('d) Evolution of NAtl mean AOU','FontSize',20)
set(gca,'FontName','Helvetica','fontsize',12,'TickDir','in','TickLength',[.008 .01],'XMinorTick','on','YMinorTick','on',...
    'xtick',0:100:500,'XGrid','on','Ygrid','on','layer','top','LineWidth',2,'ytick',-5000:1000:0,'XTickLabel',0:100:500)
%set(gca,'Position',[0.45,0.11,0.42,0.3412]);
hold on 
l1=xline(140,':','CO_2 peak','LabelHorizontalAlignment','left','FontWeight','bold');
l2=xline(280,':','End of mitigation','LabelHorizontalAlignment','left','FontWeight','bold');


%OMEGA
%subaxis(2,3,3,'Spacing',s1,'Padding',p1,'Margin',m1,'PaddingBottom',2*pb,'PaddingRight',0,'PaddingTop',15*p1,'PaddingLeft',0);rampup_sst=plot(1:140,co2annual(1:140));
subplot(5,1,5)
omg=load('/Volumes/LaCie_Leonardo/NorESM/Initial_figs/area_and_vol_NAtl_omegac.mat');
addpath /Volumes/LaCie_Leonardo/NorESM/all_ramps
panel5=pcolor(1:481,-depth,omg.ocn_volseries_depth);shading interp

caxis([0 2])
colorscheme=redblue;
colorscheme2=[colorscheme(1:50,:);colorscheme(71:end,:)];
colorscheme=flip(colorscheme2,2);
colormap(gca,colorscheme)

h=colorbar;
h.Label.String='\Omega_c';
%h.TickLabels(end)={['\geq' num2str(max(caxis))]};
h.Location='eastoutside';
h.FontSize=12;

xlim([0 480])
ylim([-5000 0])
xlabel('Time (yr)')
ylabel('Depth (m)')
title('e) Evolution of NAtl mean \Omega_c','FontSize',20)
set(gca,'FontName','Helvetica','fontsize',12,'TickDir','in','TickLength',[.008 .01],'XMinorTick','on','YMinorTick','on',...
    'xtick',0:100:500,'XGrid','on','Ygrid','on','layer','top','LineWidth',2,'ytick',-5000:1000:0)
%set(gca,'Position',[0.45,0.11,0.42,0.3412]);
hold on
l1=xline(140,':','CO_2 peak','LabelHorizontalAlignment','left','FontWeight','bold');
l2=xline(280,':','End of mitigation','LabelHorizontalAlignment','left','FontWeight','bold');

set(gcf,'color','w')
print('/Volumes/LaCie_Leonardo/NorESM/Initial_figs/fig2_NAtl.png','-dpng','-r300')


  
