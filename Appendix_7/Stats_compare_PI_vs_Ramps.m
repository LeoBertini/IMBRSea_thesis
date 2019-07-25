%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%STATISTICAL ANALYSES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% - Detrends last 30yr of Extension time of Biogeochem field
% - Vertical interpolation in 200m intervals
% - Performs AD normality tests 
% - Performs paired t-tests between detrended PI and Ext
% - Computes percentage change of Biogeochem fields
% - Saves stippling masks
% - Saves table of relative number of gird points with significant differences per depth stratum

% Author: Leonardo Bertini
% Modified: 22nd April 2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% CODE
clear all 
close all
warning('off','all')

%loading pre industrial mean and sd

 %% loading depths
   
   depth=load('/Volumes/LaCie_Leonardo/NorESM/all_ramps/Depth_Levels.mat'); %get depth array from one of the files in the folder
    depth=depth.depths;
   
   
   %% imporitng grid file lat and long
   
   longitude=ncread('/Volumes/LaCie_Leonardo/NorESM/all_ramps/grid_gx1v6.nc','plon') ;
   latitude=ncread('/Volumes/LaCie_Leonardo/NorESM/all_ramps/grid_gx1v6.nc','plat') ;

   %RESHAPING SO THAT LAT AND LONG GRIDS HAVE THE ATLANTIC IN THE MIDDLE
   longitude_bit1=longitude(200:end,:,:);
   longitude_bit2=longitude(1:199,:,:);
   longitude=[longitude_bit1; longitude_bit2];
   
   latitude_bit1=latitude(200:end,:,:);
   latitude_bit2=latitude(1:199,:,:);
   latitude=[latitude_bit1; latitude_bit2];

varlist={'ph','o2','omegac','detoc','templvl','AOU','o2sat'}; 
for varIDX=7%length(varlist)
   
    %% Importing Pre-Industrial 30-yr Standard Deviation
 
    PI_sd=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/Pre_industrial/30yrMeans/PI_%s_30yr_sd.mat',varlist{varIDX}));
    PI_sd=PI_sd.PI_sd;
    
    PI_sd_bit1=PI_sd(200:end,:,:); %reshaping -- position 200 is where Pacific is
    PI_sd_bit2=PI_sd(1:199,:,:);
    PI_sd=[PI_sd_bit1 ; PI_sd_bit2];
    
    %% Importing Pre-Industrial 30-yr Detrended series
    
    PI_detrend=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/Pre_industrial/30yrMeans/PI_%s_30yr_detrended.mat',varlist{varIDX}));
    PI_detrend=PI_detrend.PI_detrended;
    
    
    for l=1:size(PI_detrend,1) 
        for c=1:size(PI_detrend,2)
            for k=1:length(depth)
                if sum(isnan(PI_detrend{l,c,k}))==length(PI_detrend{l,c,k}) %whenever series have 30 invalid datapoints, assing a single NAN to that cell element 
                   PI_detrend{l,c,k}=NaN; 
                end
            end
        end
    end
    
PI_detrend_bit1=PI_detrend(200:end,:,:); %reshaping -- position 200 is where Pacific is
PI_detrend_bit2=PI_detrend(1:199,:,:);
PI_detrend=[PI_detrend_bit1 ; PI_detrend_bit2];

sprintf('%s PI sd and PI detrend imported',varlist{varIDX})
sprintf('\n')

clear *bit*

%swapping dimensions
longitude=permute(longitude,[2,1]);
latitude=permute(latitude,[2,1]);

PI_sd=permute(PI_sd,[2,1,3]);
PI_detrend=permute(PI_detrend,[2,1,3]);

%% Interpolation to make regular grid in the vertical for PI cell array
depth_interval=100;
tic
[dummystruc,PI_detrend_interp_series,~]=interp_vertical_cell(PI_detrend,depth,depth_interval);
toc
%depth_interp=min(depth):100:max(depth);

%% Interpolation to make regular grid in the vertical for PR individual files
cellarraypath='/Volumes/LaCie_Leonardo/NorESM/all_ramps/filtered/'; %this is where depth-filtered files are
varname=varlist{varIDX};
yearwindow=30;

tic
[dummystruc2,PR_interp_series,depth_interp]=interp_vertical_individual_files(cellarraypath,varname,yearwindow,depth,depth_interval);
toc

%%
%%depth_interp=depth;

%% STATS
sprintf('%s Initializing stats',varlist{varIDX})
sprintf('\n')

%Declaring parameters for stats   
alpha=0.05;
groupnames1=repmat({'Pre Industrial'},1,30); %Assingning categories for the datapoints to be able to do anova with 2 separate groups
groupnames2=repmat({'Mitigation End'},1,30); %Assingning categories for the datapoints to be able to do anova with 2 separate groups
groupnames3=repmat({'Extension Middle'},1,30); %Assingning categories for the datapoints to be able to do anova with 2 separate groups
groupnames4=repmat({'Extension End'},1,30);
groupnames=[groupnames1 , groupnames2, groupnames3, groupnames4]; %when doing anove, ALWAYS PLACE PI first


%Preallocating space
Ptest_individualgrid=zeros(size(latitude,1),size(latitude,2),size(depth_interp,1)) ; %flag for doimain being continent or bottom
Percentage_change_mean=zeros(size(latitude,1),size(latitude,2),size(depth_interp,1)) ;
Percentage_change_hidrogen_ion=zeros(size(latitude,1),size(latitude,2),size(depth_interp,1)) ;
Variance_Percentage_change=zeros(size(latitude,1),size(latitude,2),size(depth_interp,1)) ;

Elem_sgnif=zeros(size(depth_interp,1),1) ;
Elem_notsgnif=zeros(size(depth_interp,1),1) ;
Elem_notnormal=zeros(size(depth_interp,1),1) ;

%% Statistics on Vertically Interpolated fields (Now the PI and PR fields are spaced evenly in the vertical) 
%         load PI_interp_series
%         load PR_interp_series
tic
for k=1:size(depth_interp,1)
%     %loading filtered series and
%     %extracting last 30yrs for each grid-cell at each depth
%     load(sprintf('/Volumes/LaCie_Leonardo/NorESM/all_ramps/filtered/new_series_%s_k_%d.mat',varlist{varIDX},k)) 
    sprintf('Var = %s depth level %d\n',varlist{varIDX},k)

    for l=1:size(PR_interp_series,1)
        for c=1:size(PR_interp_series,2)
            
            if (sum(~isnan(PR_interp_series{l,c,k}(:)))==90) %if the cell actually contains a time series of 90 points (3 windows of 30)
                
                %PRE INDUSTRIAL
                PI_last30=PI_detrend_interp_series{l,c,k};            %selecting last 30years of Pre_industrial
                PI_last30_mean=nanmean(PI_detrend_interp_series{l,c,k}); %Pre-Industrial 30yr Mean
                
                %End of mitigation
                EndMitigation30=PR_interp_series{l,c,k}(1:30);
                dummy3=detrend(EndMitigation30);                   %removing trend
                EndMitigation30_mean=mean(EndMitigation30);        %Mitigation end  30yr mean
                EndMitigation30=dummy3 + EndMitigation30_mean;
                
                %Extension-Middle
                Ext_middle30=PR_interp_series{l,c,k}(31:60);    %selecting last 30 years of extension time
                dummy2=detrend(Ext_middle30);             %removing trend
                Ext_middle30_mean=mean(Ext_middle30);        %Extension middle  30yr mean
                Ext_middle30=dummy2 + Ext_middle30_mean;     %adding increments without thend to the mean value to get detrended series
                
                %Extension-END
                Ext_last30=PR_interp_series{l,c,k}(end-29:end);    %selecting last 30 years of extension time
                dummy1=detrend(Ext_last30);             %removing trend
                Ext_last30_mean=mean(Ext_last30);        %Extension time 30yr mean
                Ext_last30=dummy1 + Ext_last30_mean;     %adding increments without thend to the mean value to get detrended series
                
                % MEAN PERCENTAGE CHANGES PAIRWISE
                Percentage_change_mean_pair1(l,c,k)=percentagechangeNorESM(PI_last30_mean,EndMitigation30_mean);
                Percentage_change_mean_pair2(l,c,k)=percentagechangeNorESM(PI_last30_mean,Ext_middle30_mean);
                Percentage_change_mean_pair3(l,c,k)=percentagechangeNorESM(PI_last30_mean,Ext_last30_mean);
                
                %DELTAS
                if ~isnan(PI_last30_mean)
                delta_pair1(l,c,k)=EndMitigation30_mean-PI_last30_mean;
                delta_pair2(l,c,k)=Ext_middle30_mean-PI_last30_mean;
                delta_pair3(l,c,k)=Ext_last30_mean-PI_last30_mean;
                else
                delta_pair1(l,c,k)=NaN;
                delta_pair2(l,c,k)=NaN;
                delta_pair3(l,c,k)=NaN;
                end
                
                
                if varIDX==1 %if it's ph also do percentage change in hidrogen ion
                    Percentage_change_hidrogen_ion_pair1(l,c,k)= (((10^(-EndMitigation30_mean))-(10^(-PI_last30_mean)))/(10^(-PI_last30_mean)))*100; %calculating Percentage change of hidrogen ion concentration
                    Percentage_change_hidrogen_ion_pair2(l,c,k)= (((10^(-Ext_middle30_mean))-(10^(-PI_last30_mean)))/(10^(-PI_last30_mean)))*100; %calculating Percentage change of hidrogen ion concentration
                    Percentage_change_hidrogen_ion_pair3(l,c,k)= (((10^(-Ext_last30_mean))-(10^(-PI_last30_mean)))/(10^(-PI_last30_mean)))*100; %calculating Percentage change of hidrogen ion concentration
                
                else
                Percentage_change_hidrogen_ion_pair1(l,c,k)=NaN;
                Percentage_change_hidrogen_ion_pair2(l,c,k)=NaN;
                Percentage_change_hidrogen_ion_pair3(l,c,k)=NaN;
                end
                
                
                % % % % VARIANCE PERCENTAGE CHANGE
                % % % sd_preIND=nanstd(PI_last30);
                % % % sd_extension=nanstd(Ext_last30);
                % % % Variance_Percentage_change(l,c,k)= (((sd_extension^2)-(sd_preIND^2))/(sd_preIND^2))*100;
                % % %
                % % % figure
                % % % plot(PI_last30)
                % % % hold on
                % % % plot(PR_last30)
                % % % legend('Pre Ind','Extention')
                
                
                %% Testing for normality and then doing a 1-way anova to determine which differences are significant
                Pflags_pair1(l,c,k) = ttestNorESM(PI_last30,EndMitigation30); %%PI versus Extension END
                Pflags_pair2(l,c,k) = ttestNorESM(PI_last30,Ext_middle30); %%PI versus Extension END
                Pflags_pair3(l,c,k) = ttestNorESM(PI_last30,Ext_last30); %%PI versus Extension END

                %                             figure
                %                             boxplot([PI_last30 EndMitigation30 Ext_middle30 Ext_last30],groupnames)
                %                             title (sprintf('line %d col %d depth %d',l,c,k))
                %                             xlabel('Model Phase')
                %                             ylabel('pH')
                
                else %if does not contain area of atlantic ocean or falls onto continents or bottom
                
                Percentage_change_mean_pair1(l,c,k)=NaN;   Pflags_pair1(l,c,k)=NaN;
                Percentage_change_mean_pair2(l,c,k)=NaN;   Pflags_pair2(l,c,k)=NaN;
                Percentage_change_mean_pair3(l,c,k)=NaN;   Pflags_pair3(l,c,k)=NaN;
                
                Percentage_change_hidrogen_ion_pair1(l,c,k)=NaN;
                Percentage_change_hidrogen_ion_pair2(l,c,k)=NaN;
                Percentage_change_hidrogen_ion_pair3(l,c,k)=NaN;
                
                delta_pair1(l,c,k)=NaN;
                delta_pair2(l,c,k)=NaN;
                delta_pair3(l,c,k)=NaN;
 
            end  
        end
    end
    
%     pcolor(Ptest_individualgrid(:,:,40))
%     shading flat
%     colorbar
%     
    %% Percentages over whole domain for each depth level
    
    % Points with significant difference
    Elem_sgnif_pair1(k,1)= numel(find(Pflags_pair1(:,:,k)==100))/numel(find(~isnan(Pflags_pair1(:,:,k))));
    Elem_sgnif_pair2(k,1)= numel(find(Pflags_pair2(:,:,k)==100))/numel(find(~isnan(Pflags_pair2(:,:,k))));
    Elem_sgnif_pair3(k,1)= numel(find(Pflags_pair3(:,:,k)==100))/numel(find(~isnan(Pflags_pair3(:,:,k))));
    
   
    % Points with no significant difference
    Elem_notsgnif_pair1(k,1)=numel(find(Pflags_pair1(:,:,k)==-100))/numel(find(~isnan(Pflags_pair1(:,:,k))));
    Elem_notsgnif_pair2(k,1)=numel(find(Pflags_pair2(:,:,k)==-100))/numel(find(~isnan(Pflags_pair2(:,:,k))));
    Elem_notsgnif_pair3(k,1)=numel(find(Pflags_pair3(:,:,k)==-100))/numel(find(~isnan(Pflags_pair3(:,:,k))));

    % Points discarded because normality test showed series were not normally distributed 
    Elem_notnormal_pair1(k,1)=numel(find(Pflags_pair1(:,:,k)==-200))/numel(find(~isnan(Pflags_pair1(:,:,k))));
    Elem_notnormal_pair2(k,1)=numel(find(Pflags_pair2(:,:,k)==-200))/numel(find(~isnan(Pflags_pair2(:,:,k))));
    Elem_notnormal_pair3(k,1)=numel(find(Pflags_pair3(:,:,k)==-200))/numel(find(~isnan(Pflags_pair3(:,:,k))));
end

TABLE_STATS_pair1=table(Elem_sgnif_pair1,Elem_notsgnif_pair1,Elem_notnormal_pair1,depth_interp); %table with depth strata and their respective fractions
TABLE_STATS_pair2=table(Elem_sgnif_pair2,Elem_notsgnif_pair2,Elem_notnormal_pair2,depth_interp); %table with depth strata and their respective fractions
TABLE_STATS_pair3=table(Elem_sgnif_pair3,Elem_notsgnif_pair3,Elem_notnormal_pair3,depth_interp); %table with depth strata and their respective fractions

sprintf('Saving all results')
save (sprintf('/Volumes/LaCie_Leonardo/NorESM/StatsResults/%s_stats_updated_v1D_vertinterp.mat',varlist{varIDX}),...
    'PI_detrend_interp_series', 'PR_interp_series',...
    'delta_pair1','delta_pair2','delta_pair3',...
    'Percentage_change_mean_pair1', 'Percentage_change_mean_pair2','Percentage_change_mean_pair3',...
    'Percentage_change_hidrogen_ion_pair1','Percentage_change_hidrogen_ion_pair2','Percentage_change_hidrogen_ion_pair3',...
    'Pflags_pair1','Pflags_pair2','Pflags_pair3',...
    'TABLE_STATS_pair1', 'TABLE_STATS_pair2', 'TABLE_STATS_pair3')
toc
end


% figure
% pcolor(delta_pair3(:,:,6));shading flat; colorbar
% caxis([-0.25 0.25 ])


% % %% adjusting colormap
% % redblueleo=zeros(10,3);
% % %shades of red
% % redblueleo(end,:)=[1,0.0514285714285714,0.0514285714285714];
% % redblueleo(end-1,:)=[1,0.214285714285714,0.214285714285714];
% % redblueleo(end-2,:)=[1,0.357142857142857,0.357142857142857];
% % redblueleo(end-3,:)=[1,0.500000000000000,0.500000000000000];
% % redblueleo(end-4,:)=[1,0.642857142857143,0.642857142857143];
% % 
% % %shades of blue
% % redblueleo(1,:)=[0.0714285714285714,0.0714285714285714,1];
% % redblueleo(2,:)=[0.214285714285714,0.214285714285714,1];
% % redblueleo(3,:)=[0.357142857142857,0.357142857142857,1];
% % redblueleo(4,:)=[0.500000000000000,0.500000000000000,1];
% % redblueleo(5,:)=[0.642857142857143,0.642857142857143,1];
% % 
% % 
% % 
% % %% Making GIF images
% % 
% % % lat = linspace(min(latitude,[],'all'), max(latitude,[],'all'), size(latitude,2)); % ... and your real lat and lon here
% % % lon = linspace(min(longitude,[],'all'), max(longitude,[],'all'), size(longitude,1));
% % % [LT, LN] = meshgrid(lat, lon);
% % 
% % 
% % addpath /Volumes/LaCie_Leonardo/NorESM/scripts_jerry
% % 
% % 
% % 
% % for k = 1 : length(depth_interp)
% % 
% % fig=figure;
% % m_proj('Robinson','lon',[-80 30],'lat',[-5 70])
% % 
% % fielddummy=squeeze(Percentage_change_mean(:,:,k));
% % 
% % for smoothnumber=1:7
% %     fielddummy=smooth2d(fielddummy);
% % end
% % 
% % m_pcolor(longitude,latitude,fielddummy); shading interp
% % 
% % %stippling
% % mask=Ptest_individualgrid(:,:,k)==100;%%this is the stippling mask
% % hold on
% % m_scatter(longitude(mask),latitude(mask),2,[.5 .5 .5])
% % 
% % set(gca, 'color', 'k');
% % m_coast('patch',[.7 .7 .7],'edgecolor','k');
% % m_grid('box','fancy','tickdir','out');
% % colormap(redblueleo);
% %  h=colorbar('northoutside');
% %  h.Ticks=-2:0.4:2;
% %  title(h,'Percentage Change of Mean pH (%)','fontsize',12);
% %  caxis([-2  2]) 
% % h.TickLabels={'\leq-2';'-1.6';'-1.2';'-0.8';'-0.4';'0';'0.4';'0.8';'1.2';'1.6';'\geq2'};
% % patch([0.92 0.92 0.4 0.4],[0.10 -0.07 -0.07 0.10],'w','LineWidth',2); 
% % anotation=sprintf('Depth = %.2f m',depth_interp(k));
% % text(0.65, 0.02 ,anotation,'fontsize',10,'FontWeight','bold','color','k', 'vertical','middle','horizontal','center');
% % 
% %  drawnow
% %     frame = getframe(fig);
% %     im{k} = frame2im(frame);
% %     hold all
% %     
% % pause( 0.5 )
% % end
% % close;
% % 
% % filename = 'Percentage_change_mean_interp.gif'; % Specify the output file name
% % 
% %  for k = 1 : length(depth_interp)-1
% %     [A,map] = rgb2ind(im{k},256);
% %     if k == 1
% %         imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
% %     else
% %         imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
% %     end
% %  end
% %  
% %  close all
% % 
% % 
% % for k = 1 : length(depth_interp)
% % 
% %     
% % fig=figure;
% % m_proj('Robinson','lon',[-80 30],'lat',[-5 70])
% % 
% % fielddummy=squeeze(Percentage_change_hidrogen_ion(:,:,k));
% % 
% % for smoothnumber=1:7
% %     fielddummy=smooth2d(fielddummy);
% % end
% % 
% % m_pcolor(longitude,latitude,fielddummy); shading interp
% % 
% % %stippling
% % mask=Ptest_individualgrid(:,:,k)==100;%%this is the stippling mask
% % hold on
% % m_scatter(longitude(mask),latitude(mask),2,[.5 .5 .5])
% % 
% % set(gca, 'color', 'k');
% % m_coast('patch',[.7 .7 .7],'edgecolor','k');
% % m_grid('box','fancy','tickdir','out');
% % colormap(redblueleo);
% %  h=colorbar('northoutside');
% %  title(h,'Percentage Change of [ H^+ ] (%)','fontsize',12);
% %  caxis([-50  50]) 
% % 
% % h.TickLabels={'\leq-50';'-40';'-30';'-20';'-10';'0';'10';'20';'30';'40';'\geq50'};
% % patch([0.92 0.92 0.4 0.4],[0.10 -0.07 -0.07 0.10],'w','LineWidth',2); 
% % anotation=sprintf('Depth = %.2f m',depth_interp(k));
% % text(0.65, 0.02 ,anotation,'fontsize',10,'FontWeight','bold','color','k', 'vertical','middle','horizontal','center');
% %  
% % 
% % drawnow
% %     frame = getframe(fig);
% %     im{k} = frame2im(frame);
% %     hold all
% %     
% % pause( 0.5 )
% % end
% % close;
% % 
% % filename = 'Percentage_change_hydrogen_interp.gif'; % Specify the output file name
% % 
% %  for k = 1 : length(depth_interp)-1
% %     [A,map] = rgb2ind(im{k},256);
% %     if k == 1
% %         imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
% %     else
% %         imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
% %     end
% %  end
% % 
% % close all
% % 
% % 
% % 
% % for k = 1 : length(depth_interp)
% %     
% % fig=figure;
% % m_proj('Robinson','lon',[-80 30],'lat',[-5 70])
% % 
% % fielddummy=squeeze(Variance_Percentage_change(:,:,k));
% % 
% % for smoothnumber=1:7
% %     fielddummy=smooth2d(fielddummy);
% % end
% % 
% % m_pcolor(longitude,latitude,fielddummy); shading interp
% % 
% % %stippling
% % mask=Ptest_individualgrid(:,:,k)==100;%%this is the stippling mask
% % hold on
% % m_scatter(longitude(mask),latitude(mask),2,[.5 .5 .5])
% % 
% % 
% % set(gca, 'color', 'k');
% % m_coast('patch',[.7 .7 .7],'edgecolor','k');
% % m_grid('box','fancy','tickdir','out');
% % colormap(redblueleo);
% %  h=colorbar('northoutside');
% %  titlebar=('Percentage Change of pH \sigma^2 (%)');
% %  title(h,titlebar,'fontsize',12);
% %  caxis([-100  100]) 
% % h.TickLabels={'\leq-100';'-80';'-60';'-40';'-20';'0';'20';'40';'60';'80';'\geq100'};
% % patch([0.92 0.92 0.4 0.4],[0.10 -0.07 -0.07 0.10],'w','LineWidth',2); 
% % anotation=sprintf('Depth = %.2f m',depth_interp(k));
% % text(0.65, 0.02 ,anotation,'fontsize',10,'FontWeight','bold','color','k', 'vertical','middle','horizontal','center');
% %   
% % drawnow
% %     frame = getframe(fig);
% %     im{k} = frame2im(frame);
% %     hold all
% %     
% % pause( 0.5 )
% % end
% % close;
% % 
% % filename = 'Percentage_change_variance_interp.gif'; % Specify the output file name
% % 
% %  for k = 1 : length(depth_interp)-1
% %     [A,map] = rgb2ind(im{k},256);
% %     if k == 1
% %         imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
% %     else
% %         imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1);
% %     end
% %  end
% % 
% % close all


