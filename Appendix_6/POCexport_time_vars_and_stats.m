%% CODE

clear all
close all
warning('off','all') %deactivate warnings for speed

%% loading some data prior to variable loop
load('/Volumes/LaCie_Leonardo/NorESM/all_ramps/time_array.mat') %loading time array

%% Loading domain mask
domain_mask=load ('/Volumes/LaCie_Leonardo/NorESM/Mask_regions/Mask_corrected.mat');
domain_mask=domain_mask.domain_mask;

%% loading  depths
depth=load('/Volumes/LaCie_Leonardo/NorESM/all_ramps/Depth_Levels.mat'); %get depth array from one of the files in the folder
depth=depth.depths;

%% VARIABLE LOOP
varlist={'epc100'};
      
 varIDX=1;
    
    %% Importing Pre-Industrial 5-yr mean
    
    PI_mean=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/Pre_industrial/30yrMeans/PI_%s_first_5yr_mean.mat',varlist{varIDX}));
    PI_mean=PI_mean.PI_mean;
    
    PI_mean_bit1=PI_mean(200:end,:,:); %reshaping -- position 200 is where Pacific is
    PI_mean_bit2=PI_mean(1:199,:,:);
    PI_mean=[PI_mean_bit1 ; PI_mean_bit2];
    
    PI_mean=permute(PI_mean,[2,1,3]);%tranposing variables from (320X384 to 384X320) to match domain_mask
    
    %% Importing Pre-Industrial 30-yr Standard Deviation
    
    PI_sd=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/Pre_industrial/30yrMeans/PI_%s_30yr_sd.mat',varlist{varIDX}));
    PI_sd=PI_sd.PI_sd;
    
    PI_sd_bit1=PI_sd(200:end,:,:); %reshaping -- position 200 is where Pacific is
    PI_sd_bit2=PI_sd(1:199,:,:);
    PI_sd=[PI_sd_bit1 ; PI_sd_bit2];
    
    PI_sd=permute(PI_sd,[2,1,3]); %tranposing variables from (320X384 to 384X320) to match domain_mask
    
    clear *bit*
    
    %% declaring new variables
    
    tic
    
    dummyTOD=NaN(size(domain_mask,1),size(domain_mask,2)); %allocating memory space
    TOD=cell(1,3);
    TOD(:,:)={dummyTOD}; %allocating memory space. This is TOD based on a very tight envelope of +-1sd, normal envelope +- 2Sd and broader and +-3Sd envelope
     
    TpeakMAX=NaN(size(domain_mask,1),size(domain_mask,2));
    PeakVALUE=NaN(size(domain_mask,1),size(domain_mask,2));
    TimeLAG=NaN(size(domain_mask,1),size(domain_mask,2));
       
    dummyTrecovery=NaN(size(domain_mask,1),size(domain_mask,2));
    Trecovery=cell(1,3);
    Trecovery(:,:)={dummyTrecovery}; %allocating memory space. This is TOD based on a very tight envelope of +-1sd, normal envelope +- 2Sd and broader and +-3Sd envelope
      
    MitigationStart=140; %marks when Ramp-down began
    gap=9; %used in time-window for Tod and Trec
    
    k=100;
        
        fprintf('Variable %s | Depth Level=%1d \n',varlist{varIDX},k);
        ans2=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/all_ramps/filtered/new_series_%s_k_%d.mat',varlist{varIDX},k)); %%%Loading filtered series for that particular k level
        ans2=ans2.new_series;
        
        %% calculating Tpeak and value associated, which is the maximum absolute difference
        
        sprintf('calculating Time Lowest')
        
        for l=1:size(domain_mask,1) %setting latitudinal limits
            for c=1:size(domain_mask,2) %setting longitudinal limits
                
                if domain_mask(l,c)==2 && ~isnan(PI_mean(l,c)) %if it's Atlantic and PI_pH exists
                    
                    seriesPH_interpolated=ans2{l,c};
                    %[valmin col]=min(seriesPH_interpolated);%finding the minimum in the series
                    [valmax , colmax]=max(abs(seriesPH_interpolated-PI_mean(l,c)));
                    TpeakMAX(l,c)=colmax;
                    PeakVALUE(l,c)=seriesPH_interpolated(colmax);
 
                else
                    TpeakMAX(l,c)=NaN;
                    PeakVALUE(l,c)=NaN;
                    TimeLAG(l,c)=NaN;
                end
            end
        end
        
     
        %% Time of Departure (ToD)
        
        sprintf('Calculating Time of Departure')
        
        for l=1:size(domain_mask,1) %setting longitudinal limits
            for c=1:size(domain_mask,2)
                if domain_mask(l,c)==2 && ~isnan(PI_mean(l,c)) %if it's Atlantic and PI_values exist
                    
                    %loading pre-ind mean value for depth k
                    pre_ind_gridcell=PI_mean(l,c);
                    
                    %loading pre-ind mean SD for depth k
                    pre_ind_sd=PI_sd(l,c);
                    
                    
                    PI_ub =   PI_mean(l,c) + [1*PI_sd(l,c) 2*PI_sd(l,c) 3*PI_sd(l,c)]; %upper boundaries
                    PI_lwb = PI_mean(l,c) - [1*PI_sd(l,c) 2*PI_sd(l,c) 3*PI_sd(l,c)]; %lower boundaries
                    
                    
                    %                     %lower and upper boundaries for pre-industrial mean for depth k
                    %                     PI_lwb=PI_mean(l,c)-2*pre_ind_sd; %lower boundary times 2 SD
                    %                     PI_ub=PI_mean(l,c)+2*pre_ind_sd; %upper boundary times 2 SD
                    
                    %importing time series
                    series_complete=ans2{l,c};
                    
                    for indx=1:length(PI_ub)
                        
                        %if series leaves any of the boundaries, mark time position when this happens
                        [colTODup]=(find(series_complete>PI_ub(indx))); %from above
                        [colTODdown]=(find(series_complete<PI_lwb(indx))); %from below
                        
                        noreturnup=[]; noreturndown=[];
                        
                        if ~isempty(colTODup) || ~isempty(colTODdown) %if indeed there is a point where it leaves the envelope
                            
                            if ~isempty(colTODup) %leaving envelope from the top
                                dummyTOP=diff(colTODup)==1;
                                Mtop = movsum(dummyTOP,[0 gap]);
                                noreturnup=find(Mtop==10,1,'first');
                            end
                            
                            if ~isempty(colTODdown) %leaving envelope from below
                                dummydown=diff(colTODdown)==1;
                                Mdown = movsum(dummydown,[0 gap]);
                                noreturndown=find(Mdown==10,1,'first');
                            end
                            
                            if ~isempty(noreturnup) && ~isempty(noreturndown) && ~isempty(colTODup) && ~isempty(colTODdown)  %if leaves both from below and from the top
                                dummyTOD(l,c)=min([colTODdown(noreturndown); colTODup(noreturnup)]); %take the longest time
                                
                            elseif ~isempty(noreturnup) && isempty(noreturndown) %if leaves only from the top
                                dummyTOD(l,c)=colTODup(noreturnup);
                                
                            elseif isempty(noreturnup) && ~isempty(noreturndown) % if it leaves only from below
                                dummyTOD(l,c)=colTODdown(noreturndown);
                            end
                        end
                        TOD{1,indx}(l,c)=dummyTOD(l,c);
                        clear dummyTOP Mtop noreturnup dummydown Mdown noreturndown
                    end
                else
                    TOD{1,1}(l,c)= NaN;
                    TOD{1,2}(l,c)= NaN;
                    TOD{1,3}(l,c)= NaN;
                end
                clear dummyTOP Mtop noreturnup dummydown Mdown noreturndown
            end
        end
        
  
        %% Time lag
        sprintf('Calculating Time lag')
        
        for l=1:size(domain_mask,1) %setting latitudinal limits
            for c=1:size(domain_mask,2) %setting longitudinal limits
                
                if domain_mask(l,c)==2 && ~isnan(PI_mean(l,c))%if it's Atlantic and PI_pH exists
                    
                    TimeLAG(l,c)=TpeakMAX(l,c)-MitigationStart; %subtract year of last minimum from the year the mitigation phase started
                    
                else
                    TimeLAG(l,c)=NaN;
                end
            end
        end
        
        
        
        %% Trecovery
        sprintf('Calculating Trecovery')
        
        %flags to be used
        badflagslow=-9999; %for slow recovery
        badflagsnegslope=-7777; %for negative slope (declining trend)
        badflagtooearly=-4444; %for early recovery
        bottomflag=-1111; %bottom
        
        reg_window=100; %number of year-points to be considered when doing Trec linear regression
        tic
        
        for l=1:size(domain_mask,1) %setting latitudinal limits
            for c=1:size(domain_mask,2) %setting longitudinal limits
                
                if domain_mask(l,c)==2 && ~isnan(PI_mean(l,c)) && ~isnan(PI_sd(l,c)) %if it's Atlantic and PI_values exist
                    
                    %building the x of the time series for our regression
                    time_recovery_regress=time_lowest_regress(end-reg_window:end);
                    
                    %loading pre-ind mean value for depth k
                    pre_ind_gridcell=PI_mean(l,c);
                    
                    %loading pre-ind mean SD for depth k
                    pre_ind_sd=PI_sd(l,c);
                    
                    %                     %lower and upper boundaries for pre-industrial mean for depth k
                    %                     PI_lwb=PI_mean(l,c)-2*pre_ind_sd; %lower boundary (2 times SD)
                    %                     PI_ub=PI_mean(l,c)+2*pre_ind_sd; %upper boundary (2 times SD)
                    
                    PI_ub =   PI_mean(l,c) + [1*PI_sd(l,c) 2*PI_sd(l,c) 3*PI_sd(l,c)]; %upper boundaries
                    PI_lwb = PI_mean(l,c) - [1*PI_sd(l,c) 2*PI_sd(l,c) 3*PI_sd(l,c)]; %lower boundaries
                    
                    %importing time series
                    series_regress=ans2{l,c}; %load timeseries
                    
                    for indx=1:length(PI_ub)
                        
                        
                        %% if series dives into the envelope from above upper boundary
                        if PeakVALUE(l,c)>PI_ub(indx) %series_regress(140)>PI_ub
                            
                            %finding whether there is a first point after the start
                            %of CO2atm ramp down to enter time series from above
                            point_envelope_up=find(series_regress(TpeakMAX(l,c):end)<PI_ub(indx));
                            
                            if ~isempty(point_envelope_up) %in case there are points that enter the envelope take the index difference to find the first 10 consecutive points
                                dummy2up=diff(point_envelope_up)==1;
                                M2up = movsum(dummy2up,[0 gap]);
                                noreturn2up=find(M2up==10,1,'first');
                                
                                if ~isempty(noreturn2up) %if there are 10 consecutive points
                                    trecovery_dummy_up=point_envelope_up(noreturn2up)+MitigationStart-1;
                                    Trecovery{1,indx}(l,c)=trecovery_dummy_up;
                                    
                                else  %if there are no 10 consecutive points then do regression
                                    %calculate polyfit based on extension ramp points
                                    
                                    curve=fit(time_recovery_regress,series_regress(end-reg_window:end)','poly1'); %polynomial fit
                                    polinomyal_up=[curve.p1 curve.p2-PI_ub(indx)]; %coefficients
                                    trecovery_dummy_up=round(roots(polinomyal_up)); %solving and rounding polynomial roots
                                    Trecovery{1,indx}(l,c)=trecovery_dummy_up; %Trec
                                    
                                    %EXPLANATION FOR A 3rd ORDER POLYNOM
                                    %curve=fit(time_recovery_regress,series_regress(end-reg_window:end)','poly3');
                                    
                                    %equation for any Y is given by the coeficients of the curve
                                    %Y=(curve.p1)*x^3+ curve.p2)*x^2+ (curve.p3)*x+ curve.p4 ;
                                    %if I want to solve for Y=PI_lwb I just need to addapt the equation
                                    %PI_lwb = (curve.p1*x)^3+ (curve.p2*x)^2+ curve.p3*x+ curve.p4;
                                    %0 = (curve.p1*x)^3+ (curve.p2*x)^2+ curve.p3*x+ curve.p4 -PI_lwb
                                    %then the 'roots' would correspond to an Y=PI_lwb
                                    % pol=[curve.p1 curve.p2 curve.p3 curve.p4];
                                    % polyval(polinomial,140:700);
                                    %
                                    % polinomyal_up=[curve.p1 curve.p2 curve.p3 curve.p4-PI_ub]; %specifying the coefficients 'note the last coeficient incorporates PI_lwb
                                    % dummyroots_up=roots(polinomyal_up); %finding the roots (in this case 3) and picking the earliest one which its imaginary part is zero
                                    %
                                    % [indxroot_up]=find(imag(dummyroots_up)==0); %get only the roots where imaginary part is zero.
                                    % dummyroots_up=round(dummyroots_up(indxroot_up));
                                    %
                                    % trecovery_dummy_up=min(dummyroots_up(dummyroots_up>200 & dummyroots_up<1000)); %then take only the real part (because it's still real+000i), round it and then take the earliest time between 280 and 1000
                                    % Trecovery(l,c)=trecovery_dummy_up; %this finds the intercept after the start of the extension and it has to bee less than 1000yrs otherwise no recovery
                                    
                                end
                                
                            elseif isempty(point_envelope_up) || isempty(noreturn2up) %if there are no points to enter envelope do a regression using reg window of the series anyway
                                curve=fit(time_recovery_regress,series_regress(end-reg_window:end)','poly1');
                                polinomyal_up=[curve.p1 curve.p2-PI_ub(indx)];
                                trecovery_dummy_up=round(roots(polinomyal_up));
                                Trecovery{1,indx}(l,c)=trecovery_dummy_up;
                            end
                        end
                        
                        
                        %% if series dives into the envelope from lower boundary
                        if PeakVALUE(l,c)<PI_lwb(indx) %PI_ub series_regress(140)<PI_lwb
                            %finding whether there is a first point after the start
                            %of CO2atm ramp down to enter time series from above
                            point_envelope_down=find(series_regress(TpeakMAX(l,c):end)>PI_lwb(indx));
                            noreturn2down=[]; trecovery_dummy_down=[];
                            
                            if ~isempty(point_envelope_down) %in case there are points that enter the envelope take the index difference to find the first 10 consecutive points
                                dummy2down=diff(point_envelope_down)==1;
                                M2down = movsum(dummy2down,[0 gap]);
                                noreturn2down=find(M2down==10,1,'first');
                                
                                if ~isempty(noreturn2down) %if there are 10 consecutive points
                                    trecovery_dummy_down=point_envelope_down(noreturn2down)+MitigationStart-1;
                                    Trecovery{1,indx}(l,c)=trecovery_dummy_down;
                                    
                                else%if there are no 10 consecutive points then do regression
                                    %calculate polyfit based on extension ramp points
                                    curve=fit(time_recovery_regress,series_regress(end-reg_window:end)','poly1');
                                    polinomyal_down=[curve.p1 curve.p2-PI_lwb(indx)];
                                    trecovery_dummy_down=round(roots(polinomyal_down));
                                    Trecovery{1,indx}(l,c)=trecovery_dummy_down;
                                end
                                
                            elseif isempty(point_envelope_down) || isempty(noreturn2down) %if there are no points to enter envelope do a regression using reg window of the series anyway
                                
                                curve=fit(time_recovery_regress,series_regress(end-reg_window:end)','poly1');
                                polinomyal_down=[curve.p1 curve.p2-PI_lwb(indx)];
                                trecovery_dummy_down=round(roots(polinomyal_down));
                                Trecovery{1,indx}(l,c)=trecovery_dummy_down;
                                
                            end
                        end
                        
                        %                     if series_regress(140)>PI_lwb && series_regress(140)<PI_ub %if already inside the envelope then early recovery
                        %                         Trecovery(l,c)=-4444;
                        %                     end
                        
                    end
                elseif domain_mask(l,c)==2 && isnan(PI_mean(l,c)) %if is atlantic but PI mean does not exist
                    Trecovery{1,1}(l,c)=bottomflag;
                    Trecovery{1,2}(l,c)=bottomflag;
                    Trecovery{1,3}(l,c)=bottomflag;
                end
            end
        end
  
        
%         figure
%         pcolor(Trecovery);shading flat
%         caxis([120 520])
%         colorbar
    
    sprintf('Saving all results \n')
    save(sprintf('/Volumes/LaCie_Leonardo/NorESM/all_ramps/filtered/Results/%s_results_updated_v1D_poly1_NEW_NEW.mat',varlist{varIDX}), 'Trecovery', 'TOD', 'TimeLAG', 'PeakVALUE', 'TpeakMAX')
    
    
    
    
    %% STATISTICAL ANALYSES
    sprintf('calculating stats \n')
    
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
            if sum(isnan(PI_detrend{l,c}))==length(PI_detrend{l,c}) %whenever series have 30 invalid datapoints, assing a single NAN to that cell element
                PI_detrend{l,c}=NaN;
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

PI_sd=permute(PI_sd,[2,1]);
PI_detrend=permute(PI_detrend,[2,1]);

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
Ptest_individualgrid=zeros(size(latitude,1),size(latitude,2)) ; %flag for doimain being continent or bottom
Percentage_change_mean=zeros(size(latitude,1),size(latitude,2)) ;
Variance_Percentage_change=zeros(size(latitude,1),size(latitude,2)) ;

Elem_sgnif=zeros(1,1);
Elem_notsgnif=zeros(1,1) ;
Elem_notnormal=zeros(1,1) ;

%% Statistics on Vertically Interpolated fields (Now the PI and PR fields are spaced evenly in the vertical) 

      %loading filtered series and
%     %extracting last 30yrs for each grid-cell at each depth
    
PR_series=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/all_ramps/filtered/new_series_%s_k_%d.mat',varlist{varIDX},k)); 
PR_series=PR_series.new_series; 

sprintf('Var = %s depth level %d\n',varlist{varIDX},k)

    for l=1:size(PR_series,1)
        l
        for c=1:size(PR_series,2)
            
            if (sum(~isnan(PR_series{l,c}(:)))==481) %if the cell actually contains a time series of 90 points (3 windows of 30)
               
                %PRE INDUSTRIAL
                PI_last30=PI_detrend{l,c};            %selecting last 30years of Pre_industrial
                PI_last30_mean=mean(PI_detrend{l,c}); %Pre-Industrial 30yr Mean
                
                %End of mitigation
                EndMitigation30=PR_series{l,c}(251:280);
                dummy3=detrend(EndMitigation30);                   %removing trend
                EndMitigation30_mean=mean(EndMitigation30);        %Mitigation end  30yr mean
                EndMitigation30=dummy3 + EndMitigation30_mean;
                
                %Extension-Middle
                Ext_middle30=PR_series{l,c}(351:380);    %selecting last 30 years of extension time
                dummy2=detrend(Ext_middle30);             %removing trend
                Ext_middle30_mean=mean(Ext_middle30);        %Extension middle  30yr mean
                Ext_middle30=dummy2 + Ext_middle30_mean;     %adding increments without thend to the mean value to get detrended series
                
                %Extension-END
                Ext_last30=PR_series{l,c}(end-29:end);    %selecting last 30 years of extension time
                dummy1=detrend(Ext_last30);             %removing trend
                Ext_last30_mean=mean(Ext_last30);        %Extension time 30yr mean
                Ext_last30=dummy1 + Ext_last30_mean;     %adding increments without thend to the mean value to get detrended series
                
                % MEAN PERCENTAGE CHANGES PAIRWISE
                addpath /Volumes/LaCie_Leonardo/NorESM
                Percentage_change_mean_pair1(l,c)=percentagechangeNorESM(PI_last30_mean,EndMitigation30_mean);
                Percentage_change_mean_pair2(l,c)=percentagechangeNorESM(PI_last30_mean,Ext_middle30_mean);
                Percentage_change_mean_pair3(l,c)=percentagechangeNorESM(PI_last30_mean,Ext_last30_mean);
                
              %DELTAS
              delta_pair1(l,c)=EndMitigation30_mean-PI_last30_mean;
              delta_pair2(l,c)=Ext_middle30_mean-PI_last30_mean;
              delta_pair3(l,c)=Ext_last30_mean-PI_last30_mean;
                
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
                Pflags_pair1(l,c) = ttestNorESM(PI_last30,EndMitigation30); %%PI versus Extension END
                Pflags_pair2(l,c) = ttestNorESM(PI_last30,Ext_middle30); %%PI versus Extension END
                Pflags_pair3(l,c) = ttestNorESM(PI_last30,Ext_last30); %%PI versus Extension END

                %                             figure
                %                             boxplot([PI_last30 EndMitigation30 Ext_middle30 Ext_last30],groupnames)
                %                             title (sprintf('line %d col %d depth %d',l,c,k))
                %                             xlabel('Model Phase')
                %                             ylabel('pH')
                
                else %if does not contain area of atlantic ocean or falls onto continents or bottom
                
                Percentage_change_mean_pair1(l,c)=NaN;   Pflags_pair1(l,c)=NaN;
                Percentage_change_mean_pair2(l,c)=NaN;   Pflags_pair2(l,c)=NaN;
                Percentage_change_mean_pair3(l,c)=NaN;   Pflags_pair3(l,c)=NaN;
                
              delta_pair1(l,c)=NaN;
              delta_pair2(l,c)=NaN;
              delta_pair3(l,c)=NaN;
           
 
            end  
        end
    end
    
    
    figure
    pcolor(Percentage_change_mean_pair3);shading flat
    colorbar
    caxis([-20 20])
    
    figure
    pcolor(Pflags_pair3);shading flat
    colormap(parula(3))
    colorbar
    
    
    
       % Points with significant difference
    Elem_sgnif_pair1(1,1)= numel(find(Pflags_pair1(:,:)==100))/numel(find(~isnan(Pflags_pair1(:,:))));
    Elem_sgnif_pair2(1,1)= numel(find(Pflags_pair2(:,:)==100))/numel(find(~isnan(Pflags_pair2(:,:))));
    Elem_sgnif_pair3(1,1)= numel(find(Pflags_pair3(:,:)==100))/numel(find(~isnan(Pflags_pair3(:,:))));
    
   
    % Points with no significant difference
    Elem_notsgnif_pair1(1,1)=numel(find(Pflags_pair1(:,:)==-100))/numel(find(~isnan(Pflags_pair1(:,:))));
    Elem_notsgnif_pair2(1,1)=numel(find(Pflags_pair2(:,:)==-100))/numel(find(~isnan(Pflags_pair2(:,:))));
    Elem_notsgnif_pair3(1,1)=numel(find(Pflags_pair3(:,:)==-100))/numel(find(~isnan(Pflags_pair3(:,:))));

    % Points discarded because normality test showed series were not normally distributed 
    Elem_notnormal_pair1(1,1)=numel(find(Pflags_pair1(:,:)==-200))/numel(find(~isnan(Pflags_pair1(:,:))));
    Elem_notnormal_pair2(1,1)=numel(find(Pflags_pair2(:,:)==-200))/numel(find(~isnan(Pflags_pair2(:,:))));
    Elem_notnormal_pair3(1,1)=numel(find(Pflags_pair3(:,:)==-200))/numel(find(~isnan(Pflags_pair3(:,:))));
    
    
TABLE_STATS_pair1=table(Elem_sgnif_pair1,Elem_notsgnif_pair1,Elem_notnormal_pair1,100); %table with depth strata and their respective fractions
TABLE_STATS_pair2=table(Elem_sgnif_pair2,Elem_notsgnif_pair2,Elem_notnormal_pair2,100); %table with depth strata and their respective fractions
TABLE_STATS_pair3=table(Elem_sgnif_pair3,Elem_notsgnif_pair3,Elem_notnormal_pair3,100); %table with depth strata and their respective fractions


sprintf('Saving all results')
save (sprintf('/Volumes/LaCie_Leonardo/NorESM/StatsResults/%s_stats_updated_v1D.mat',varlist{varIDX}),...
    'delta_pair1','delta_pair2','delta_pair3',...
    'Percentage_change_mean_pair1', 'Percentage_change_mean_pair2','Percentage_change_mean_pair3',...
    'Pflags_pair1','Pflags_pair2','Pflags_pair3',...
    'TABLE_STATS_pair1', 'TABLE_STATS_pair2', 'TABLE_STATS_pair3')
