%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TIME VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculates for all biogeochem filtered variables (pH O2 OMEGAc T and Detoc):
% ToD
% TPeak
% PeakVALUE
% TimeLAG
% Trecovery
% Saves everyhting into a single mat file 'VARIABLE_results_updated_v1D_poly1.mat'

% Author: Leonardo Bertini
% Modified: 22nd April 2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
varlist={'ph';'o2';'omegac';'templvl';'detoc';'AOU'};
tic      
for varIDX=7 %length(varlist)
    
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
    
    
    
    dummyTOD=NaN(size(domain_mask,1),size(domain_mask,2),size(depth,1)); %allocating memory space
    TOD=cell(1,3);
    TOD(:,:)={dummyTOD}; %allocating memory space. This is TOD based on a very tight envelope of +-1sd, normal envelope +- 2Sd and broader and +-3Sd envelope
    
    TpeakMAX=NaN(size(domain_mask,1),size(domain_mask,2),size(depth,1));
    PeakVALUE=NaN(size(domain_mask,1),size(domain_mask,2),size(depth,1));
    TimeLAG=NaN(size(domain_mask,1),size(domain_mask,2),size(depth,1));
    
    dummyTrecovery=NaN(size(domain_mask,1),size(domain_mask,2),size(depth,1));
    Trecovery=cell(1,3);
    Trecovery(:,:)={dummyTrecovery}; %allocating memory space. This is TOD based on a very tight envelope of +-1sd, normal envelope +- 2Sd and broader and +-3Sd envelope
    
  
    MitigationStart=140; %marks when Ramp-down began
    gap=9; %used in time-window for Tod and Trec
    
    for  k=1:length(depth)
        
        fprintf('Variable %s | Depth Level=%1d \n',varlist{varIDX},k);
        ans2=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/all_ramps/filtered/new_series_%s_k_%d.mat',varlist{varIDX},k)); %%%Loading filtered series for that particular k level
        ans2=ans2.new_series;
        
        %% calculating Tpeak and value associated, which is the maximum absolute difference
        
        sprintf('calculating Time Lowest')
        
        for l=1:size(domain_mask,1) %setting latitudinal limits
            for c=1:size(domain_mask,2) %setting longitudinal limits
                
                if domain_mask(l,c)==2 && ~isnan(PI_mean(l,c,k)) %if it's Atlantic and PI_pH exists
                    
                    seriesPH_interpolated=ans2{l,c};
                    %[valmin col]=min(seriesPH_interpolated);%finding the minimum in the series
                    [valmax , colmax]=max(abs(seriesPH_interpolated-PI_mean(l,c,k)));
                    TpeakMAX(l,c,k)=colmax;
                    PeakVALUE(l,c,k)=seriesPH_interpolated(colmax);
 
                else
                    TpeakMAX(l,c,k)=NaN;
                    PeakVALUE(l,c,k)=NaN;
                    TimeLAG(l,c,k)=NaN;
                end
            end
        end
     
        %% Time of Departure (ToD) - with 3 different envelopes (1sd 2sd and 3sd)
        
        sprintf('Calculating Time of Departure \n with 3 differente envelope widths 1, 2 and 3 sd')

        for l=1:size(domain_mask,1) %setting longitudinal limits
            for c=1:size(domain_mask,2)
                if domain_mask(l,c)==2 && ~isnan(PI_mean(l,c,k)) %if it's Atlantic and PI_values exist
                    
                    %loading pre-ind mean value for depth k
                    pre_ind_gridcell=PI_mean(l,c,k);
                    
                    %loading pre-ind mean SD for depth k
                    pre_ind_sd=PI_sd(l,c,k);
                    
                    
                    PI_ub =   PI_mean(l,c,k) + [1*PI_sd(l,c,k) 2*PI_sd(l,c,k) 3*PI_sd(l,c,k)]; %upper boundaries
                    PI_lwb = PI_mean(l,c,k) - [1*PI_sd(l,c,k) 2*PI_sd(l,c,k) 3*PI_sd(l,c,k)]; %lower boundaries
                    
                    
                    %lower and upper boundaries for pre-industrial mean for depth k
                    %PI_lwb=PI_mean(l,c,k)-2*pre_ind_sd; %lower boundary times 2 SD
                    %PI_ub=PI_mean(l,c,k)+2*pre_ind_sd; %upper boundary times 2 SD
                    
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
                                dummyTOD(l,c,k)=min([colTODdown(noreturndown); colTODup(noreturnup)]); %take the longest time
                                
                            elseif ~isempty(noreturnup) && isempty(noreturndown) %if leaves only from the top
                                dummyTOD(l,c,k)=colTODup(noreturnup);
                                
                            elseif isempty(noreturnup) && ~isempty(noreturndown) % if it leaves only from below
                                dummyTOD(l,c,k)=colTODdown(noreturndown);
                            end
                        end  
                     TOD{1,indx}(l,c,k)=dummyTOD(l,c,k);
                     clear dummyTOP Mtop noreturnup dummydown Mdown noreturndown
                    end
                    
                else
                    TOD{1,1}(l,c,k)= NaN;
                    TOD{1,2}(l,c,k)= NaN;
                    TOD{1,3}(l,c,k)= NaN;

                    clear dummyTOP Mtop noreturnup dummydown Mdown noreturndown
                end
                
            end
        end
        
        
%         figure
%         subplot(1,3,1)
%         pcolor(TOD{1,1}(:,:,1));shading flat
%         colorbar
%         caxis([0 20])
%         
%         subplot(1,3,2)
%         pcolor(TOD{1,2}(:,:,1));shading flat
%         colorbar
%         caxis([0 20])
%         
%         subplot(1,3,3)
%         pcolor(TOD{1,3}(:,:,1));shading flat
%         colorbar
%         caxis([0 20])
        
%          %% Time of Departure (ToD) - 2sd
%         
%         sprintf('Calculating Time of Departure - 2sd')
% 
%         for l=1:size(domain_mask,1) %setting longitudinal limits
%             for c=1:size(domain_mask,2)
%                 if domain_mask(l,c)==2 && ~isnan(PI_mean(l,c,k)) %if it's Atlantic and PI_values exist
%                     
%                     %loading pre-ind mean value for depth k
%                     pre_ind_gridcell=PI_mean(l,c,k);
%                     
%                     %loading pre-ind mean SD for depth k
%                     pre_ind_sd=PI_sd(l,c,k);
%                     
%                     %lower and upper boundaries for pre-industrial mean for depth k
%                     PI_lwb=PI_mean(l,c,k)-2*pre_ind_sd; %lower boundary times 2 SD
%                     PI_ub=PI_mean(l,c,k)+2*pre_ind_sd; %upper boundary times 2 SD
%                     
%                     %importing time series
%                     series_complete=ans2{l,c};
%                     
%                     %if series leaves any of the boundaries, mark time position when this happens
%                     [colTODup]=(find(series_complete>PI_ub)); %from above
%                     [colTODdown]=(find(series_complete<PI_lwb)); %from below
%                     
%                     noreturnup=[]; noreturndown=[];
% 
%                     if ~isempty(colTODup) || ~isempty(colTODdown) %if indeed there is a point where it leaves the envelope
%                         
%                         if ~isempty(colTODup) %leaving envelope from the top
%                             dummyTOP=diff(colTODup)==1;
%                             Mtop = movsum(dummyTOP,[0 gap]);
%                             noreturnup=find(Mtop==10,1,'first');
%                         end
%                         
%                         if ~isempty(colTODdown) %leaving envelope from below
%                             dummydown=diff(colTODdown)==1;
%                             Mdown = movsum(dummydown,[0 gap]);
%                             noreturndown=find(Mdown==10,1,'first');
%                         end
%                         
%                         if ~isempty(noreturnup) && ~isempty(noreturndown) && ~isempty(colTODup) && ~isempty(colTODdown)  %if leaves both from below and from the top
%                             TOD(l,c,k)=min([colTODdown(noreturndown); colTODup(noreturnup)]); %take the longest time
%                             
%                         elseif ~isempty(noreturnup) && isempty(noreturndown) %if leaves only from the top
%                             TOD(l,c,k)=colTODup(noreturnup);
%                             
%                         elseif isempty(noreturnup) && ~isempty(noreturndown) % if it leaves only from below
%                             TOD(l,c,k)=colTODdown(noreturndown);
%                         end
%                     end
%                     
%                 else
%                     TOD(l,c,k)= NaN;
%                 end
%                 clear dummyTOP Mtop noreturnup dummydown Mdown noreturndown
%             end
%         end
%         
        
        
        
        
        
        
  
        %% Time lag 
        sprintf('Calculating Time lag')
        
        for l=1:size(domain_mask,1) %setting latitudinal limits
            for c=1:size(domain_mask,2) %setting longitudinal limits
                
                if domain_mask(l,c)==2 && ~isnan(PI_mean(l,c,k))%if it's Atlantic and PI_pH exists
                    
                    TimeLAG(l,c,k)=TpeakMAX(l,c,k)-MitigationStart; %subtract year of last minimum from the year the mitigation phase started
                    
                else
                    TimeLAG(l,c,k)=NaN;
                end
            end
        end
        
        
        
        %% Trecovery with 3 different envelopes (1sd 2sd and 3sd)
        sprintf('Calculating Trecovery')
        
        %flags to be used
        badflagslow=-9999; %for slow recovery
        badflagsnegslope=-7777; %for negative slope (declining trend)
        badflagtooearly=-4444; %for early recovery
        bottomflag=-1111; %bottom
        
        reg_window=100; %number of year-points to be considered when doing Trec linear regression
        
        for l=1:size(domain_mask,1) %setting latitudinal limits
            for c=1:size(domain_mask,2) %setting longitudinal limits
                
                if domain_mask(l,c)==2 && ~isnan(PI_mean(l,c,k)) && ~isnan(PI_sd(l,c,k)) %if it's Atlantic and PI_values exist
                    
                    %building the x of the time series for our regression
                    time_recovery_regress=time_lowest_regress(end-reg_window:end);
                    
                    %loading pre-ind mean value for depth k
                    pre_ind_gridcell=PI_mean(l,c,k);
                    
                    %loading pre-ind mean SD for depth k
                    pre_ind_sd=PI_sd(l,c,k);
                    
                    
                    PI_ub =   PI_mean(l,c,k) + [1*PI_sd(l,c,k) 2*PI_sd(l,c,k) 3*PI_sd(l,c,k)]; %upper boundaries
                    PI_lwb = PI_mean(l,c,k) - [1*PI_sd(l,c,k) 2*PI_sd(l,c,k) 3*PI_sd(l,c,k)]; %lower boundaries
                    
                    %                     %lower and upper boundaries for pre-industrial mean for depth k
                    %                     PI_lwb=PI_mean(l,c,k)-2*pre_ind_sd; %lower boundary (2 times SD)
                    %                     PI_ub=PI_mean(l,c,k)+2*pre_ind_sd; %upper boundary (2 times SD)
                    %
                    %importing time series
                    series_regress=ans2{l,c}; %load timeseries
                    
                    
                    for indx=1:length(PI_ub)
                        
                        %% if series dives into the envelope from above upper boundary
                        if PeakVALUE(l,c,k)>PI_ub(indx) %series_regress(140)>PI_ub
                            
                            %finding whether there is a first point after the start
                            %of CO2atm ramp down to enter time series from above
                            point_envelope_up=find(series_regress(TpeakMAX(l,c,k):end)<PI_ub(indx));
                            noreturn2up=[]; trecovery_dummy_up=[];
                            
                            if ~isempty(point_envelope_up) %in case there are points that enter the envelope take the index difference to find the first 10 consecutive points
                                dummy2up=diff(point_envelope_up)==1;
                                M2up = movsum(dummy2up,[0 gap]);
                                noreturn2up=find(M2up==10,1,'first');
                                
                                if ~isempty(noreturn2up) %if there are 10 consecutive points
                                    trecovery_dummy_up=point_envelope_up(noreturn2up)+MitigationStart-1;
                                    Trecovery{1,indx}(l,c,k)=trecovery_dummy_up;
                                    
                                else  %if there are no 10 consecutive points then do regression
                                    %calculate polyfit based on extension ramp points
                                    
                                    curve=fit(time_recovery_regress,series_regress(end-reg_window:end)','poly1'); %polynomial fit
                                    polinomyal_up=[curve.p1 curve.p2-PI_ub(indx)]; %coefficients
                                    trecovery_dummy_up=round(roots(polinomyal_up)); %solving and rounding polynomial roots
                                    Trecovery{1,indx}(l,c,k)=trecovery_dummy_up; %Trec
                                    
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
                                    % Trecovery(l,c,k)=trecovery_dummy_up; %this finds the intercept after the start of the extension and it has to bee less than 1000yrs otherwise no recovery
                                    
                                end
                                
                            elseif isempty(point_envelope_up) || isempty(noreturn2up) %if there are no points to enter envelope do a regression using reg window of the series anyway
                                curve=fit(time_recovery_regress,series_regress(end-reg_window:end)','poly1');
                                polinomyal_up=[curve.p1 curve.p2-PI_ub(indx)];
                                trecovery_dummy_up=round(roots(polinomyal_up));
                                Trecovery{1,indx}(l,c,k)=trecovery_dummy_up;
                            end
                        end
                        
                        
                        %% if series dives into the envelope from lower boundary
                        if PeakVALUE(l,c,k)<PI_lwb(indx) %PI_ub series_regress(140)<PI_lwb
                            %finding whether there is a first point after the start
                            %of CO2atm ramp down to enter time series from above
                            point_envelope_down=find(series_regress(TpeakMAX(l,c,k):end)>PI_lwb(indx));
                            noreturn2down=[]; trecovery_dummy_down=[];
                            
                            if ~isempty(point_envelope_down) %in case there are points that enter the envelope take the index difference to find the first 10 consecutive points
                                dummy2down=diff(point_envelope_down)==1;
                                M2down = movsum(dummy2down,[0 gap]);
                                noreturn2down=find(M2down==10,1,'first');
                                
                                if ~isempty(noreturn2down) %if there are 10 consecutive points
                                    trecovery_dummy_down=point_envelope_down(noreturn2down)+MitigationStart-1;
                                    Trecovery{1,indx}(l,c,k)=trecovery_dummy_down;
                                    
                                else%if there are no 10 consecutive points then do regression
                                    %calculate polyfit based on extension ramp points
                                    curve=fit(time_recovery_regress,series_regress(end-reg_window:end)','poly1');
                                    polinomyal_down=[curve.p1 curve.p2-PI_lwb(indx)];
                                    trecovery_dummy_down=round(roots(polinomyal_down));
                                    Trecovery{1,indx}(l,c,k)=trecovery_dummy_down;
                                end
                                
                            elseif isempty(point_envelope_down) || isempty(noreturn2down) %if there are no points to enter envelope do a regression using reg window of the series anyway
                                
                                curve=fit(time_recovery_regress,series_regress(end-reg_window:end)','poly1');
                                polinomyal_down=[curve.p1 curve.p2-PI_lwb(indx)];
                                trecovery_dummy_down=round(roots(polinomyal_down));
                                Trecovery{1,indx}(l,c,k)=trecovery_dummy_down;
                                
                            end
                            
                            
                            
                        end
                        
                        %                     if series_regress(140)>PI_lwb && series_regress(140)<PI_ub %if already inside the envelope then early recovery
                        %                         Trecovery(l,c,k)=-4444;
                        %                     end
                    end
                    
                elseif domain_mask(l,c)==2 && isnan(PI_mean(l,c,k)) %if is atlantic but PI mean does not exist
                    Trecovery{1,1}(l,c,k)=bottomflag;
                    Trecovery{1,2}(l,c,k)=bottomflag;
                    Trecovery{1,3}(l,c,k)=bottomflag;
                end
                
            end
        end
    
        
        
        
        
        
        
        
%         %% Trecovery 2sd
%         sprintf('Calculating Trecovery')
%         
%         %flags to be used
%         badflagslow=-9999; %for slow recovery
%         badflagsnegslope=-7777; %for negative slope (declining trend)
%         badflagtooearly=-4444; %for early recovery
%         bottomflag=-1111; %bottom
%         
%         reg_window=100; %number of year-points to be considered when doing Trec linear regression
%         tic
%         
%         for l=1:size(domain_mask,1) %setting latitudinal limits
%             for c=1:size(domain_mask,2) %setting longitudinal limits
%                 
%                 if domain_mask(l,c)==2 && ~isnan(PI_mean(l,c,k)) && ~isnan(PI_sd(l,c,k)) %if it's Atlantic and PI_values exist
%                     
%                     %building the x of the time series for our regression
%                     time_recovery_regress=time_lowest_regress(end-reg_window:end);
%                     
%                     %loading pre-ind mean value for depth k
%                     pre_ind_gridcell=PI_mean(l,c,k);
%                     
%                     %loading pre-ind mean SD for depth k
%                     pre_ind_sd=PI_sd(l,c,k);
%                     
%                     %lower and upper boundaries for pre-industrial mean for depth k
%                     PI_lwb=PI_mean(l,c,k)-2*pre_ind_sd; %lower boundary (2 times SD)
%                     PI_ub=PI_mean(l,c,k)+2*pre_ind_sd; %upper boundary (2 times SD)
%                     
%                     %importing time series
%                     series_regress=ans2{l,c}; %load timeseries
%  
%                     %% if series dives into the envelope from above upper boundary
%                     if PeakVALUE(l,c,k)>PI_ub %series_regress(140)>PI_ub 
%                         
%                         %finding whether there is a first point after the start
%                         %of CO2atm ramp down to enter time series from above
%                         point_envelope_up=find(series_regress(TpeakMAX(l,c,k):end)<PI_ub);
%                         noreturn2up=[]; trecovery_dummy_up=[];
%                         
%                         if ~isempty(point_envelope_up) %in case there are points that enter the envelope take the index difference to find the first 10 consecutive points
%                             dummy2up=diff(point_envelope_up)==1;
%                             M2up = movsum(dummy2up,[0 gap]);
%                             noreturn2up=find(M2up==10,1,'first');
%                             
%                             if ~isempty(noreturn2up) %if there are 10 consecutive points
%                                 trecovery_dummy_up=point_envelope_up(noreturn2up)+MitigationStart-1;
%                                 Trecovery(l,c,k)=trecovery_dummy_up;
%                                 
%                             else  %if there are no 10 consecutive points then do regression
%                                   %calculate polyfit based on extension ramp points
%                                 
%                                 curve=fit(time_recovery_regress,series_regress(end-reg_window:end)','poly1'); %polynomial fit
%                                 polinomyal_up=[curve.p1 curve.p2-PI_ub]; %coefficients
%                                 trecovery_dummy_up=round(roots(polinomyal_up)); %solving and rounding polynomial roots
%                                 Trecovery(l,c,k)=trecovery_dummy_up; %Trec 
%                                 
%                                 %EXPLANATION FOR A 3rd ORDER POLYNOM 
%                                 %curve=fit(time_recovery_regress,series_regress(end-reg_window:end)','poly3');
%                                 
%                                 %equation for any Y is given by the coeficients of the curve
%                                 %Y=(curve.p1)*x^3+ curve.p2)*x^2+ (curve.p3)*x+ curve.p4 ;
%                                 %if I want to solve for Y=PI_lwb I just need to addapt the equation
%                                 %PI_lwb = (curve.p1*x)^3+ (curve.p2*x)^2+ curve.p3*x+ curve.p4;
%                                 %0 = (curve.p1*x)^3+ (curve.p2*x)^2+ curve.p3*x+ curve.p4 -PI_lwb
%                                 %then the 'roots' would correspond to an Y=PI_lwb
%                                 % pol=[curve.p1 curve.p2 curve.p3 curve.p4];
%                                 % polyval(polinomial,140:700);               
%                                 %
%                                 % polinomyal_up=[curve.p1 curve.p2 curve.p3 curve.p4-PI_ub]; %specifying the coefficients 'note the last coeficient incorporates PI_lwb
%                                 % dummyroots_up=roots(polinomyal_up); %finding the roots (in this case 3) and picking the earliest one which its imaginary part is zero
%                                 %
%                                 % [indxroot_up]=find(imag(dummyroots_up)==0); %get only the roots where imaginary part is zero.
%                                 % dummyroots_up=round(dummyroots_up(indxroot_up));
%                                 %
%                                 % trecovery_dummy_up=min(dummyroots_up(dummyroots_up>200 & dummyroots_up<1000)); %then take only the real part (because it's still real+000i), round it and then take the earliest time between 280 and 1000
%                                 % Trecovery(l,c,k)=trecovery_dummy_up; %this finds the intercept after the start of the extension and it has to bee less than 1000yrs otherwise no recovery
%                                 
%                             end
%                             
%                         elseif isempty(point_envelope_up) || isempty(noreturn2up) %if there are no points to enter envelope do a regression using reg window of the series anyway
%                             curve=fit(time_recovery_regress,series_regress(end-reg_window:end)','poly1');
%                             polinomyal_up=[curve.p1 curve.p2-PI_ub];
%                             trecovery_dummy_up=round(roots(polinomyal_up));
%                             Trecovery(l,c,k)=trecovery_dummy_up;
%                         end
%                     end
%                     
%                     
%                     %% if series dives into the envelope from lower boundary
%                     if PeakVALUE(l,c,k)<PI_lwb %PI_ub series_regress(140)<PI_lwb
%                         %finding whether there is a first point after the start
%                         %of CO2atm ramp down to enter time series from above
%                         point_envelope_down=find(series_regress(TpeakMAX(l,c,k):end)>PI_lwb);
%                         noreturn2down=[]; trecovery_dummy_down=[];
%                         
%                         if ~isempty(point_envelope_down) %in case there are points that enter the envelope take the index difference to find the first 10 consecutive points
%                             dummy2down=diff(point_envelope_down)==1;
%                             M2down = movsum(dummy2down,[0 gap]);
%                             noreturn2down=find(M2down==10,1,'first');
%                             
%                             if ~isempty(noreturn2down) %if there are 10 consecutive points
%                                 trecovery_dummy_down=point_envelope_down(noreturn2down)+MitigationStart-1;
%                                 Trecovery(l,c,k)=trecovery_dummy_down;
%                                 
%                             else%if there are no 10 consecutive points then do regression
%                                 %calculate polyfit based on extension ramp points
%                                 curve=fit(time_recovery_regress,series_regress(end-reg_window:end)','poly1');
%                                 polinomyal_down=[curve.p1 curve.p2-PI_lwb];
%                                 trecovery_dummy_down=round(roots(polinomyal_down));
%                                 Trecovery(l,c,k)=trecovery_dummy_down;  
%                             end
%                             
%                         elseif isempty(point_envelope_down) || isempty(noreturn2down) %if there are no points to enter envelope do a regression using reg window of the series anyway
% 
%                             curve=fit(time_recovery_regress,series_regress(end-reg_window:end)','poly1');
%                             polinomyal_down=[curve.p1 curve.p2-PI_lwb];
%                             trecovery_dummy_down=round(roots(polinomyal_down));
%                             Trecovery(l,c,k)=trecovery_dummy_down;
%                             
%                         end
%                     end
%                     
% %                     if series_regress(140)>PI_lwb && series_regress(140)<PI_ub %if already inside the envelope then early recovery
% %                         Trecovery(l,c,k)=-4444;
% %                     end
%                     
%                     
%                 elseif domain_mask(l,c)==2 && isnan(PI_mean(l,c,k)) %if is atlantic but PI mean does not exist
%                     Trecovery(l,c,k)=bottomflag;
%                 end
%             end
%         end
%     end
    
    sprintf('Saving all results')
    save(sprintf('/Volumes/LaCie_Leonardo/NorESM/all_ramps/filtered/Results/%s_results_updated_v1D_poly1_NEW_NEW.mat',varlist{varIDX}), 'Trecovery', 'TOD', 'TimeLAG', 'PeakVALUE', 'TpeakMAX')
    
end
toc

 %% flagging
                    
%if trecovery_dummy happens before mitigation start
%then flag early
%                     if trecovery_dummy<MitigationStart+gap %%%%REMEBER
%                         Trecovery(l,c,k)=badflagtooearly;
%                     end
%
%if regression time turns out to be higher
%than 1000 yrs AND slope of regression is
%positive it means very slow recovery
%                     if  trecovery_dummy>1000 && p1>0
%                         Trecovery(l,c,k)=badflagslow;
%                     end

%if regression curve has negative slope,
%means it will never recover.
%                     if isempty(trecovery_dummy)
%                         Trecovery(l,c,k)=badflagsnegslope;
%                     elseif trecovery_dummy<MitigationStart
%                      Trecovery(l,c,k)=badflagtooearly;
%                     else
%                         Trecovery(l,c,k)=trecovery_dummy;
%                     end



end
