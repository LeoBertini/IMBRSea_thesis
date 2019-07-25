%% EXPORT PRODUCTION PRE-IUNDUSTRIAL

clear all 
close all

   %variables in biogechem ncfiles are ...
    ... o2=ncread(FN, 'o2'); %units = 'mol O2 m-3' --> to convert to ?mol/kg --> (o2/sigma)*E6 
    ... thickness is 'pddpo', 
    ... sigma= ncread(FN, 'sigma'); %units = 'kg m-3'
    ... omgegac (calcite saturation state)
    ... detoc=ncread(FN,'detoc'); %units = 'mol C m-3' --> detritus
    ... export production at 100 m = ncread(FN,'epc100'); %units ---> mol C m-2 s-1 

%path=['/shared/projects/uniklima/globclim/nun008/NorESM/PI/'];

%##TO USE IN GRUNCH#
% cd ('/shared/projects/uniklima/globclim/nun008/NorESM/PI/')
% listings=dir(['/shared/projects/uniklima/globclim/nun008/NorESM/PI/' 'N*.nc']);
% load('/shared/projects/uniklima/globclim/nun008/NorESM/PI/Depth_Levels.mat')

%##TO USE WITH PERSONAL HARD DRIVE#
cd('/Volumes/LaCie_Leonardo/NorESM/Pre_industrial')
listings=dir(['/Volumes/LaCie_Leonardo/NorESM/Pre_industrial/' 'N*.nc']);
depth=load('/Volumes/LaCie_Leonardo/NorESM/all_ramps/Depth_Levels.mat');
depths=depth.depths;

%%
varlist={'epc100'}; %list of variables to have their data interpolated and then monthly data reduced to annual

%%
filenames_org = {listings.name};
filenames_org= filenames_org';

jump=11; %this is the 'jump window of 'jump+1'--> reads a sequence of 12 files
indx_start=1; %this is where we indicate the first file that's going to be read 
indx_end=indx_start+jump;
times=1;

varIDX=1;


while times<=30 %the amount of years
        
        
        for a=indx_start:indx_end %--> here the loop for all the .nc files starts with a and ends with b
            
            FN=['/Volumes/LaCie_Leonardo/NorESM/Pre_industrial/' char(filenames_org(a))]; %converts into characther array with complete path of target nc file to be open
            
            %FN=['/shared/projects/uniklima/globclim/nun008/NorESM/PI/' char(filenames_org(a))]; %converts into characther array with complete path of target nc file to be open
            %ncdisp(FN);
            
            %show file number
            fprintf('File number=%.f\n',a)
            
            %importing variables from nc file
            
            VarField=ncread(FN, varlist{varIDX});
            
            %figure
            %pcolor(VarField');shading flat
            
            
            final_struc.File_Name{a}=listings(a).name;
            final_struc.Field{a}=VarField;
            
        end
             
            fprintf('All montlhy files for year %d loaded in structure \n',times)
            fprintf('\n')
            
        
            
        storage_struc=final_struc.Field(indx_start:end);    
            
        %computing the means
        fprintf('Calculating file of annual mean\n')
        
        mean_test=zeros(320,384); %allocating space
        
        tic
        for i=1:320
            for j=1:384
               
                    for h=1:length(storage_struc)  %for each file in the cell array
                        
                        %pH
                        test(h)=storage_struc{h}(i,j);
                        mean_test(i,j)=nanmean(test(:));
                        
                    end

            end
        end
        toc
        
%         figure
%         pcolor(mean_test);shading flat


        elapsedTime = toc;
        fprintf('Elapsed time = %.2f min \n',elapsedTime/60);
        
        name1=filenames_org(indx_start);
        name1=name1{1,1}(42:45);
        name2=filenames_org(indx_end);
        name2=name2{1,1}(42:45);
        
        name_file2=sprintf('/Volumes/LaCie_Leonardo/NorESM/Pre_industrial/Partials/mean_PI_%s_%s_%s.mat',varlist{varIDX},name1, name2);
        save (name_file2,'mean_test', '-v7.3')
        
        indx_start=indx_end+1; %this is where we indicate the first file that's going to be read in the next loop
        indx_end=indx_start+jump;
        times=times+1;
        
        clear storage_struc final_struc
        
    end

%% %EXPORT PRODUCTION PRE-IUNDUSTRIAL TRENDS

cd /Volumes/LaCie_Leonardo/NorESM/Pre_industrial/Partials/

folder='/Volumes/LaCie_Leonardo/NorESM/Pre_industrial/Partials/';
    
    listings=dir(sprintf('/Volumes/LaCie_Leonardo/NorESM/Pre_industrial/Partials/mean_PI_%s_*.mat',varlist{varIDX})); %list of all partial files
    
    filenames_org = {listings.name};
    filenames_org= filenames_org';
    
%     
%     for a=length(filenames_org)-4:length(filenames_org) %length(filenames_org) (last 5 yrs of PI)
%         partial_avg(length(filenames_org)-a+1)=load(char(filenames_org{a,1})); %each mat file is loaded into this structure
%     end
%     
    
    
    for a=1:5 %length(filenames_org) (first 5 yrs of PI)
        partial_avg(a)=load(char(filenames_org{a,1})); %each mat file is loaded into this structure
    end
    
 %computing 5yr mean
    
    PI_mean=zeros(320,384);
    whos
    for l=1:320
        for c=1:384
            
            %sprintf('%s line %d column %d',varlist{varIDX},l,c)
               
                for i=1:length(partial_avg)
                    seriestoaverage(i)=partial_avg(i).mean_test(l,c);
                end
                
                 PI_mean(l,c)= nanmean(seriestoaverage);
                 %this is the 5yr mean to be ploted as center of envelope
        end
    end
    
    PI_mean(PI_mean==0)=NaN;
    
%     figure
%     pcolor(PI_mean');shading flat
    
    fprintf('Saving mean based on first 5yr of PI  \n')
    fprintf('\n')
    save (sprintf('/Volumes/LaCie_Leonardo/NorESM/Pre_industrial/30yrMeans/PI_%s_first_5yr_mean.mat',varlist{varIDX}), 'PI_mean', '-v7.3')
    
%30yr trends
     clear a b partial_avg
    
    fprintf('Now Calculating detrending and 30yr mean and 30yr sd \n')
    
    
    for a=1:length(filenames_org) %all the 30 files of the PI
        partial_avg(a)=load(char(filenames_org{a,1})); %each mat file is loaded into this structure
    end
    
    %computing time series detrending for each grid cell
    
    PI_detrended=zeros(320,384);
    PI_detrended=num2cell(PI_detrended);
    PI_sd=zeros(320,384);
    PI_30mean=zeros(320,384);
    seriestofilter=NaN(1,30);
    
    whos

    
    for l=1:320
        for c=1:384
            
           sprintf('%s line %d column %d',varlist{varIDX},l,c)
            
                for i=1:length(partial_avg)
                    seriestofilter(i)=partial_avg(i).mean_test(l,c);
                end
                
                PI_sd(l,c)=nanstd(seriestofilter);
                [X, T]=detrend_NAN(seriestofilter); % X is the increments to be added or removed to form a detrended series
                % T is the actual trend line
                PI_30mean(l,c)=nanmean(seriestofilter);
                PI_detrended(l,c)= {X + PI_30mean(l,c)} ; %this is the 30yr mean to be added back to the detrended series to get corrected for time trend
                
%                 figure
%                 plot(seriestofilter)
%                 hold on
%                 plot(PI_detrended{l,c})
%                 hold on
%                 plot(T)
%                 legend ('Original','Detrended data','Trendline')
                
        end
    end
    
    PI_sd(PI_sd==0)=NaN;
    
    fprintf('Now Saving Detrended data and 30yr mean and 30 yr sd \n')
    save (sprintf('/Volumes/LaCie_Leonardo/NorESM/Pre_industrial/30yrMeans/PI_%s_30yr_detrended.mat',varlist{varIDX}), 'PI_detrended', '-v7.3')
    save (sprintf('/Volumes/LaCie_Leonardo/NorESM/Pre_industrial/30yrMeans/PI_%s_30yr_mean.mat',varlist{varIDX}), 'PI_30mean', '-v7.3')
    save (sprintf('/Volumes/LaCie_Leonardo/NorESM/Pre_industrial/30yrMeans/PI_%s_30yr_sd.mat',varlist{varIDX}), 'PI_sd', '-v7.3')
    clear partial_avg
    

%% FILTERING CARBON EXPORT 100m

cd /Volumes/LaCie_Leonardo/NorESM/all_ramps

%% Loading domain mask

domain_mask=load ('/Volumes/LaCie_Leonardo/NorESM/Mask_regions/Mask_corrected.mat');
domain_mask=domain_mask.domain_mask;


%% List of RAMP files

listings_ramps= dir('/Volumes/LaCie_Leonardo/NorESM/all_ramps/N*.nc'); %list all nc files of interest
filenames_ramps= string({listings_ramps.name});
filenames_ramps='/Volumes/LaCie_Leonardo/NorESM/all_ramps/'+ filenames_ramps;
filenames_ramps=filenames_ramps';  %converts into string array

% arranging into ramp up then ramp down
filenames_ramps=[filenames_ramps(141:end);filenames_ramps(1:140)];

%% List of POST-RAMP files

listings_extension= dir('/Volumes/LaCie_Leonardo/NorESM/post_ramps/Annually_Biogeochem/*.nc'); %list all nc files of interest
filenames_extension_bit1 =string({listings_extension.folder});
filenames_extension_bit2= string({listings_extension.name});
filenames_extension =filenames_extension_bit1+'/'+filenames_extension_bit2;
filenames_extension=filenames_extension';
clear filenames_extension_bit*
    
% Importing Pre-Industrial 30-yr mean

PI_mean=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/Pre_industrial/30yrMeans/PI_%s_30yr_mean.mat',varlist{varIDX}));
PI_mean=PI_mean.PI_30mean;

PI_mean_bit1=PI_mean(200:end,:,:); %reshaping -- position 200 is where Pacific is
PI_mean_bit2=PI_mean(1:199,:,:);
PI_mean=[PI_mean_bit1 ; PI_mean_bit2];

% Importing Pre-Industrial 30-yr Standard Deviation

PI_sd=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/Pre_industrial/30yrMeans/PI_%s_30yr_sd.mat',varlist{varIDX}));
PI_sd=PI_sd.PI_sd;

PI_sd_bit1=PI_sd(200:end,:,:); %reshaping -- position 200 is where Pacific is
PI_sd_bit2=PI_sd(1:199,:,:);
PI_sd=[PI_sd_bit1 ; PI_sd_bit2];

clear *bit*

% figure
% pcolor(PI_mean(:,:)')
% shading flat
% colorbar

ncvarname='epc100';

for i=1:length(filenames_ramps)
    
    %ncdisp(filenames_ramps(i,:))
    
    if i<=140 %for the ramp up
        timeimport=ncread(filenames_ramps(i,:), 'time'); %importing var
        startdate=datenum(datetime(0001,01,01)); %days since 0001-01-01 00:00
        time(i)= addtodate(startdate,timeimport,'day');
        
    else %for the ramp down
        timeimport=ncread(filenames_ramps(i,:), 'time'); %importing var
        startdate=datenum(datetime(0141,01,01)); %days since 0141-01-01 00:00'
        time(i)= addtodate(startdate,timeimport,'day');
    end
    
    
    DataRamp=ncread(filenames_ramps(i,:),ncvarname); %importing ramps
    DataRamp_bit1=DataRamp(200:end,:,:); %reshaping -- position 200 is where Pacific is
    DataRamp_bit2=DataRamp(1:199,:,:);
    DataRamp=[DataRamp_bit1 ; DataRamp_bit2];
    
    DataRamp=permute(DataRamp,[2,1,3]); %this is to make the data orientation the same as the grid mask
    
    %saving everything in a structure
    RampStruk(i)=struct('FileName',filenames_ramps(i,:),'YearINDEX',time(i),'YearStr',datestr(time(i)),'Field_Ramp',DataRamp); %saving all variables into structure each year being a column
end

%% IMPORTING POST_RAMP_VALUES

for i=1:length(filenames_extension)
    
    %ncdisp(filenames_extension(i,:))
    
    
    timeimport=ncread(filenames_extension(i,:), 'time'); %importing var
    startdate=datenum(datetime(0280,01,01));   %days since 0280-01-01 00:00
    
    YearINDEXpost= addtodate(startdate,timeimport,'day');
    YearStrpost=datestr(YearINDEXpost);
    
    DataPostRamp=ncread(filenames_extension(i,:), ncvarname); %importing ramps
    DataPostRamp1=DataPostRamp(200:end,:,:); %reshaping -- position 200 is where Pacific is
    DataPostRamp2=DataPostRamp(1:199,:,:);
    DataPostRamp=[DataPostRamp1 ; DataPostRamp2];
    DataPostRamp=permute(DataPostRamp,[2,1,3]);
    
    %saving everything in a structure
    PostRampStruk(i)=struct('FileName', filenames_extension(i),'YearINDEX',YearINDEXpost,'YearStrpost',YearStrpost, 'Field_PostRamp', DataPostRamp); %saving all variables into structure each year being a column
end

clearvars -except  PI_mean PI_sd  ...
    RampStruk PostRampStruk domain_mask...
    k depth filenames_ramps filenames_extension ncvarname varlist varIDX

PI_mean=permute(PI_mean,[2,1,3]);
PI_sd=permute(PI_sd,[2,1,3]);   

time_lowest_regress=[RampStruk.YearINDEX PostRampStruk.YearINDEX]; %concatenating all the years to use in the regression
    time_lowest_regress=str2num(datestr(time_lowest_regress,'YYYY'))-str2num(datestr(RampStruk(1).YearINDEX,'YYYY')); %transforming them into positive integers
    time_lowest_regress=time_lowest_regress+1;

    %% Filtering all the series

    k=100;
    fprintf('Depth Level=100m \n');
    
    %declaring varibales to be used inside loop
    
    series_interpolated=zeros;
    new_series=cell(384,320); %this is the individual series for depth k
    gap=9;

    for l=1:size(domain_mask,1) %setting latitudinal limits
        for c=1:size(domain_mask,2) %setting longitudinal limits
            
            if domain_mask(l,c)==2 && ~isnan(PI_mean(l,c)) %if it's Atlantic and PI_pH exists
                
                %importing time series %RAMP VALUES
                for i=1:length(RampStruk)
                    if ~isnan(RampStruk(i).Field_Ramp(l,c)) %if Ramp Value exists
                        RAMPlinevalues(i)=RampStruk(i).Field_Ramp(l,c); %load timeseries
                    end
                end
                
                %POST RAMP VALUES
                for n=1:length(PostRampStruk)
                    if ~isnan(PostRampStruk(n).Field_PostRamp(l,c)) %If Post ramp value exists
                        POSTRAMPlinevalues(n)=PostRampStruk(n).Field_PostRamp(l,c);  %load timeseries
                    end
                end
                
                series_original=[RAMPlinevalues POSTRAMPlinevalues]; %concatenating ramps and post ramps
                series_interpolated=series_original; %making a copy of the original series to interpolate
                %hampel filter
                hampelout=hampel(series_original);
                %removing remaining spikes
                final_fitered=medfilt1(hampelout,3);
                final_fitered=medfilt1(final_fitered,3);
                
                
%                                     figure
%                                     plot(series_original,'linewidth',2)
%                                     hold on
%                                     plot(final_fitered,'linewidth',2)
%                                     hold on
%                                     legend ('Original pH series','MedianFilter(Hampel)')
%                                     title(sprintf ('line=%d column=%d',l,c))
%                                     ylim([0 20*10^(-8)])

%                 
                difference=series_original-final_fitered;
                absolut_diff=abs(difference);
                threshold=2*PI_sd(l,c);
                
                indexes=find(absolut_diff>threshold); %points in the series where the difference is higher than threshold
                
                %declaring variables before loop
                jump=1; %initial window
                
                idx=1;
                loop_number=0;
                
                while loop_number<3
                    
                    for idx=1:length(indexes)
                        
                        if indexes(idx)~=1 && indexes(idx)~=480
                            
                            jump=1; %pula mais um pra frente e um pra traz
                            back=(indexes(idx)-1-jump);
                            forth=(indexes(idx)+1+jump);
                            
                            %enquanto back e forth forem elementos de
                            %indexes, continue atualizando jump e back e
                            %forth em seguida
                            
                            
                            while ismember(back,indexes) && ismember(forth,indexes)
                                jump=jump+1;
                                back=(indexes(idx)-1-jump);
                                forth=(indexes(idx)+1+jump);
                            end
                            
                            %calcula a interpolacao agora espelhada se back nao for <=0
                            if back>0  && forth<480
                                windowY=[series_interpolated(indexes(idx)-1-jump) series_interpolated(indexes(idx)+1+jump)];
                                windowX=[time_lowest_regress((indexes(idx)-1-jump)) time_lowest_regress((indexes(idx)+1+jump))];
                                series_interpolated(indexes(idx))=interp1(windowX,windowY,time_lowest_regress(indexes(idx)));
                                
                                difference(indexes(idx))=series_interpolated(indexes(idx))-final_fitered(indexes(idx)); %difference between original and interpolated
                                absolut_diff(indexes(idx))=abs(difference(indexes(idx)));
                                
                            else
                                continue
                            end
                            
                            
                        else
                            
                            if indexes(idx)==1 && length(idx)>idx %caso o elemento seja o primeiro ou o ultimo da time series
                                idx=idx+1;
                            else
                                break
                                
                                if indexes(idx)==480
                                    idx=idx-1;
                                end
                                
                                
                            end
                        end
                    end
                    
                    loop_number=loop_number+1;
                    indexes=find(absolut_diff>threshold);
                end
                
                
                
                
                %For point where the filter could not remove, just susbstitute using the previous values from the hampel filter
                indexes=find(absolut_diff>threshold); %points in the series where the difference is higher than threshold
                series_interpolated(indexes)=final_fitered(indexes);
                series_interpolated=hampel(series_interpolated);
                new_series(l,c)={series_interpolated}; %saving all the interpolated data
                
                
                
%                                     figure
%                                     plot(series_original,'linewidth',3)
%                                     hold on
%                                     plot(series_interpolated,'linewidth',2)
%                                     hold on
%                                     plot(final_fitered,'linewidth',2)
%                                     legend ('Original pH series', 'Linear Interp(Original Series)', 'MedianFilter(Hampel)')
%                                     title(sprintf ('line=%d column=%d',l,c))
%                                     ylim([0 20*10^(-8)])
%                                     
%                                     
%                                     figure
%                                     pcolor(PI_mean(:,:,k));shading flat
%                                     hold on
%                                     scatter(c,l,'filled','red')
%                                     set(gca, 'color','k' );
                                      %axis([50 250 250 350])
                                
                                    
%                                     figure
%                                     addpath     /Volumes/LaCie_Leonardo/NorESM
%                                     field1=permute(PI_mean(:,:,k),[3,1,2]);
% 
%                                     %unflipping
%                                     bit1=field1(:,:,120:end); %bit1 (200:end,:,:); %reshaping -- position 200 is where Pacific is
%                                     bit2=field1(:,:,1:119);   %bit2 (1:199,:,:);
%                                     field2=cat(3,bit1,bit2);
%                                     field3=squeeze(field2);
%                                     
%                                     m_proj('Lambert','lon',[-100 20], 'lat', [-10 75])
%                                     m_pcolor(shift_atl(plon),shift_atl(plat),shift_atl(field3))
%                                     caxis([7.8 8.2])
%                                     m_coast('patch',[.7 .7 .7],'edgecolor','k');
%                                     m_grid('backgroundcolor','k','box','fancy','tickdir','out'); %if brown [160/255 82/255 45/255]
%                                     hold on
%                                     h2=m_line(plon(c,l),plat(c,l),'marker','o','color','r','linewi',2,...
%                                         'linest','none','markersize',8,'markerfacecolor','w');
% 


                                    
                                    
                
            else
                new_series(l,c)={NaN};
                
            end
        end
    end
    
    %all_new_series(:,:,k)=new_series;
    %saving individual filtered series
    fprintf('saving individual filtered series for Depth Level=%1d \n',k);
    save(sprintf('/Volumes/LaCie_Leonardo/NorESM/all_ramps/filtered/new_series_%s_k_%d.mat',varlist{varIDX},k),'new_series')
    

fprintf('\n');
fprintf('####All individual series have been filtered and saved for variable %s ###',varlist{varIDX});
fprintf('\n');

% figure
% plot(new_series{l,c})
% ylim([0 20*10^(-8)])

    
