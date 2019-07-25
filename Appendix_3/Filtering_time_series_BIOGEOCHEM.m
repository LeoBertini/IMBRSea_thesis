
%%%THIS IS NOW SET FOR BIOGEOCHEM VARIABLES
tic
clear all
close all
warning('off','all')

%% Directories

%cd /shared/projects/uniklima/globclim/nun008/NorESM/ocn/
cd /Volumes/LaCie_Leonardo/NorESM/all_ramps

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

%% Loading domain mask

domain_mask=load ('/Volumes/LaCie_Leonardo/NorESM/Mask_regions/Mask_corrected.mat');
domain_mask=domain_mask.domain_mask;


%% loading  depths

depth=load('/Volumes/LaCie_Leonardo/NorESM/all_ramps/Depth_Levels.mat'); %get depth array from one of the files in the folder
depth=depth.depths;


%% VARIABLE LOOP

varlist={'ph';'o2';'omegac';'detoc';'POtracer'};

for varIDX=5 %1:length(varlist) %this is set for all 

%% Importing Pre-Industrial 30-yr mean

PI_mean=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/Pre_industrial/30yrMeans/PI_%s_30yr_mean.mat',varlist{varIDX}));
PI_mean=PI_mean.PI_30mean;
PI_mean(PI_mean==0)=NaN;

PI_mean_bit1=PI_mean(200:end,:,:); %reshaping -- position 200 is where Pacific is
PI_mean_bit2=PI_mean(1:199,:,:);
PI_mean=[PI_mean_bit1 ; PI_mean_bit2];

%% Importing Pre-Industrial 30-yr Standard Deviation

PI_sd=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/Pre_industrial/30yrMeans/PI_%s_30yr_sd.mat',varlist{varIDX}));
PI_sd=PI_sd.PI_sd;

PI_sd_bit1=PI_sd(200:end,:,:); %reshaping -- position 200 is where Pacific is
PI_sd_bit2=PI_sd(1:199,:,:);
PI_sd=[PI_sd_bit1 ; PI_sd_bit2];

clear *bit*

% figure
% pcolor(PI_mean(:,:,55)')
% shading flat
% colorbar


%% IMPORTING VARIABLES FROM RAMP FILES

if strcmp(varlist{varIDX},'ph')==1
    ncvarname='phlvl';
elseif strcmp(varlist{varIDX},'o2')==1
    ncvarname='o2lvl';
elseif strcmp(varlist{varIDX},'omegac')==1
    ncvarname='omegaclvl';
elseif strcmp(varlist{varIDX},'detoc')==1
    ncvarname='detoclvl';
elseif strcmp(varlist{varIDX},'POtracer')==1
end

for i=1:length(filenames_ramps)
    i
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
    
    if ~(strcmp(varlist{varIDX},'POtracer')==1)
    DataRamp=ncread(filenames_ramps(i,:),ncvarname); %importing ramps
    DataRamp_bit1=DataRamp(200:end,:,:); %reshaping -- position 200 is where Pacific is
    DataRamp_bit2=DataRamp(1:199,:,:);
    DataRamp=[DataRamp_bit1 ; DataRamp_bit2];
    
    else
    o2=ncread(filenames_ramps(i,:),'o2lvl');
    po4=ncread(filenames_ramps(i,:),'po4lvl');
    DataRamp=172*po4+o2; %importing ramps
    DataRamp_bit1=DataRamp(200:end,:,:); %reshaping -- position 200 is where Pacific is
    DataRamp_bit2=DataRamp(1:199,:,:);
    DataRamp=[DataRamp_bit1 ; DataRamp_bit2];
    end
    
    
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

    if ~(strcmp(varlist{varIDX},'POtracer')==1)
    DataPostRamp=ncread(filenames_extension(i,:), ncvarname); %importing ramps
    DataPostRamp1=DataPostRamp(200:end,:,:); %reshaping -- position 200 is where Pacific is
    DataPostRamp2=DataPostRamp(1:199,:,:);
    DataPostRamp=[DataPostRamp1 ; DataPostRamp2];
    DataPostRamp=permute(DataPostRamp,[2,1,3]);
    
    else
    o2=ncread(filenames_extension(i,:),'o2lvl');
    po4=ncread(filenames_extension(i,:),'po4lvl');
    DataPostRamp=172*po4+o2; %importing ramps
    DataPostRamp1=DataPostRamp(200:end,:,:); %reshaping -- position 200 is where Pacific is
    DataPostRamp2=DataPostRamp(1:199,:,:);
    DataPostRamp=[DataPostRamp1 ; DataPostRamp2];
    DataPostRamp=permute(DataPostRamp,[2,1,3]);
    end
     
    %saving everything in a structure
    PostRampStruk(i)=struct('FileName', filenames_extension(i),'YearINDEX',YearINDEXpost,'YearStrpost',YearStrpost, 'Field_PostRamp', DataPostRamp); %saving all variables into structure each year being a column
end

clearvars -except  PI_mean PI_sd  ...
    RampStruk PostRampStruk domain_mask...
    k depth filenames_ramps filenames_extension ncvarname varlist varIDX


%% tranposing variables from (320X384 to 384X320) to match domain_mask

PI_mean=permute(PI_mean,[2,1,3]);
PI_sd=permute(PI_sd,[2,1,3]);
toc

% 500 sec up to this point

    time_lowest_regress=[RampStruk.YearINDEX PostRampStruk.YearINDEX]; %concatenating all the years to use in the regression
    time_lowest_regress=str2num(datestr(time_lowest_regress,'YYYY'))-str2num(datestr(RampStruk(1).YearINDEX,'YYYY')); %transforming them into positive integers
    time_lowest_regress=time_lowest_regress+1;
    
%     %Saving Tlowest regress to be used in the Time-scales script
%     fprintf('Saving time-array for regression \n')
%     save('/Volumes/LaCie_Leonardo/NorESM/all_ramps/time_array.mat','time_lowest_regress')

%% Filtering all the series

%all_new_series=cell(384,320,70); %this is where all the filtered series for all depths are gonna be stored


for  k=1:length(depth)
    
    fprintf('Depth Level=%1d \n',k);
    
    %declaring varibales to be used inside loop
    
    series_interpolated=zeros;
    new_series=cell(384,320); %this is the individual series for depth k
    gap=9;
    
    
    for l=1:size(domain_mask,1) %setting latitudinal limits
        for c=1:size(domain_mask,2) %setting longitudinal limits
            
            %if domain_mask(l,c)==2 && ~isnan(PI_mean(l,c,k)) %if it's Atlantic and PI_pH exists
             
            %including arctic and southern ocean
            if (domain_mask(l,c)==2 || domain_mask(l,c)==8) && ~isnan(PI_mean(l,c,k)) 
   
            
            %importing time series %RAMP VALUES
                for i=1:length(RampStruk)
                    if ~isnan(RampStruk(i).Field_Ramp(l,c,k)) %if Ramp Value exists
                        RAMPlinevalues(i)=RampStruk(i).Field_Ramp(l,c,k); %load timeseries
                    end
                end
                
                %POST RAMP VALUES
                for n=1:length(PostRampStruk)
                    if ~isnan(PostRampStruk(n).Field_PostRamp(l,c,k)) %If Post ramp value exists
                        POSTRAMPlinevalues(n)=PostRampStruk(n).Field_PostRamp(l,c,k);  %load timeseries
                    end
                end
                
                series_original=[RAMPlinevalues POSTRAMPlinevalues]; %concatenating ramps and post ramps
                series_interpolated=series_original; %making a copy of the original series to interpolate
                %hampel filter
                hampelout=hampel(series_original);
                %removing remaining spikes
                final_fitered=medfilt1(hampelout,3);
                final_fitered=medfilt1(final_fitered,3);
                
                %
                %                     figure
                %                     plot(seriesPH_original,'linewidth',2)
                %                     hold on
                %                     plot(final_fitered,'linewidth',2)
                %                     hold on
                %                     legend ('Original pH series','MedianFilter(Hampel)')
                %                     title(sprintf ('line=%d column=%d',l,c))
                
                difference=series_original-final_fitered;
                absolut_diff=abs(difference);
                threshold=2*PI_sd(l,c,k);
                
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
                
                
%                 
%                                     figure
%                                     plot(series_original,'linewidth',3)
%                                     hold on
%                                     plot(series_interpolated,'linewidth',2)
%                                     hold on
%                                     plot(final_fitered,'linewidth',2)
%                                     legend ('Original pH series', 'Linear Interp(Original Series)', 'MedianFilter(Hampel)')
%                                     title(sprintf ('line=%d column=%d',l,c))
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
    
end

fprintf('\n');
fprintf('####All individual series have been filtered and saved for variable %s ###',varlist{varIDX});
fprintf('\n');

clear new_series

end

fprintf('Initialize script for calculation of time variables \n');
