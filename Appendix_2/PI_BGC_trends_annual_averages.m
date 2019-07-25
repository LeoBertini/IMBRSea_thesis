clear all
close all

%THIS IS SET FOR the following variables
variablelist={'ph';'o2';'omegac';'detoc';'AOU';'templvl';'salnlvl'};

cd /media/nun008/Seagate/Leo_Files/PI_NorESM/Pre_Industrial_Means

for varIDX=4:5%:length(variablelist) %starting with all but salinity this time
    
    %path=['/Volumes/LaCie_Leonardo/NorESM/Pre_industrial/'];
    
    folder='/media/nun008/Seagate/Leo_Files/PI_NorESM/Pre_Industrial_Means';
    
    listings=dir(sprintf('/media/nun008/Seagate/Leo_Files/PI_NorESM/Pre_Industrial_Means/mean_PI_%s_*.mat',variablelist{varIDX})); %list of all partial files
    
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
    
    
    
    
    % loading depths
    depth=load('/media/nun008/Seagate/Leo_Files/PI_NorESM/Pre_Industrial_Means/Depth_Levels.mat'); %get depth array from one of the files in the folder
    depth=depth.depths;
    
    
    %computing 5yr mean
    
    PI_mean=zeros(320,384,70);
    whos
    
    for l=1:320
        for c=1:384
            
            sprintf('%s line %d column %d',variablelist{varIDX},l,c)
            
            for k=1:length(depth)
                
                
                for i=1:length(partial_avg)
                    seriestoaverage(i)=partial_avg(i).mean_test(l,c,k);
                end
                
                 PI_mean(l,c,k)= nanmean(seriestoaverage);
                 %this is the 5yr mean to be ploted as center of envelope
                
            end
        end
    end
    
    
    % %     figure
    % %     plot(seriestofilter)
    % %     hold on
    % %     plot(hampelout)
    
    % %                      figure
    % %                      plot(test)
    % %                      hold on
    % %                      plot(T)
    % %                      hold on
    % %                      plot(detrended_series{l,c,k},'LineWidth',2)
    % %                      legend('Original','Trend line', 'Detrended')
    %
    
    PI_mean(PI_mean==0)=NaN;
    
    fprintf('Saving mean based on first 5yr of PI  \n')
    fprintf('\n')
    save (sprintf('/media/nun008/Seagate/Leo_Files/PI_NorESM/Pre_Industrial_Means/30yrMeans/PI_%s_first_5yr_mean.mat',variablelist{varIDX}), 'PI_mean', '-v7.3')
    
    %Calculating detrending and 30yr sd
    clear a b partial_avg
    
    fprintf('Now Calculating detrending and 30yr mean and 30yr sd \n')
    
    
    for a=1:length(filenames_org) %all the 30 files of the PI
        partial_avg(a)=load(char(filenames_org{a,1})); %each mat file is loaded into this structure
    end
    
    
    %computing time series detrending for each grid cell
    
    PI_detrended=zeros(320,384,70);
    PI_detrended=num2cell(PI_detrended);
    PI_sd=zeros(320,384,70);
    PI_30mean=zeros(320,384,70);
    seriestofilter=NaN(1,30);
    
    whos

    
    for l=1:320
        for c=1:384
            
            sprintf('%s line %d column %d',variablelist{varIDX},l,c)
            
            for k=1:length(depth)
   
                for i=1:length(partial_avg)
                    seriestofilter(i)=partial_avg(i).mean_test(l,c,k);
                end
                
                PI_sd(l,c,k)=nanstd(seriestofilter);
                [X, T]=detrend_NAN(seriestofilter); % X is the increments to be added or removed to form a detrended series
                % T is the actual trend line
                PI_30mean(l,c,k)=nanmean(seriestofilter);
                PI_detrended(l,c,k)= {X + PI_30mean(l,c,k)} ; %this is the 30yr mean to be added back to the detrended series to get corrected for time trend
                
                % figure
                % plot(seriestofilter)
                % hold on
                % plot(PI_detrended{l,c,k})
                % hold on
                % plot(T)
                % legend ('Original','Detrended data','Trendline')
                % 
                
            end
        end
    end
    
    
    PI_sd(PI_sd==0)=NaN;
    
    fprintf('Now Saving Detrended data and 30yr mean and 30 yr sd \n')
    save (sprintf('/media/nun008/Seagate/Leo_Files/PI_NorESM/Pre_Industrial_Means/30yrMeans/PI_%s_30yr_detrended.mat',variablelist{varIDX}), 'PI_detrended', '-v7.3')
    save (sprintf('/media/nun008/Seagate/Leo_Files/PI_NorESM/Pre_Industrial_Means/30yrMeans/PI_%s_30yr_mean.mat',variablelist{varIDX}), 'PI_30mean', '-v7.3')
    save (sprintf('/media/nun008/Seagate/Leo_Files/PI_NorESM/Pre_Industrial_Means/30yrMeans/PI_%s_30yr_sd.mat',variablelist{varIDX}), 'PI_sd', '-v7.3')
    clear partial_avg
    
end

