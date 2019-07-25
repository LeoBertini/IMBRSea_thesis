
%%LINUX GSW LOADING%
addpath(genpath('/Home/siv30/nun008/Desktop/CPLEX/cplex/matlab/x86-64_linux'))
addpath(genpath('/Home/siv30/nun008/Desktop/CPLEX/cplex/examples/src/matlab'))
addpath(genpath('/Home/siv30/nun008/Desktop/GSW'))

% 
% %MAC GSW LOADING %
% addpath(genpath('/Users/leonardobertini/Desktop/GSW'))
% addpath(genpath('/Applications/CPLEX_Studio129/cplex/matlab/x86-64_osx'))
% addpath(genpath('/Applications/CPLEX_Studio129/cplex/examples/src/matlab'))

varlist={'templvl';'salnlvl';'o2'};
varname='AOU';
varIDX=1;

%Need to grab pre-ind physics files and concatenate with ramp and ext series

folderpath=('/media/nun008/Seagate/Leo_Files/PI_NorESM/Pre_Industrial_Means/');

% % listings1=dir(sprintf('/media/nun008/Seagate/Leo_Files/PI_NorESM/Pre_Industrial_Means/mean_PI_%s_*.mat',varlist{varIDX}));
% % 
% % for i=1:length(listings1)
% % dummy=load([folderpath listings1(i).name]);
% % dummy=dummy.mean_test;
% % % bit1=dummy(200:end,:,:);
% % % bit2=dummy(1:199,:,:);
% % % dummy=[bit1 ; bit2];
% % % dummy=permute(dummy,[2,1,3]);
% % piDATA_Temperature(i,1)={dummy}; %this is where all the Temperature pi_data is...
% % end
% % 
% % %now build salinity 
% % 
% % listings2=dir(sprintf('/media/nun008/Seagate/Leo_Files/PI_NorESM/Pre_Industrial_Means/mean_PI_%s_*.mat',varlist{varIDX+1}));
% % 
% % for i=1:length(listings2)
% % dummy=load([folderpath listings2(i).name]);
% % dummy=dummy.mean_test;
% % % bit1=dummy(200:end,:,:);
% % % bit2=dummy(1:199,:,:);
% % % dummy=[bit1 ; bit2];
% % % dummy=permute(dummy,[2,1,3]);
% % piDATA_Salinity(i,1)={dummy}; %this is where all the Salinity pi_data is...
% % end
% % 
% % 
% % %now build oxygen 
% % 
% % listings3=dir(sprintf('/media/nun008/Seagate/Leo_Files/PI_NorESM/Pre_Industrial_Means/mean_PI_%s_*.mat',varlist{varIDX+2}));
% % 
% % for i=1:length(listings3)
% % dummy=load([folderpath listings3(i).name]);
% % dummy=dummy.mean_test;
% % % bit1=dummy(200:end,:,:);
% % % bit2=dummy(1:199,:,:);
% % % dummy=[bit1 ; bit2];
% % % dummy=permute(dummy,[2,1,3]);
% % piDATA_Oxygen(i,1)={dummy}; %this is where all the Oxygen pi_data is...
% % end

%now that we have imported T S O2 

%% loading  depths

depth=load('/media/nun008/Seagate/Leo_Files/PI_NorESM/Pre_Industrial_Means/Depth_Levels.mat'); %get depth array from one of the files in the folder
depth=depth.depths;

%%
 plon=ncread('/media/nun008/Seagate/Leo_Files/Scripts/grid_gx1v6.nc','plon') ;
 plat=ncread('/media/nun008/Seagate/Leo_Files/Scripts/grid_gx1v6.nc','plat') ;
 

   %% calculate AOU for each PI_year and save --> done for whole globe PI
%    
%    for i=1:length(listings1)
%    tic 
%     sprintf('Calculating PI AOU annual means - Handling ### year %d ### data \n',i)
%        for l=1:320
%            for c=1:384
%                for k=1:length(depth)
%                        
%                    Temp=piDATA_Temperature{i,1}(l,c,k);
%                    Sal=piDATA_Salinity{i,1}(l,c,k);
%                    O2in=piDATA_Oxygen{i,1}(l,c,k); %units [mol m-3] needs conversion into [umol kg-1]
%                    
%                    if ~isnan(Temp) && ~isnan(Sal) && ~isnan(O2in)
%                        
%                        z=gsw_z_from_depth(depth(k));
%                        pressure=gsw_p_from_z(z,plat(l,c));
%                        
%                        Ptemp=gsw_pt0_from_t(Sal,Temp,pressure);
%                        Ctemp=gsw_CT_from_pt(Sal,Ptemp);
%                        
%                        O2sat=gsw_O2sol(Sal,Ctemp,pressure,plon(l,c),plat(l,c)); %units O2sol = solubility of oxygen [umol/kg ]
%                        
%                        
%                        density=gsw_rho(Sal,Ctemp,pressure); %[kg m-3]
%                        O2in_new=(O2in/density)*(10^6); %[umol/kg]
%                        
%                        mean_test(l,c,k)=O2sat-O2in_new; %[umol/kg] this is the annual mean AOU field of the PI
%                        
%                    else
%                        mean_test(l,c,k)=NaN;
%                    end
%                    
%                end
%            end
%        end
% 
%        tag=listings1(i).name(17:25);      
%        fprintf('saving individual AOU series for year =%1d \n',i);
%        save(sprintf('/media/nun008/Seagate/Leo_Files/PI_NorESM/Pre_Industrial_Means/mean_PI_%s_%s.mat',varname,tag),'mean_test')
%    toc
%    end

%    clear piDATA_* dummy

   %% Now onto the ramp files

   plon=plon';
   plat=plat';
   
   folderpath='/media/nun008/Seagate/Leo_Files/Ramps/Joint_Physics/';
   folderpathoxy='/media/nun008/Seagate/Leo_Files/Ramps/Joint_BGC/';
   
   for k=8:length(depth)
       tic
       listings1=dir(sprintf('/media/nun008/Seagate/Leo_Files/Ramps/Joint_Physics/new_series_%s_k_%d.mat',varlist{varIDX},k));
       dummy_T=load([folderpath listings1.name]);
       dummy_T=dummy_T.new_series;
       
       listings2=dir(sprintf('/media/nun008/Seagate/Leo_Files/Ramps/Joint_Physics/new_series_%s_k_%d.mat',varlist{varIDX+1},k));
       dummy_Sal=load([folderpath listings2.name]);
       dummy_Sal=dummy_Sal.new_series;
       
       listings3=dir(sprintf('/media/nun008/Seagate/Leo_Files/Ramps/Joint_BGC/new_series_%s_k_%d.mat',varlist{varIDX+2},k));
       dummy_Oxy=load([folderpathoxy listings3.name]);
       dummy_Oxy=dummy_Oxy.new_series;
       

       for l=1:384
           
           sprintf('depth level = %d | line = %d', k,l)
           
           for c=1:320
               
               %if there are 480 points in the time series, then calculate
               %AOU
               if numel(~isnan(dummy_Oxy{l,c}))==481 && numel(~isnan(dummy_T{l,c}))==481 && numel(~isnan(dummy_Sal{l,c}))==481
                   
                   for i=1:length(dummy_Oxy{l,c})
                       
                       Temp=dummy_T{l,c}(i);
                       Sal=dummy_Sal{l,c}(i);
                       O2in=dummy_Oxy{l,c}(i);
                       
                       z=gsw_z_from_depth(depth(k));
                       pressure=gsw_p_from_z(z,plat(l,c));
                       
                       Ptemp=gsw_pt0_from_t(Sal,Temp,pressure);
                       Ctemp=gsw_CT_from_pt(Sal,Ptemp);
                       
                       O2sat=gsw_O2sol(Sal,Ctemp,pressure,plon(l,c),plat(l,c)); %units O2sol = solubility of oxygen [umol/kg ]
                       
                       
                       density=gsw_rho(Sal,Ctemp,pressure); %[kg m-3]
                       O2in_new=(O2in/density)*(10^6); %[umol/kg]
                       dummy_AOU(i)=O2sat-O2in_new; %[umol/kg]
                   end
                   new_series(l,c)={dummy_AOU};
               else
                   new_series(l,c)={NaN};
               end
               
           end
           
       end
       
       save(sprintf('/media/nun008/Seagate/Leo_Files/Ramps/Joint_BGC/new_series_%s_k_%d.mat',varname,k),'new_series')
       toc
   end
  
   
   
   
   
   
   
   % plot(ramp_AOU{l,c})
   
%    l=150;c=150;
%    figure
%    subplot(1,3,1)
%    plot(1:481,dummy_Oxy{l,c})
%   
%   subplot(1,3,2)
%   plot(1:481,dummy_Sal{l,c})
%    
%    subplot(1,3,3)
%    plot(1:481,dummy_T{l,c})

    
   
   
