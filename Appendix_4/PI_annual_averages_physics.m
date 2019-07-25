clear all 
close all

%THIS IS NOW SET TO CALCULATE PARTIAL MEANS OF TEMPERATURE AND SALINITY
%since their PI files
%path=['/shared/projects/uniklima/globclim/nun008/NorESM/PI/'];

%##TO USE LINUX MACHINE WITH 4TB HARD DRIVE#
cd ('/media/nun008/Seagate/Leo_Files/PI_NorESM/')
listings=dir(['/media/nun008/Seagate/Leo_Files/PI_NorESM/Physics/' 'N*.nc']);
load('/media/nun008/Seagate/Leo_Files/PI_NorESM//Physics/Depth_Levels.mat')

% %##TO USE WITH PERSONAL HARD DRIVE#
% cd('/Volumes/LaCie_Leonardo/NorESM/Pre_industrial')
% listings=dir(['/Volumes/LaCie_Leonardo/NorESM/Pre_industrial/' 'N*.nc']);
% depth=load('/Volumes/LaCie_Leonardo/NorESM/all_ramps/Depth_Levels.mat'); 
% 
%% variables

varlist={'templvl';'salnlvl'};


filenames_org = {listings.name};
filenames_org= filenames_org';

 %general loop should start here to go through all the pre-industrial files
jump=11;    
indx_start=1; %this is where we indicate the first file that's going to be read 
indx_end=indx_start+jump;
times=1;

for varIDX=1%:length(varlist)

while times<=30 %the amount of years

    
 for a=indx_start:indx_end %--> here the loop for all the .nc files starts with a and ends with b
     
    FN=['/media/nun008/Seagate/Leo_Files/PI_NorESM/Physics/' char(filenames_org(a))]; %converts into characther array with complete path of target nc file to be open
    %FN=['/shared/projects/uniklima/globclim/nun008/NorESM/PI/' char(filenames_org(a))]; %converts into characther array with complete path of target nc file to be open
    %ncdisp(FN);
    
    %show file number
    fprintf('File number=%.f\n',a)
    fprintf('\n')
    fprintf(FN) 
    %importing variables from nc file
      
%     thickness=ncread(FN, 'pddpo');
   
    varField=ncread(FN,varlist{varIDX}); %units = degC
   


%this is a structure array where all the matrices will be stored with the corresponding filenames.
VAR_DATA(a).File_Name=listings(a).name;   
VAR_DATA(a).Field=varField;


 end
 
 VAR_DATA_12months=VAR_DATA((indx_start:indx_end));
 
 
 fprintf('All monthly fields imported and saved in structure\n')
 fprintf('\n')
 
%name for the files
name1=filenames_org(indx_start);
name1=name1{1,1}(39:42);
name2=filenames_org(indx_end);
name2=name2{1,1}(39:42);



%computing the means 
fprintf('Calculating means\n')  

mean_test=zeros(320,384,70); %allocating space
mean_test2=zeros(320,384,70); %allocating space


tic
 for i=1:320
     i
        for j=1:384  
        
                    for k=1:length(depths) %for each depth

                      for h=1:length(VAR_DATA_12months) %for each file in the structure array 
                     
                          %Temperature
                          test(h,k)= VAR_DATA_12months(h).Field(i,j,k); 
                          %test((h-indx_start+1),k)= VAR_DATA(h).Field(i,j,k); 
                      end
                      
                      mean_test(i,j,k)=nanmean(test(:,k)); 
                    
                      
                    end

            end
        
               
               
 end
           

toc

elapsedTime = toc;
fprintf('Elapsed time = %.2f min \n',elapsedTime/60);


name_file2=sprintf('/media/nun008/Seagate/Leo_Files/PI_NorESM/Pre_Industrial_Means/mean_PI_%s_%s_%s.mat',varlist{varIDX},name1, name2);
save (name_file2,'mean_test', '-v7.3')


indx_start=indx_end+1; %this is where we indicate the first file that's going to be read in the next loop
indx_end=indx_start+jump;
times=times+1;

clear VAR_DATA_12months

end

clear VAR_DATA
indx_start=indx_end+1; %this is where we indicate the first file that's going to be read in the next loop
indx_end=indx_start+jump;
times=times+1;

end

