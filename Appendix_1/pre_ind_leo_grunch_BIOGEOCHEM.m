clear all 
close all

   %variables in biogechem ncfiles are ...
    ... o2=ncread(FN, 'o2'); %units = 'mol O2 m-3' --> to convert to ?mol/kg --> (o2/sigma)*E6 
    ... thickness is 'pddpo', 
    ... sigma= ncread(FN, 'sigma'); %units = 'kg m-3'
    ... omgegac (calcite saturation state)
    ... detoc=ncread(FN,'detoc'); %units = 'mol C m-3' --> detritus

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
varlist={'ph';'o2';'omegac';'detoc'}; %list of variables to have their data interpolated and then monthly data reduced to annual

%%
filenames_org = {listings.name};
filenames_org= filenames_org';

 

jump=11; %this is the 'jump window of 'jump+1'--> reads a sequence of 12 files
indx_start=1; %this is where we indicate the first file that's going to be read 
indx_end=indx_start+jump;
times=1;

% jump=11; %this is the 'jump window of 'jump+1'--> reads a sequence of 12 files
% indx_start=349; %this is where we indicate the first file that's going to be read 
% indx_end=indx_start+jump;
% times=29;

for varIDX=4:length(varlist)
    
    
    while times<=30 %the amount of years
        
        
        for a=indx_start:indx_end %--> here the loop for all the .nc files starts with a and ends with b
            
            FN=['/Volumes/LaCie_Leonardo/NorESM/Pre_industrial/' char(filenames_org(a))]; %converts into characther array with complete path of target nc file to be open
            
            %FN=['/shared/projects/uniklima/globclim/nun008/NorESM/PI/' char(filenames_org(a))]; %converts into characther array with complete path of target nc file to be open
            %ncdisp(FN);
            
            %show file number
            fprintf('File number=%.f\n',a)
            
            %importing variables from nc file
            
            VarField=ncread(FN, varlist{varIDX});
            
            
            thickness=ncread(FN, 'pddpo');
            
            %now that we have the depths for each pH point, we ended up with 53
            % depth layers, which have a pH value associated.
            %However, the ramp-up and ramp-down dataset has 70 depth levels. So we have
            %to'calculate' pH for each of the 70 depth levels by interpolation in
            %order to, then, calculate the long term mean at each particular depth level.
            
            depth_bnds= ncread(FN, 'depth_bnds')'; %the depth boundaries for the grid cells used in the ramp up and down files
            
            %%%%CALCULATING THE ISOPYCNAL BOUNDARIES
            
            %%initiating varibales to allocate space
            isopycnal_bnds_upper=zeros(size(VarField,1), size(VarField,2), size(VarField,3));
            isopycnal_bnds_lower=zeros(size(VarField,1), size(VarField,2), size(VarField,3));
            
            A = isnan(VarField(:,:,:));  %upper boundary row 1 --> positions with actual values are filled with zero
            A=double(A);
            A(A==1)=NaN;
            isopycnal_bnds_upper(:,:,1) = A(:,:,1);  % 0 for surface
            isopycnal_bnds_lower(:,:,1) = 0 + thickness(:,:,1); %first lower boundary based on thickness
            
            %for isopycnal 2 onwards
            for k=2:size(A,3)
                isopycnal_bnds_upper(:,:,k) = isopycnal_bnds_lower(:,:,k-1);
                isopycnal_bnds_lower(:,:,k)= isopycnal_bnds_upper(:,:,k) + thickness(:,:,k);
            end
            
            %Example of depth of isopycnal boundaries for grid cell (1,3) across all levels
            %boundaries(:,1)= isopycnal_bnds_upper(:,3,1) ;
            %boundaries(:,2)= isopycnal_bnds_lower(:,3,1) ;
            
            
            %%%%INTERPOLATION
            
            fprintf('Interpolation begins\n')
            fprintf(listings(a).name)
            fprintf('\n')
            
            %allocating space for the interp_Var
            interp_Var=NaN(size(VarField,1), size(VarField,2), size(depth_bnds,1));
            
            for i=1:size(VarField,1) %for each level of latitude ("rows")
                for j=1:size(VarField,2) % for each level of "longitude ("columns")
                    
                    if ~isnan(VarField(i,j)) %test if the very fisrt value at the surface level is different from NaN
                        
                        z_org_field_Var=squeeze(VarField(i,j,:)); %creates an array with all the vertical values for that grid cell location in case it's not landmass
                        
                        isopycnal_bnd(:,1)= isopycnal_bnds_upper(i,j,:);
                        isopycnal_bnd(:,2)= isopycnal_bnds_lower(i,j,:);
                        
                        interp_Var(i,j,:)=z_interp_gen_2(depth_bnds,z_org_field_Var,isopycnal_bnd);
                        
                        
                        
                    end
                end
                
            end
            
            
            
            
            %this is a cell array where all the interpolated matrices will be stored.
            %this is used in the calculation of the means
            interpolated_result{a}=interp_Var;
            
            
            dummy1=interp_Var(200:end,:,:);
            dummy2=interp_Var(1:199,:,:);
            dummy=[dummy1;dummy2];
% % % % 
% % % %             k=60;
% % % %             figure
% % % %             subplot(1,2,1)
% % % %             pcolor(dummy(:,:,k)')
% % % %             shading flat
% % % %             colorbar
% % % %             %caxis([7.5 8])
% % % %             caxis([min(min(dummy(x,y,k))) max(max(dummy(x,y,k)))])
% % % %             hold on
% % % %             x=ones(110,1)*120;
% % % %             y=200:309;y=y';
% % % %             scatter(x,y,'filled','red')
% % % %             title(sprintf('Interpolated field for %s at depth %d m',varlist{varIDX},depths(k)))
% % % % 
% % % %             subplot(1,2,2)
% % % %             for i=1:length(x)
% % % %             scatter(y(i),dummy(x(i),y(i),k))
% % % %             hold on
% % % %             end
% % % %             xlabel('Latitude coordinate along transect')
% % % %             ylabel(sprintf('%s value',varlist{varIDX}))
% % % %             title(sprintf('Values of %s along transect at Depth %d m',varlist{varIDX},depths(k)))

            
           
            %this is a structure array where all the interpolated matrices will be stored with the corresponding filenames.
            final_interpolated.File_Name{a}=listings(a).name;
            final_interpolated.Field{a}=interp_Var;
            
            
            fprintf('Interpolation Complete\n')
            fprintf('\n')
            
            
        end
        
        %creates substructure to save interpolation files
        interpolated_struc=final_interpolated.Field(indx_start:end);
        
        % interpolated_struc.File_Name = final_interpolated.File_Name(indx_start:end);%%not sure why this is relevant
        % interpolated_struc.Field = final_interpolated.Field(indx_start:end);
        
        fprintf('###Saving interpolation file###\n')
        
        name1=filenames_org(indx_start);
        name1=name1{1,1}(42:45);
        name2=filenames_org(indx_end);
        name2=name2{1,1}(42:45);
        
        name_file1=sprintf('/Volumes/LaCie_Leonardo/NorESM/Pre_industrial/interpolated_PI_%s_%s_%s.mat',varlist{varIDX},name1,name2);
%         save(name_file1,'interpolated_struc','-v7.3')
%         
        
        clearvars -except filenames_org depths listings indx_start indx_end a jump times name_file1 name1 name2...
            interpolated_result interpolated_struc final_interpolated varlist varIDX
        
        
        %computing the means
        fprintf('Calculating means\n')
        
        mean_test=zeros(320,384,70); %allocating space
        
        tic
        for i=1:320
            i
            for j=1:384
               
                for k=1:length(depths) %for each depth
                    
                    for h=1:length(interpolated_struc)  %for each file in the cell array
                        
                        %pH
                        test(h,k)=interpolated_struc{h}(i,j,k);
                        mean_test(i,j,k)=nanmean(test(:,k));
                        
                    end
                    
                    
                end
                
            end
        end
        
        
        toc
        
        elapsedTime = toc;
        fprintf('Elapsed time = %.2f min \n',elapsedTime/60);
        
        name_file2=sprintf('/Volumes/LaCie_Leonardo/NorESM/Pre_industrial/Partials/mean_PI_%s_%s_%s.mat',varlist{varIDX},name1, name2);
        save (name_file2,'mean_test', '-v7.3')
        
        indx_start=indx_end+1; %this is where we indicate the first file that's going to be read in the next loop
        indx_end=indx_start+jump;
        times=times+1;
        
        clear interpolated_struc final_interpolated
        
    end
    
    jump=11; %this is the 'jump window of 'jump+1'--> reads a sequence of 12 files
    indx_start=1; %this is where we indicate the first file that's going to be read 
    indx_end=indx_start+jump;
    times=1; %set counter back
    
end

