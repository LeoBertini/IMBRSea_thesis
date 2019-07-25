%% CODE

depth=load('/Volumes/LaCie_Leonardo/NorESM/all_ramps/Depth_Levels.mat'); %get depth array from one of the files in the folder
depth=depth.depths;
depth=depth';

depth_interp=min(depth):200:max(depth);

parea=ncread('/Volumes/LaCie_Leonardo/NorESM/all_ramps/grid_gx1v6.nc','parea');
plon=ncread('/Volumes/LaCie_Leonardo/NorESM/all_ramps/grid_gx1v6.nc','plon') ;
plat=ncread('/Volumes/LaCie_Leonardo/NorESM/all_ramps/grid_gx1v6.nc','plat') ;

plat=plat';
plon=plon';
plon=plon;

addpath(genpath('/Volumes/LaCie_Leonardo/NorESM/'))


atl_ind=load('/Volumes/LaCie_Leonardo/NorESM/scripts_jerry/noresm_atl_ind_fixed.asc');

depthbds=ncread('/Volumes/LaCie_Leonardo/NorESM/Pre_industrial/N1850AERCNOC_f19_g16_CTRL_02.micom.hbgcm.1000-12.nc','depth_bnds');
for k=1:length(depthbds)
    thickness(k)=depthbds(2,k)-depthbds(1,k);
end


varlist={'templvl';'ph';'o2';'AOU';'omegac';'detoc'};

fig=figure('visible','on')
set(fig,'Units', 'pixels', 'Position', [0, 0, 764, 1000]);

for varIDX=1:5
    
    
    if strcmp(varlist{varIDX},'ph')==1
        varname='pH';
        axmaxred=100;
        axmaxblue=1000;
        positionsubplot=2;
        letter='b';
        
    elseif strcmp(varlist{varIDX},'o2')==1
        varname='DO';
        axmaxred=200;
        axmaxblue=2200;
        positionsubplot=3;
        letter='c';
        
    elseif strcmp(varlist{varIDX},'omegac')==1
        varname='\Omega_C';
        axmaxred=100;
        axmaxblue=1500;
        positionsubplot=5;
        letter='e';
        
    elseif strcmp(varlist{varIDX},'detoc')==1
        varname='POC flux';
        axmaxred=100;
        axmaxblue=2000;
        positionsubplot=6;
        letter='f';
        
    elseif strcmp(varlist{varIDX},'templvl')==1
        varname='Temperature';
        axmaxred=100;
        axmaxblue=2000;
        positionsubplot=1;
        letter='a';
        
    elseif strcmp(varlist{varIDX},'AOU')==1
        varname='AOU';
        axmaxred=200;
        axmaxblue=1500;
        positionsubplot=4;
        letter='d';
    end
    
    
    fields=load(sprintf('/Volumes/LaCie_Leonardo/NorESM/all_ramps/filtered/Results/New_time_var_results/%s_results_updated_v1D_poly1_NEW_NEW.mat',varlist{varIDX}));
    % fieldZERO=fields.TOD;
    % fieldZERO=fields.PeakVALUE;
    % fieldZERO=fields.Trecovery;
    %
    % fieldZERO(fieldZERO==-1111)=NaN; %this makes the points representing the bottom 'transparent'
    % %changing all the negative flags to avoid the concentric rings when interpolating
    % fieldZERO(fieldZERO<0)=600;
    %
    % %unflipping
    % bit1=fieldZERO(:,122:end,:); %bit1 (200:end,:,:); %reshaping -- position 200 is where Pacific is
    % bit2=fieldZERO(:,1:121,:);   %bit2 (1:199,:,:);
    % field1=cat(2,bit1,bit2);
    %
    % k=[find(depth==500) find(depth==1000) find(depth==2000)];
    
    %select and average all points within the lat lims for that particular
    %depth
    
    %% TOD
    
    for sdindx=1:3
        
        fieldZERO=fields.TOD{1,sdindx};
        %unflipping
        bit1=fieldZERO(:,122:end,:); %bit1 (200:end,:,:); %reshaping -- position 200 is where Pacific is
        bit2=fieldZERO(:,1:121,:);   %bit2 (1:199,:,:);
        field1(1,sdindx)={cat(2,bit1,bit2)};
        
        
        for i=1:70     %calculating the mean for each depth layer for a set range of latitudes
            
            tgtlayer=field1{1,sdindx}(:,:,i);
            tgtlayer=tgtlayer(:);
            indxlat=find(plat<=65 & plat>=0);
            
            
            %         numerator1(i,:)=(tgtlayer(indxlat).*parea(indxlat))*thickness(i);
            %         denominator1(i,:)=parea(indxlat)*thickness(i);
            %
            %         dummy1(i)=nansum(numerator1(i,:))/nansum(denominator1(i,:));
            
            TOD_vertprof(sdindx,i)={nanmean(tgtlayer(~isnan(tgtlayer(indxlat))))};
            TOD_sderror(sdindx,i)={nanstd(tgtlayer(~isnan(tgtlayer(indxlat))))/sqrt(length(~isnan(tgtlayer(indxlat))))};
            
        end
        
    end
    
    
    
    %% Trec
    for sdindx=1:3
        
        
        fieldZERO=fields.Trecovery{1,sdindx};
        fieldZERO(fieldZERO==-1111)=NaN; %this makes the points representing the bottom 'transparent'
        %changing all the negative flags to avoid the concentric rings when interpolating
        fieldZERO(fieldZERO<0)=600;
        
        %unflipping
        bit1=fieldZERO(:,122:end,:); %bit1 (200:end,:,:); %reshaping -- position 200 is where Pacific is
        bit2=fieldZERO(:,1:121,:);   %bit2 (1:199,:,:);
        field1(1,sdindx)={cat(2,bit1,bit2)};
        
        for i=1:70
            tgtlayer=field1{1,sdindx}(:,:,i);
            tgtlayer=tgtlayer(:);
            indxlat=find(plat<=65 & plat>=0);
            
            
            %         numerator2(i,:)=(tgtlayer(indxlat).*parea(indxlat))*thickness(i);
            %         denominator2(i,:)=parea(indxlat)*thickness(i);
            %
            %         dummy2(i)=nansum(numerator2(i,:))/nansum(denominator2(i,:));
            
            Trec_vertprof(sdindx,i)={nanmean(tgtlayer(~isnan(tgtlayer(indxlat))))};
            Trec_sderror(sdindx,i)={nanstd(tgtlayer(~isnan(tgtlayer(indxlat))))/sqrt(length(~isnan(tgtlayer(indxlat))))};
        end
    end
    
    
    %% figure
    subplot(3,2,positionsubplot)
   
    %filled area
    x1=smooth(cell2mat(TOD_vertprof(1,:)));
    x3=smooth(cell2mat(TOD_vertprof(3,:)));
    X = [x1(:).', fliplr(x3(:).')];
    Y = [-depth, fliplr(-depth)];
    shaded1=fill(X, Y, [255/255 204/255 204/255]);
    shaded1.FaceAlpha=0.4;
    shaded1.LineStyle='none';
    
    hold on
%     p1=line(x1,-depth,'linewidth',2,'color','r'); %TOD lower boundary
%     hold on
%     p3=line(x3,-depth,'linewidth',2,'color','r'); %TOD upper boundary
    hold on
    x2=smooth(cell2mat(TOD_vertprof(2,:)));
    p2=line(x2,-depth,'linewidth',2,'color','r','LineStyle','-'); %TOD mean
    
    grid on
    ax1 = gca; % current axes
    ax1.XColor = 'r';
    ax1.YColor = 'r';
    ax1.LineWidth=2;
    ax1.XLim=[0 axmaxred];
    ax1.YLim=[-4000 0];
    ax1.GridColor=[.5 .5 .5];
    ax1.YMinorTick='on';
    ax1.MinorGridColor=[.5 .5 .5];
    ax1.GridAlpha=0.4;
    ax1.FontName='Helvetica Neue';
    ax1.FontSize=16;
    ax1.YLabel.String='Depth (m)';
    ax1.XLabel.String={'ToD (yr)'};
    ax1.XTick=linspace(0,axmaxred,5);
    grid minor
        
    ax1_pos = ax1.Position; % position of first axes
    ax2 = axes('Position',ax1_pos,...
        'XAxisLocation','top',...
        'YAxisLocation','right',...
        'Color','none');
    
    
    hold on
    
    x10=smooth(cell2mat(Trec_vertprof(1,:)));
    x30=smooth(cell2mat(Trec_vertprof(3,:)));

    %filled area
    X = [x10(:).', fliplr(x30(:).')];
    Y = [-depth, fliplr(-depth)];
    shaded2=fill(X, Y, [102/255 178/255 225/255])
    shaded2.FaceAlpha=0.4;
    shaded2.LineStyle='none';

%     hold on
%     p10=line(x10,-depth,'linewidth',2,'color','b'); %TOD lower boundary
%     hold on
%     p30=line(x30,-depth,'linewidth',2,'color','b'); %TOD lower boundary
    hold on
    x20=smooth(cell2mat(Trec_vertprof(2,:)));
    p20=line(x20,-depth,'linewidth',2,'color','b','LineStyle','-'); %TREC
   
    ax2.FontName='Helvetica Neue';
    ax2.FontSize=16;
    ax2.LineWidth=2;
    
    ax2.YLabel.String='Depth (m)';
    ax2.XLabel.String='Trec (yr)';
    ax2.XColor = 'b';
    ax2.YColor = 'b';
    ax2.XLim=[0 axmaxblue];
    ax2.YLim=[-4000 0];
    ax2.GridColor=[.5 .5 .5];
    ax2.YMinorTick='on';
    ax2.MinorGridColor=[.5 .5 .5];
    ax2.GridAlpha=0.4;
    ax2.XTick=linspace(0,axmaxblue,5);
    grid minor
    grid on

    lgd=legend([p2 p20],'ToD','Trec')
    lgd.Location='southeast';
    lgd.Color=[.8 .8 .8];
    title(sprintf('%s) %s',letter,varname))
    uistack(lgd,'top')
    
    
    if positionsubplot==1 || positionsubplot==3 %|| positionsubplot==3 || positionsubplot==4
        ax2.YTick=[];
        ax2.YLabel=[];
    end
    
    if positionsubplot==2 || positionsubplot==4 %|| positionsubplot==3 || positionsubplot==4
        ax1.YTick=[];
        ax1.YLabel=[];
        
    end
    
    if positionsubplot==3 || positionsubplot==4 || positionsubplot==5 %|| positionsubplot==4
        ax2.XLabel=[];
    end
    
    if positionsubplot==1 || positionsubplot==2 || positionsubplot==3 %|| positionsubplot==4
        ax1.XLabel=[];
    end
    
    
end
set(gcf,'color','w')
fig.InvertHardcopy = 'off';
supersizeme('-')
print('/Volumes/LaCie_Leonardo/NorESM/Initial_figs/vertical_profile_time_vars.png','-dpng','-r300')

% H=annotation('textbox','String', sprintf('a) %s',varname));
% H.BackgroundColor='w';
% H.Position=[lgd.Position(1)-0.6 lgd.Position(2) lgd.Position(3)+0.1 lgd.Position(4)+0.02];
% H.HorizontalAlignment='center';
% H.FontSize=16;
% H.VerticalAlignment='middle';

