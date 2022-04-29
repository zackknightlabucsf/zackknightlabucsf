close all
%clear all
clearvars -except these_neurons
fclose all

list={
% Put in a text file where Row 1 Column 1 is the title, Row 2 column 1 is
% small file name, press tab, put in fr, press tab, stim time, press tab,
% second before, press tab, seconds after. 

%'BNSTtoDMH_FastedSucralose (ensure)'
    };
[colormap]=cbrewer('div', 'RdBu', 256);
colormap=flip(colormap);
%% basics
for i=1:length(list)
    txtfile=strcat(list{i},'.txt');
%     
    zscoreyn=1; % 1 is yes
    rawyn=1; %normally set to 1
%     [ cC, mpsthnom,c_psthnom,all_mrerank,ind_mrerank,crerank,ind_crerank,c_corrm,para,preposttimerange,cC_rerank] = cnmfe_yc_mdata_readmtxt( txtfile,zscoreyn,rawyn);
%     [ cC, mpsthnom,c_psthnom,all_mrerank,ind_mrerank,crerank,ind_crerank,c_corrm,para,preposttimerange,cC_rerank] = cnmfe_yc_mdatalowpassfilter_readmtxt( txtfile,zscoreyn,rawyn);
   [ cC, mpsthnom,c_psthnom,all_mrerank,ind_mrerank,crerank,ind_crerank,c_corrm,para,preposttimerange,...
       cC_rerank] = cnmfe_ta_mdatalowpassfilter_readmtxt_decimateto1(txtfile,zscoreyn,rawyn);

    for i2=1:size(all_mrerank,1)
        all_mrerank(i2,:)=all_mrerank(i2,:)-mean(all_mrerank(i2,1:600));
    end
    if zscoreyn
        h=heatmap_d(all_mrerank,[],[],[],'MaxColorValue',4,'MinColorValue',-4,'Colormap',colormap,'Colorbar',1);
        hold on;
        Title = list{1};
        title(Title);
        y = 1:400; 
        z = zeros(1,400)+600; 
        x = zeros(1,400)+1800;
        m = zeros(1,400)+1800;%this is to make lines at desired time
        plot(z,y,'k', 'LineWidth', 2);
        hold on;
        %plot(x,y,'k', 'LineWidth', 2); %uncomment these to plot the lines
        %hold on;
        %plot(m,y,'k', 'LineWidth', 2); hold on;
        yticks([1]); yticklabels([size(all_mrerank,1)]);


        saveas(h,strcat(list{i},'_dff_hmap.eps'),'epsc');
        save(strcat(list{i},'_dff_wspace.mat'))
        
    else
        h=heatmap_d(all_mrerank,[],[],[],'MaxColorValue',10,'MinColorValue',-10,'Colormap',colormap,'Colorbar',1);
    
    saveas(h,strcat(list{i},'_dff_hmap.eps'),'epsc');
         save(strcat(list{i},'_dff_wspace.mat'))
    end

end

