close all
clear all
fclose all
total_length = 600 + 3300;
list={
'TXT_FILE'

    };
[colormap]=cbrewer('div', 'PRGn', 256); %PRGn
colormap=flip(colormap);
%% basics
for i=1:length(list)
    txtfile=strcat(list{i},'.txt'); %combines all the neurons/columns from each mouse in the file listed in the text document
%     
    zscoreyn=1; % 1 is yes
    rawyn=1; %normally set to 1
%     [ cC, mpsthnom,c_psthnom,all_mrerank,ind_mrerank,crerank,ind_crerank,c_corrm,para,preposttimerange,cC_rerank] = cnmfe_yc_mdata_readmtxt( txtfile,zscoreyn,rawyn);
%     [ cC, mpsthnom,c_psthnom,all_mrerank,ind_mrerank,crerank,ind_crerank,c_corrm,para,preposttimerange,cC_rerank] = cnmfe_yc_mdatalowpassfilter_readmtxt( txtfile,zscoreyn,rawyn);
   [ cC, mpsthnom,c_psthnom,all_mrerank,ind_mrerank,crerank,ind_crerank,c_corrm,para,preposttimerange,...
       cC_rerank] = BCJ_cnmfe_ta_mdatalowpassfilter_readmtxt_decimateto1 (txtfile,zscoreyn,rawyn);

    for i2=1:size(all_mrerank,1)
        all_mrerank(i2,:)=all_mrerank(i2,:)-mean(all_mrerank(i2,1:300));
    end
    
    myData = all_mrerank;
    
    if zscoreyn
        % Decide how wide you want the white band:
%         % e.g. any value whose absolute value < 1 becomes 0 (white).
%         threshold = 1;
%         myData(abs(myData) < threshold) = 0;
        
        h=heatmap_d(all_mrerank,[],[],[],'MaxColorValue',4,'MinColorValue',-4,'Colormap',colormap,'Colorbar',1);
%         h=heatmap_d(myData,[],[],[],'MaxColorValue',5,'MinColorValue',-5,'Colormap',colormap,'Colorbar',1);
        hold on;
        Title = list{1};
        title(Title);
        y = 1:400; z = zeros(1,400)+600; %
%         x = zeros(1,400)+1500;
%         m = zeros(1,400)+2700;%this is to make lines at desired time
        plot(z,y,'k', 'LineWidth', 2); %plot(x,y,'k', 'LineWidth', 2)
        hold on;
%         plot(x,y,'k', 'LineWidth', 2);
        hold on;
%         plot(m,y,'k', 'LineWidth', 2); hold on;
        yticks([1]); yticklabels([size(all_mrerank,1)]);
        
        figure; 
        for i2 =1:size(all_mrerank,1)
        plot(1:total_length,all_mrerank(i2,:)); hold on; 
        end

%         figure; 
%         offset = 0;
%         for i2 =1:size(Neurons_actA,1)
%         plot(1:3600,Neurons_actA(i2,:) + offset); 
%          offset = offset + 50;  
%          hold on;
%         end
        

        saveas(h,strcat(list{i},'_zs_hmaplong.eps'),'epsc');
        savefig(strcat(list{i},'_zs_hmaplong.fig')); %save the figure as .fig
        save(strcat(list{i},'_zs_wspacelong.mat'))
        
       %[Neurons_act, Neurons_none,Neurons_inhib]=Heatmap_quant_4Brooke(strcat(list{i},'_dff_wspace'), Title, 599, 1800);
       
        save(strcat(list{i},'_dff_wspace.mat'))
  else
        h=heatmap_d(all_mrerank,[],[],[],'MaxColorValue',4,'MinColorValue',-4,'Colormap',colormap,'Colorbar',1);
    
    saveas(h,strcat(list{i},'_dff_hmap_licks.eps'),'epsc');
    savefig(strcat(list{i}, '_dff_hmap_licks.fig')); % Save the figure as .fig
    save(strcat(list{i},'_dff_wspace.mat'))
    end


end