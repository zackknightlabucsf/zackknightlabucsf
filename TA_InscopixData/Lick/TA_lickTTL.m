%% Use this code to pull out the TTLs to plot into prism    
[ttl] = csv_reshape_for_inscopix_analysis_JG('2021-10-29-16-17-08_video_green-gpio.csv','GPIO-1', 4, 9600);
temp=ttl;
%check for NaNs
positionNaN = isnan(temp);
temp(positionNaN) = 0;
ttl_logic = (ttl>2000);

%making it the same length as the calcium data traces
    ttl_new= decimate(temp,4); %4 is fr
    ttl_logicNew = (ttl_new>2000);
    
%optional plot to see if they look the same
%     figure; 
%     subplot(1,2,1);
%     plot(1:9600, ttl_logic, 'k'); hold on; ylim([-0.5, 1.5]);
%     subplot(1,2,2); 
%     plot(1:2400, ttl_logicNew, 'k'); hold on; ylim([-0.5, 1.5]);

%rescaling it; turns 1 to 20 and 0 to -5
    ttl_logicNew = ttl_logicNew*20;
    for i = 1:length(ttl_logicNew)
        if ttl_logicNew(i) == 0
            ttl_logicNew(i) = -5;
        end
    end
    
%optional plot to see if they look the same
%     figure; 
%     subplot(1,2,1);
%     plot(1:9600, ttl_logic, 'k'); hold on; ylim([-0.5, 1.5]);
%     subplot(1,2,2); 
%     plot(1:2400, ttl_logicNew, 'k'); hold on; ylim([-10, 25]);    