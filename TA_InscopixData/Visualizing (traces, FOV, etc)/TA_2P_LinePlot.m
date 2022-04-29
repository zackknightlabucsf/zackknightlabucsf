%% Input Stuffies

ToPlot = []; %put in your array here (rows=neurons, columns = time)
factor = 5; %by how many points do you want to offset each line? 

%% Plot
y_offset = zeros(1,size(all_mrerank,2))-factor; 

x = 1:size(ToPlot,2); 
figure; 
plot(x, ToPlot(1,:), 'k');
hold on;
a = 1;
for i = 2:size(ToPlot,1)
    plot(x, ToPlot(i,:)+y_offset*a, 'k');
    hold on; 
    a = a+1;
end

% Old code below if you want to plot licks over it
% hold on; %plot licks over it
% y1 = ylim;
% 
ttl_logicNew = [];
for i = 1:length(ttl_clean)
    if ttl_clean(i) == 0
        ttl_logicNew(i) = -10;
    else
        ttl_logicNew(i) = 90;
    end
end
% plot(ttl_logicNew);
% 
% title('TA124 fasted Ensure Spaced 1014');
