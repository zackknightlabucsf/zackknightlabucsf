%% Perform unsupervised kmeans clustering on single-cell data
myMatrix = []; %Put here whichever matrix you want to classify. I usually use all_mrerank
klist=2:10;%the number of clusters you want to try
eva = evalclusters(myMatrix,'kmeans','Silhouette','klist',klist);


IDX=kmeans(myMatrix,eva.OptimalK, 'Replicates', 100);

%make a figure and plot kmeans in separate subplots
figure; 
for i = 1:eva.OptimalK
   ToPlot = myMatrix(find(IDX==i),:);
   subplot(eva.OptimalK,1,i);
   h=heatmap_d(ToPlot,[],[],[],'MaxColorValue',1,'MinColorValue',-1,'Colormap',colormap,'Colorbar',1);
   hold on; 
   %Uncomment below if you want to add lines to a graph
%    for j = 1:length(a) %lines to mark where taste was
%     plot(zeros(1,50)+a(j), 1:50, 'k', 'LineWidth', 2);
%    end
end

