% Overlay the ROIs on a z-stack image of the movie (I)

%load the wspace file from the CNMFe_mdata_readmlist first;
figure; 
I ='AVG_2020-07-10-15-51-57_video-pp-BP-MC(2)_crop.tif';
imshow(I); hold on;
A = 0;
for i = 1:size(neuron.Coor,1)
    plot(neuron.Coor{i}(1,:), neuron.Coor{i}(2,:), 'k-');
    hold on;
    A = A+1;
    xCenter = mean(neuron.Coor{i}(1,:));
    yCenter = mean(neuron.Coor{i}(2,:));
    text(xCenter, yCenter, num2str(A), 'Color', 'y');
    %keyboard;
end

