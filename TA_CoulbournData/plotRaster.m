function [] = plotRaster(lickMat , tVec)
% Visualize raster plot
%lickMat = logical matrix of when licks happened
%tVec is a vector with the times, should be equal length to lickMat
hold all ;
%spikeMat = randi([0 1], [3 30]); tVec = 0:1:29;
for trialCount = 1: size(lickMat ,1)
    small_spikeMat = logical(lickMat(trialCount, :));
    spikePos = tVec(small_spikeMat);
    for spikeCount = 1: length(spikePos)
        plot([spikePos(spikeCount) spikePos(spikeCount)] , [trialCount-0.4
            trialCount+0.4], 'k') ;
    end
end
end
