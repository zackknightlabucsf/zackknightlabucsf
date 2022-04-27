%water control A2 in wspace_small in 1/
%A3 needs to be batched
%en/ghrelin A2/3 in small in 1/inscopix data

close all
clear all
fclose all
infile='JG_111121_DV2ex1-pp-BP-MC_lowPNR_wspace';
load(infile);
display_merge = true; % set to true if you want to inspect every candidate merge
view_neurons = false; % set to true if you want to inspect all neurons after quick merge routine

merge_thr = [0.6, 0.6, 0.1];  % choose thresholds for merging neurons (this will priccccmarily merge neurons redundantly found by multiple patch processes, likely in the patch-overlaps)
% merge neurons based on the correlation computed with {'A', 'S', 'C'}
cnmfe_quick_merge;            % run neuron merges

merge_thr = [0.5, 0.5, 0.1];  % choose thresholds for merging neurons (this will primarily merge neurons redundantly found by multiple patch processes, likely in the patch-overlaps)
% merge neurons based on the correlation computed with {'A', 'S', 'C'}
cnmfe_quick_merge;            % run neuron merges

merge_thr = [0.5, 0.3, 0.1];  % choose thresholds for merging neurons (this will primarily merge neurons redundantly found by multiple patch processes, likely in the patch-overlaps)
% merge neurons based on the correlation computed with {'A', 'S', 'C'}
cnmfe_quick_merge;            % run neuron merges

merge_thr = [0.3, 0.5, 0.1];  % choose thresholds for merging neurons (this will primarily merge neurons redundantly found by multiple patch processes, likely in the patch-overlaps)
% merge neurons based on the correlation computed with {'A', 'S', 'C'}
cnmfe_quick_merge;            % run neuron merges

merge_thr = [0.6, 0.01, 0.01];  % choose thresholds for merging neurons (this will primarily merge neurons redundantly found by multiple patch processes, likely in the patch-overlaps)
% merge neurons based on the correlation computed with {'A', 'S', 'C'}
cnmfe_quick_merge;            % run neuron merges
% 

% merge_thr = [0.3, 0.3, 0.1];  % choose thresholds for merging neurons (this will primarily merge neurons redundantly found by multiple patch processes, likely in the patch-overlaps)
% % merge neurons based on the correlation computed with {'A', 'S', 'C'}
% cnmfe_quick_merge;            % run neuron merges
% 

neuron.Coor = neuron.get_contours(0.8); % energy  within the contour is 80% of the total

% cnmfe_yc_deleteNeurons_v3(neuron, [1:size(neuron.C,1)], neuron.C, Cn,pnr);
cnmfe_yc_deleteNeurons_locranked_v3(neuron, [1:size(neuron.C,1)], neuron.C, Cn,pnr);
neuron.Coor = neuron.get_contours(0.8); % energy  within the contour is 80% of the total
clear neuron_full neuron_bk neuron_full neuron_small neuron_ds RESULTS data Y
save(strcat(infile,'_small'));