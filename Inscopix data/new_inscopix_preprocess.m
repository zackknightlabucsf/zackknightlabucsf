%Inscopix Preprocessing for nVista 3.0 and nVoke 2.0 .isxd files
%requires new data processing software
%Adapted from example API by James Grove, 6/17/19

%Brief explanation:
%1. First, you will be given a mean projection of the video for cropping
%draw a square and double click it when done
%2. This cropped video will then be temporally+spatially downsampled
%3. Then it will spatially bandpass filtered
        % essentially a smoothed frame will be subtracted from a less smoothed frame
        % smoothing is done using 2D Gaussian kernels by default standard deviations
        % this removes out of focus cells/fluorescence (highpass)
        % and smoothes the video (lowpass)
%4. Then you will be given another mean projection to crop for motion corr
        %draw a polygon to envelop the area and double click it when done
%5. This video will then be motion corrected and cropped
%6. Then it will be saved to a .tiff file for cnmfe
%7. Double check files in ImageJ
%Further explanation: https://support.inscopix.com/inscopix-data-processing-121-user-guide-html

close all; clear all;
rec_names={ %list of videos to process and export in .isxd (must NOT include .isxd)
'JG_111121_DV2ex1'
'JG_102621_DV2ex2'
'JG_111121_DV2ex2'
'JG_111121_DVDL102'
'JG_111121_DVDL103'
'JG_111121_DVDL101'
'JG_111021_DV2ex1'
'JG_110921_DVDL102'
'JG_110921_DVDL103'
'JG_110821_DV2ex5'
'JG_110821_DV2ex2'
'JG_110821_DVDL102'
'JG_110821_DVDL103'
'JG110521_DVDL102'
'JG_110321_DV2ex5'
'JG_110321_DV2ex6'
'JG_110321_DV2ex2'
'JG_110221_DV2ex5'
'JG_110221_DV2ex2'
'JG_110121_DV2ex5'
'JG_110121_DV2ex2'
'JG_103021_DV2ex5'
'JG_103021_DV2ex2'
'JG_102821_DV2ex2'
'JG_102721_DV2ex5'
'JG_102721_DV2ex2'
'JG_102621_DV2ex5'
};

sdf=2; %spatial downsample
tdf=2; %temporal downsample
output_dir=cd; %output directory (by default "cd", or current directory)

%% Start up
addpath('C:\Program Files\Inscopix\Data Processing'); %path to inscopix data processing package

%% Generate the recording file paths.
rec_files = cellfun(@(x) fullfile(fileparts(which([x, '.isxd'])), [x, '.isxd']), rec_names, 'UniformOutput', false);

% Preprocess the recordings by spatially downsampling and temporally downsampling
pp_files = isx.make_output_file_paths(rec_files, output_dir, 'pp');
disp('starting video projection');
project_files = isx.make_output_file_paths(rec_files, output_dir, 'mean1');

for p=1:length(pp_files)
    isx.project_movie(rec_files{p},project_files{p}, 'stat_type', 'mean');
end

disp('please crop the videos');
roi={};
for p=1:length(pp_files)
    movie = isx.Movie.read(project_files{p});
    image = movie.get_frame_data(0);
    clear movie
    figure
    imagesc(image);
    colormap hot;
    caxis([max(prctile(image,0)) max(prctile(image,99))])
    h=imrect;
    rect=wait(h);
    rect=round(rect);
    newrect=[rect(2),rect(1),rect(2)+rect(4),rect(3)+rect(1)];
    roi{p}=newrect;
end
close all
disp('starting preprocessing');
for p=1:length(pp_files)
    isx.preprocess(rec_files{p}, pp_files{p}, 'spatial_downsample_factor', sdf, 'temporal_downsample_factor', tdf,'crop_rect',roi{p});
end
%% Spatial bandpass
disp('starting spatial bandpass')
bp_files = isx.make_output_file_paths(pp_files, output_dir, 'BP');

for p=1:length(pp_files)
    isx.spatial_filter(pp_files{p}, bp_files{p}, 'low_cutoff', 0.005, 'high_cutoff', 0.500);
end

%% Motion correction
disp("starting motion correction")
mean_proj_file = isx.make_output_file_paths(bp_files, output_dir, 'mean2');
mc_files = isx.make_output_file_paths(bp_files, output_dir, 'MC');
crop_rect_file =  isx.make_output_file_paths(bp_files, output_dir, 'cropped');
translation_files = isx.make_output_file_paths(mc_files, output_dir, 'translations', 'ext', '.csv');
roi={};
for p=1:length(pp_files) %for loop used bc parfor was having problems
    isx.project_movie(bp_files{p}, mean_proj_file{p}, 'stat_type', 'mean');
end

for p=1:length(pp_files)
    movie = isx.Movie.read(bp_files{p});
    image = movie.get_frame_data(0);
    clear movie
    figure
    imagesc(image);
    colormap hot;
    [~,xi,yi]=roipoly;
    roi{p}=round([xi,yi]);
end
close all

for p=1:length(pp_files) %for loop used bc parfor was having problems
    isx.motion_correct( ...
        bp_files{p}, mc_files{p}, 'max_translation', 20, 'reference_file_name', mean_proj_file{p}, ...
        'low_bandpass_cutoff', 0, 'high_bandpass_cutoff', 0, 'roi', roi{p}, ...
        'output_translation_files', translation_files{p}, ...
        'output_crop_rect_file', crop_rect_file{p});
end

%% Export to tiff
disp('exporting to tiff')

for p=1:length(rec_names)
    rec_names{p}=[rec_names{p}, '-pp-BP-MC']
end

rec_files = cellfun(@(x) fullfile(fileparts(which([x, '.isxd'])), [x, '.isxd']), rec_names, 'UniformOutput', false);

for p=1:length(rec_files)
    output=fullfile(fileparts(which([rec_names{p}])), [rec_names{p}, '.tiff']);
    isx.export_movie_to_tiff(rec_files{p}, output);
end
   

disp('ready for cnmfe')
