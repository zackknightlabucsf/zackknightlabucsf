%Inscopix Preprocessing for nVista 3.0 and nVoke 2.0 .isxd files
%requires new data processing software

%Don't need startup_all

close all; clear all;
rec_names={ %list of videos to process and export in .isxd (must NOT include .isxd)

'2022-03-24-12-32-06_video' %TA146_fastedStimSuc_0324
'2022-03-24-13-34-09_video' %ND05_fastedStimSuc_0324
'2022-03-24-12-34-31_video' %TA145_fastedStimSuc_0324
'2022-03-24-13-39-22_video' %TA186_fastedStimSuc_0324
'2022-03-24-14-40-12_video' %TA175_fastedStimSuc_0324
'2022-03-24-15-39-40_video' %TA176_fastedStimSuc_0324


};

sdf=2; %spatial downsample
tdf=5; %temporal downsample
output_dir=cd; %output directory (by default "cd", or current directory)
already_preprocessed=false; %if true, the files should be in .isxd format

%% Start up
addpath('C:\Program Files\Inscopix\Data Processing');

%% Generate the recording file paths.
rec_files = cellfun(@(x) fullfile(fileparts(which([x, '.isxd'])), [x, '.isxd']), rec_names, 'UniformOutput', false);
%rec_files = cellfun(@(x) fullfile(fileparts(which([x, '.tiff'])), [x, '.tiff']), rec_names, 'UniformOutput', false);

%% Preprocess the recordings by spatially downsampling and temporally downsampling
if(already_preprocessed)
    pp_files=rec_files
else
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
end
disp('starting preprocessing');
for p=1:length(pp_files)
    isx.preprocess(rec_files{p}, pp_files{p}, 'spatial_downsample_factor', sdf, 'temporal_downsample_factor', tdf,'crop_rect',roi{p});
end
disp('starting spatial bandpass')
bp_files = isx.make_output_file_paths(pp_files, output_dir, 'BP');

for p=1:length(pp_files)
    isx.spatial_filter(pp_files{p}, bp_files{p}, 'low_cutoff', 0.005, 'high_cutoff', 0.500);
end

disp('starting motion correction')
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

