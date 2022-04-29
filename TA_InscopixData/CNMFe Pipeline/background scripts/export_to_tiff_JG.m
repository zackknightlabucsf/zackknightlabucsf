%for exporting files
clear all;close all;
addpath('C:\Program Files\Inscopix\Data Processing');

rec_names={ %list of videos to process and export in .isxd (must NOT include .isxd)
'2019-05-23-17-52-07_video-pp-BP-MC'
};
output_dir=cd; %output directory (by default "cd", or current directory)

rec_files = cellfun(@(x) fullfile(fileparts(which([x, '.isxd'])), [x, '.isxd']), rec_names, 'UniformOutput', false);

for p=1:length(rec_files)
    output=fullfile(fileparts(which([rec_names{p}])), [rec_names{p}, '.tiff']);
    isx.export_movie_to_tiff(rec_files{p}, output);
end
   