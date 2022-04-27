%This script writes photometry data to .mat files stored within immediate subfolders
%And saves the .mat files to save_dir
%Written by: James Grove, on March 4, 2020

clear all; close all;
save_dir='C:\Users\ZK\Desktop\PRLH Analysis\FD Empty Bottle';%place to save mat files to

%%%%%%%%%%%%

%find current directory
base=pwd;
% Get a list of all files and folders in this folder.
files = dir;
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags);
subFolders(1)=[];subFolders(1)=[];
% Print folder names to command window.
for k = 1 : length(subFolders)
  cd([subFolders(k).folder '/' subFolders(k).name])
  if(~isempty(dir('*.tev')))
          names=split(cd, "\");
          d=TDTbin2mat(cd,'TYPE',{'epocs','streams'});
          cd(save_dir)
          save(names{length(names)},'d');
          fprintf(['wrote matfile for ' subFolders(k).name])
          fprintf('\n')
  else
      %find current directory
      base2=pwd;
      % Get a list of all files and folders in this folder.
      files2 = dir;
      % Get a logical vector that tells which is a directory.
      dirFlags2 = [files2.isdir];
      % Extract only those that are directories.
      subFolders2 = files2(dirFlags2);
      subFolders2(1)=[];subFolders2(1)=[];
      % Print folder names to command window.
      for k2 = 1 : length(subFolders2)
          cd([subFolders2(k2).folder '/' subFolders2(k2).name])
          if(~isempty(dir('*.tev')))
              try
              names=split(cd, "\");
              d=TDTbin2mat(cd,'TYPE',{'epocs','streams'});
              cd(save_dir)
              save(names{length(names)},'d');
              fprintf(['wrote matfile for ' subFolders2(k2).name])
              fprintf('\n')
              end
          end
      end
  end
end
cd(base)