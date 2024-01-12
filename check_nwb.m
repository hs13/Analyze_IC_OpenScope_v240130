% requested by Ahad: Could you look at sub-630506_ses-1192952692_ogen.nwb just to confirm everything looks correct?
addpath(genpath('d:\Users\USER\Documents\MATLAB\matnwb'))
addpath('C:\Users\USER\GitHub\Analyze_IC_OpenScope_v230821')
addpath(genpath('C:\Users\USER\GitHub\helperfunctions'))

datadir = 'S:\OpenScopeData\00248_v230821\';

% whichnwbsession = 'sub-630506';
% nwbfiles = cat(1, dir([datadir whichnwbsession filesep '*.nwb']), dir([datadir whichnwbsession filesep '*' filesep '*.nwb']));
% % take filename  with shortest length or filename that does not contain probe
% [~, fileind] = min(cellfun(@length, {nwbfiles.name}));
% nwbspikefile = fullfile([nwbfiles(fileind).folder filesep nwbfiles(fileind).name]);
% % nwbspikefile = string(nwbspikefile);
% disp(nwbspikefile)

nwbspikefile = "G:\My Drive\RESEARCH\sub-625554_ses-1181330601-acq-FINAL_ogen.nwb";
nwb = nwbRead(nwbspikefile); 

