addpath(genpath('/Users/hyeyoung/Documents/CODE/matnwb'))
addpath(genpath('d:\Users\USER\Documents\MATLAB\matnwb'))

datadir = 'S:\OpenScopeData\00248_v240130\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') );

probes = {'A', 'B', 'C', 'D', 'E'};

pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
load('')
nwbfiles = cat(1, dir([datadir nwbsessions{ises} filesep '*.nwb']), dir([datadir nwbsessions{ises} filesep '*' filesep '*.nwb']));
nwblfpfiles = nwbfiles(contains({nwbfiles.name}, 'probe'));

%% load lfp data for each probe
% for probeind = 0:5
fileind = find(contains({nwbfiles.name}, ['probe-' num2str(probeind)]));
nwblfpfile = fullfile([nwbfiles(fileind).folder filesep nwbfiles(fileind).name]);

nwblfp = nwbRead(nwblfpfile); 

% lfp = nwblfp.acquisition.get('probe_2_lfp_data');
% lfpdata = lfp.data.load();

lfp_probe = nwblfp.acquisition.get('probe_2_lfp');
lfp = lfp_probe.electricalseries.get('probe_2_lfp_data');
lfpdata = lfp.data.load(); % nprobes * ntimepoints
lfptimestaps = lfp.timestamps.load(); % ntimepoints * 1, units in seconds


%% figure out which electrode is where


%% align to vis stim

%% align to opto stim