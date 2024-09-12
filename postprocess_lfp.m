addpath(genpath('/Users/hyeyoung/Documents/CODE/matnwb'))

datadir = 'S:\OpenScopeData\00248_v240130\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') );

probes = {'A', 'B', 'C', 'D', 'E'};

nwbfiles = cat(1, dir([datadir nwbsessions{ises} filesep '*.nwb']), dir([datadir nwbsessions{ises} filesep '*' filesep '*.nwb']));
nwblfpfiles = nwbfiles(contains({nwbfiles.name}, 'probe'));

% for iprobe = 0:5
fileind = find(contains({nwbfiles.name}, ['probe-' num2str(iprobe)]));
nwblfpfile = fullfile([nwbfiles(fileind).folder filesep nwbfiles(fileind).name]);

nwblfp = nwbRead(nwblfpfile); 

%lfp_probe = nwblfp.acquisition.get('probe_2_lfp');
lfp = nwblfp.acquisition.get('probe_2_lfp_data');
lfpdata = lfp.data.load();
lfelectrodes = lfp.data.load();

nwblfp.processing
