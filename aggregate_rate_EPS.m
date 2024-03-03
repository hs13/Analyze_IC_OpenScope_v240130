% info to save
% Ravg (per trialtype)
% nwbsessions, sesgenotype
% neuinfo: session, RS, FS, area, OT
% neugroup: IC1-encoders, IC2-encoders, LC1-encoders, LC2-encoders, RE1-faith, RE2-faith, indenc1, indenc2, indenc3, indenc4, non-responsive


addpath('C:\Users\USER\GitHub\Analyze_IC_OpenScope_v240130')
addpath(genpath('C:\Users\USER\GitHub\helperfunctions'))

datadir = 'S:\OpenScopeData\00248_v240130\';
nwbdir = dir(datadir);
nwbsessions = {nwbdir.name};
nwbsessions = nwbsessions( contains(nwbsessions, 'sub-') | contains(nwbsessions, 'sub_') );

ICblocks = {'ICwcfg1_presentations','ICwcfg0_presentations','ICkcfg1_presentations','ICkcfg0_presentations'};
Ravgagg = struct();
neuinfoagg = struct();
neugroupsoi = {'all', 'ICencoder1', 'ICencoder2', 'LCencoder1', 'LCencoder2', ...
    'RElfaith1', 'RElfaith2', 'indenc1', 'indenc2', 'indenc3', 'indenc4', ...
    'indresp1', 'indresp2', 'indresp3', 'indresp4', 'nonresponsive'};
neugroupagg = struct();
for ises = 1:numel(nwbsessions)
    clearvars -except datadir nwbsessions Ravgagg neuinfoagg neugroupagg
    pathpp = [datadir 'postprocessed' filesep nwbsessions{ises} filesep];
    load([pathpp 'qc_units.mat'])
    load([pathpp 'visresponses.mat'])
    load([pathpp 'postprocessed.mat'])
end