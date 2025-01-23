%% load all data
%{
nwbsessions = {'sub-619293', 'sub-619296', 'sub-620333', 'sub-620334', ...
    'sub-625545', 'sub-625554', 'sub-625555', 'sub-630506', ...
    'sub-631510', 'sub-631570', 'sub-633229', 'sub-637484'};

trialorderacc = cell(numel(nwbsessions),1);
cnt = 0;
for ises = 1:numel(nwbsessions)
    cnt = cnt+1;
    mousedate = nwbsessions{ises};
    fprintf('%s %d\n', mousedate, ises)
    pathpp = ['S:\OpenScopeData\00248_v240130\postprocessed' filesep mousedate filesep];
    load([pathpp 'spkcnt_ICwcfg1_hireptt_V1RS_lmlv.mat'], 'trialorder')
    load([pathpp 'UMAP_V1RS_lmlvslopes.mat'])

    trialorderacc{ises} = trialorder;
    if cnt ==1
        UMAPorigagg = UMAPorig;
        UMAPorigallagg = UMAPorigall;
        UMAPorig_unsupagg = UMAPorig_unsup;
        UMAPorigall_unsupagg = UMAPorigall_unsup;
        UMAPlmlvacc = UMAPlmlv;
        UMAPlmlvallacc = UMAPlmlvall;
        UMAPlmlv_unsupacc = UMAPlmlv_unsup;
        UMAPlmlvall_unsupacc = UMAPlmlvall_unsup;
    else
        UMAPorigagg = cat(1, UMAPorigagg, UMAPorig);
        UMAPorigallagg = cat(1, UMAPorigallagg, UMAPorigall);
        UMAPorig_unsupagg = cat(1, UMAPorig_unsupagg, UMAPorig_unsup);
        UMAPorigall_unsupagg = cat(1, UMAPorigall_unsupagg, UMAPorigall_unsup);
        UMAPlmlvacc = cat(1, UMAPlmlvacc, UMAPlmlv);
        UMAPlmlvallacc = cat(1, UMAPlmlvallacc, UMAPlmlvall);
        UMAPlmlv_unsupacc = cat(1, UMAPlmlv_unsupacc, UMAPlmlv_unsup);
        UMAPlmlvall_unsupacc = cat(1, UMAPlmlvall_unsupacc, UMAPlmlvall_unsup);
    end

end

save('G:\My Drive\RESEARCH\logmean_logvar\OpenScope_UMAP_V1RS_lmlvslopes.mat', ...
    'lmlvslope_list', 'UMAPorigagg', 'UMAPorigallagg', 'UMAPorig_unsupagg', 'UMAPorigall_unsupagg', ...
    'UMAPlmlvacc', 'UMAPlmlvallacc', 'UMAPlmlv_unsupacc', 'UMAPlmlvall_unsupacc')
%}
%%
if ismac
    drivepath = '/Users/hyeyoung/Library/CloudStorage/GoogleDrive-shinehyeyoung@gmail.com/My Drive/';
    codepath = '/Users/hyeyoung/Documents/CODE/';
else
    drivepath = 'G:/My Drive/';
    codepath = 'C:\Users\USER\GitHub\';
end
addpath([codepath 'helperfunctions'])

load([drivepath 'RESEARCH/logmean_logvar/OpenScope_UMAP_V1RS_lmlvslopes.mat'])


