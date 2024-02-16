function opm_check_badch(sub,ses,project_dir)
%5 if running
% sub = '016'
% ses = '001'
%project_dir = 'R:\DRS-KidsOPM\Paediatric_OPM_Notts\';
addpath('R:\DRS-KidsOPM\Paediatric_OPM_Notts\fieldtrip-20220906')
ft_defaults
%%
%restoredefaultpath
cleaning_only = 0;
close all
clc
clearvars -except cleaning_only sub ses project_dir
run = 'run-001';

addpath(['R:\DRS-KidsOPM\Paediatric_OPM_Notts\fieldtrip-20220906'])
script_dir = mfilename('fullpath');fname = mfilename;script_dir = script_dir(1:end-length(fname));
addpath(script_dir)
addpath([script_dir,'Beamformer',filesep,''])
ft_defaults;
datadir = [project_dir,'Data',filesep,'BIDS',filesep];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exp_type = '_task-faces_circles';
filename = ['sub-',sub,'_ses-',ses,exp_type,'_',run];
path_main = [datadir,'sub-',sub,filesep,'ses-',ses,filesep];

path_ICA = [datadir,'derivatives',filesep,'ICA',filesep,'sub-',sub,filesep];
files_ICA = [path_ICA,filename];

path_cleaning = [datadir,'derivatives',filesep,'cleaning',filesep,'sub-',sub,filesep];
files_cleaning = [path_cleaning,filename];

path_VEs = [datadir,'derivatives',filesep,'VEs',filesep,'sub-',sub,filesep];
files_VE_AAL = [filename,'_VE_AAL'];
files_VE_cont = [filename,'_VE_cont'];

path_meg_data = [path_main,'meg',filesep];

path_meshes = [datadir,'derivatives',filesep,'sourcespace',filesep,'sub-',sub,filesep];
files_meshes = ['sub-',sub,'_meshes.mat'];
files_AAL_centroids = ['sub-',sub,'_AAL_centroids.nii.gz'];
files_AAL_regions = ['sub-',sub,'_AAL_regions.nii.gz'];
files_VOI = ['sub-',sub,'_AAL_VOI.mat'];
files_brain = ['sub-',sub,'_brain.nii'];

path_mri = [path_main,'anat',filesep];
files_mri = ['sub-',sub,'_anat.nii'];
S.mri_file = [path_mri,files_mri];

path_helmet = [datadir,'derivatives',filesep,'helmet',filesep,'sub-',sub,filesep];
files_helmet_info = dir([path_helmet,'*.mat']);files_helmet_info=files_helmet_info.name;

files_voxlox = ['sub-',sub,'_AAL_locs.mat'];
files_channels = [filename,'_channels.tsv'];

files_events = [filename,'_events.tsv'];

path_Tstat = [datadir,'derivatives',filesep,'Tstats',filesep,'sub-',sub,filesep]
files_Tstat = [filename,'_pseudoT_'];

path_TFS = [datadir,'derivatives',filesep,'VEs',filesep,'sub-',sub,filesep]
files_TFS_vis = [filename,'_TFS_vis'];
files_TFS_cont = [filename,'_TFS_cont']

%read meg data
cd(path_meg_data)
read_info = readlines([filename,'_meg_read_info.txt']);
Size = strsplit(read_info(1));Size = [str2num(Size(2)),str2num(Size(4))];

Precision = strsplit(read_info(2));Precision = Precision(2);

Ordering = strsplit(read_info(3));Ordering = Ordering(2);

FileID = fopen([path_meg_data,filename,'_meg.dat'],'r');

data=fread(FileID,Size,lower(Precision),Ordering)';
fclose(FileID);

% fs and other info
fID = fopen([path_meg_data,filename,'_meg.json']);
raw = fread(fID,inf);
json_info = jsondecode(char(raw'));
fs = json_info.SamplingFrequency;

% trigger info
event_table = readtable([path_meg_data,filename,'_events.tsv'],'FileType','text','delimiter','\t');
all_start_samps = [round(event_table(startsWith(event_table.type(:),'Circles'),:).sample)];

%helmet info
Helmet_info = [];
load([path_helmet,files_helmet_info])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Preproc
% Mean correct
data = data - mean(data,1);

hp = 1;
lp = 150;
[b,a] = butter(4,2*[hp lp]/fs);
data_f = [filtfilt(b,a,data)]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get rid of bad channels
%% Get rid of bad channels
if ~exist([path_meg_data,files_channels(1:end-4),'_proc.tsv'],'file')
    %     error("Use 'get_good_channels.m for this dataset first")
    ch_table = readtable([path_meg_data,files_channels],'FileType','text','Delimiter','tab');
    ch_table.isx = endsWith(ch_table.name,'X');
    ch_table.isy = endsWith(ch_table.name,'Y');
    ch_table.isz = endsWith(ch_table.name,'Z');
    ch_table.slot_no = zeros(height(ch_table),1);
    % sanity check
    if sum(sum([ch_table.isx,ch_table.isy,ch_table.isz],2)) ~= height(ch_table)
        error('Channel orientation [x,y,z] labels might be wrong!')
    end
    [ch_table] = Bad_Channels(data_f',ch_table,fs);
    writetable(ch_table,[path_meg_data,files_channels(1:end-4),'_proc.tsv'],...
        'WriteRowNames',true,'Delimiter','tab','FileType','text')
else
    ch_table = readtable([path_meg_data,files_channels(1:end-4),'_proc.tsv'],...
        'Delimiter','tab','FileType','text');
end

%%
% remove from data matrix
disp("Removing bad channels")
bad_chans_data = find(startsWith(ch_table.status,'bad'));

ch_table(bad_chans_data,:) = [];
data_f(bad_chans_data,:) = [];
Nchans = size(data_f,1);
save([path_cleaning 'sub-' sub '_Nchans.mat'],'Nchans')
end