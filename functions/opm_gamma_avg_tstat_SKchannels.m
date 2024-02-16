function opm_gamma_avg_tstat_SKchannels(sub,ses,project_dir)
%5 if running
close all
% sub = '125'
% ses = '001'
project_dir = 'R:\DRS-KidsOPM\Paediatric_OPM_Notts_AdultData\';
addpath('R:\DRS-KidsOPM\Paediatric_OPM_Notts\fieldtrip-20220906')
ft_defaults
%%
%restoredefaultpath
cleaning_only = 0;
close all
clc

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
files_ICA = [path_ICA,filename '_SK'];

path_cleaning = [datadir,'derivatives',filesep,'cleaning',filesep,'sub-',sub,filesep];
files_cleaning = [path_cleaning,filename];

path_VEs = [datadir,'derivatives',filesep,'VEs',filesep,'sub-',sub,filesep];
files_VEs = [filename,'_VE_unchopped_SK'];

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
files_Tstat = [filename,'_pseudoT_SK'];

path_TFS = [datadir,'derivatives',filesep,'VEs',filesep,'sub-',sub,filesep];
files_TFS = [filename,'_TFS_SK'];

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

% Notch filter
for harms = [50,100,150]
    Wo = harms/(fs/2);  BW = Wo/35;
    [b,a] = iirnotch(Wo,BW);
    disp('Applying Notch filter')
    data = filter(b,a,data,[],1);
end
Nchans = size(data,2);

% bandpass filter for viewing
disp('Applying 1-150 Hz bandpass filter')

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

%% Somewhat hacky way of getting slot numbers for later. Change ICA comp visualiser?

precision_limit = 1e-6;
pos = [ch_table.Px,ch_table.Py,ch_table.Pz];
for sl_i = 1:size(Helmet_info.lay.pos,1)
    detected_ind = find(sqrt(sum((repmat(Helmet_info.sens_pos(sl_i,:),height(ch_table),1) - pos/100).^2,2)) < precision_limit);
    if ~isempty(detected_ind)
        ch_table.slot_no(detected_ind) = sl_i;
    end
end
%%
% remove from data matrix
disp("Removing bad channels")
xchans = find(ch_table.isx);
ch_table.status(xchans)={'bad'};
sk_slot_nos = load('R:\DRS-KidsOPM\Paediatric_OPM_Notts_AdultData\Data\BIDS\derivatives\SK_slot_nos.mat');

for j = 1:size(ch_table,1)
    if ~ismember(ch_table.slot_no(j),sk_slot_nos.slot_no)
        ch_table.status(j) = {'bad'};
    end
end
bad_chans_data = find(startsWith(ch_table.status,'bad'));
ch_table(bad_chans_data,:) = [];
data_f(bad_chans_data,:) = [];

%% sensor info
S.sensor_info.pos = [ch_table.Px,ch_table.Py,ch_table.Pz]./100;
S.sensor_info.ors = [ch_table.Ox,ch_table.Oy,ch_table.Oz];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

xfm = load(sprintf('sub-%s_ses-%s%s_%s_sens2template_transform.txt',sub,ses,exp_type,run));
xfm = reshape(xfm,4,4);
xfm_rot = xfm(1:3,1:3);
xfm_translation = xfm(1:3,4)';
S.sensor_info.pos = (xfm_rot*S.sensor_info.pos' + xfm_translation')';
S.sensor_info.ors = (xfm_rot*S.sensor_info.ors')';
%% Mean field correction matrix
N = S.sensor_info.ors; % orientation matrix (N_sens x 3)
S.M = eye(length(N)) - N*pinv(N);

%% Epoch data
% segment using trigger
disp("Chopping data into epochs")
epoch_length = 3;
trig_offset = -1;
events_table = readtable([path_meg_data,files_events],'FileType','text','Delimiter','tab');
circles_events = startsWith(events_table.type,'Circles');
start_samples = events_table.sample(circles_events) + trig_offset*fs;
end_samples = start_samples + epoch_length*fs-1;
%check whether last trial is fully available and discard if not
unfinished_trls = find(end_samples > size(data_f,2));
start_samples(unfinished_trls) = [];
end_samples(unfinished_trls) = [];
circles_events(unfinished_trls) = [];
data_f_mat = zeros(size(data_f,1),epoch_length*fs,size(end_samples,1));
for tr_i = 1:size(end_samples,1)
    data_f_mat(:,:,tr_i) = data_f(:,start_samples(tr_i):end_samples(tr_i));
end

clear data_f
%%%%%

%% Put data in FT format
addpath('C:\Users\ppynr2\OneDrive - The University of Nottingham\Documents\GitHub\Gamma_circles_development\functions')
disp("Converting to FieldTrip format")
[data_strct] = makeFTstruct(data_f_mat,fs,ch_table,S.sensor_info);
%% trial info
data_strct = ft_checkdata(data_strct, 'datatype', 'raw', 'hassampleinfo', 'yes');

trl_mat = sampleinfo2trl(data_strct);
trial_info = table(trl_mat(:,1),trl_mat(:,2),trl_mat(:,3),'VariableNames',{'start','end','offset'});
trial_info.type = events_table.type(circles_events);
removefields(data_strct,'sampleinfo')
%% resample for viewing and mark artefacts
disp("Check for bad trials")
if ~exist([files_cleaning,'_vis_artfcts.mat'],'file')
    resamplefs = 150; %Hz
    cfg            = [];
    cfg.resamplefs = resamplefs;
    cfg.detrend    = 'no';
    data_preproc_150   = ft_resampledata(cfg, data_strct);

    %%% Get rid of bad trials
    cfg_art          = [];
    cfg_art.viewmode = 'vertical';
    cfg_art = ft_databrowser(cfg_art,data_preproc_150);
    clear data_preproc_150

    vis_artfcts = cfg_art.artfctdef.visual.artifact * ...
        (data_strct.fsample/resamplefs);
    save([files_cleaning,'_vis_artfcts.mat'],'vis_artfcts')
else
    load([files_cleaning,'_vis_artfcts.mat'],'vis_artfcts')
end

% automatic artifact rejection
thresh_val = 3;
auto_artfcts = get_bad_segments(data_f_mat,thresh_val);
fprintf('Found %d artifacts using a threshold of %d std. deviations.\n',...
    size(auto_artfcts,1),thresh_val);

% combine artifacts and reject bad trials
cfg = [];
cfg.artfctdef.visual.artifact = [vis_artfcts;auto_artfcts];
cfg.artfctdef.reject  = 'complete';
data_vis_clean = ft_rejectartifact(cfg,data_strct);

fprintf('\nRejected %d of %d epochs of length %1.2f s.\n',...
    size(data_strct.trial,2)-size(data_vis_clean.trial,2),size(data_strct.trial,2),epoch_length);

good_trials = false(height(trial_info),1);
for clean_trial_ind = 1:size(data_vis_clean.cfg.trl,1)
    good_trials(sum(data_vis_clean.cfg.trl(clean_trial_ind,1:2)==trial_info{:,1:2},2)==2) = true;
end
% only include good trials in info table
trial_info = trial_info(good_trials,:);
%% ICA
% lay = Helmet_info.lay;
% disp("ICA artifact rejection")
% if ~exist([files_ICA,'_bad_ICA_comps.mat'],'file') || ~exist([files_ICA,'_ICA_data.mat'],'file')
% 
%     % Resample for faster ICA
%     cfg            = [];
%     cfg.resamplefs = 150;
%     cfg.detrend    = 'no';
%     data_ica_150   = ft_resampledata(cfg, data_vis_clean);
% 
%     % Run ICA on 150 Hz data or load previous unmixing matrix
%     cfg            = [];
%     if ~exist([files_ICA,'_ICA_data.mat'],'file')
%         cfg.method = 'runica';
%     else
%         load([files_ICA,'_ICA_data.mat'],'comp150')
%         cfg.unmixing   = comp150.unmixing;
%         cfg.topolabel  = comp150.topolabel;
%     end
%     comp150    = ft_componentanalysis(cfg, data_ica_150);
% 
%     % Inspect components for rejection or load file with bad component list
%     if ~exist([files_ICA,'_bad_ICA_comps.mat'],'file')
%         disp("Choose bad components")
%         [bad_comps] = plot_ICA_comps(comp150,ch_table,lay,[]);
% 
%         save([files_ICA,'_bad_ICA_comps.mat'],'bad_comps')
%         close(gcf)
%     else
%         disp("Loading saved bad components")
%         load([files_ICA,'_bad_ICA_comps.mat'],'bad_comps')
%     end
% 
%     % only keep unmixing matrix and topolabel for component removal
%     tokeep = {'unmixing','topolabel'};
%     fns=fieldnames(comp150);
%     toRemove = fns(~ismember(fns,tokeep));
%     comp150 = rmfield(comp150,toRemove);
%     save([files_ICA,'_ICA_data.mat'],'comp150')
% 
%     clear data_ica_150
% else
%     disp("Loading bad coponents, topographies and old unmixing matrix")
%     load([files_ICA,'_bad_ICA_comps.mat'],'bad_comps')
%     load([files_ICA,'_ICA_data.mat'],'comp150')
% end
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     Perform ICA on original 1200 Hz data using unmixing matrix and topolabel
% cfg            = [];
% cfg.unmixing   = comp150.unmixing;
% cfg.topolabel  = comp150.topolabel;
% comp1200        = ft_componentanalysis(cfg, data_vis_clean);
% 
% % Plot comps again to confirm they are correct
% disp("Confirm bad ICA components")
% % [bad_comps] = plot_ICA_comps(comp1200,ch_table,lay,bad_comps)
% save([files_ICA,'_bad_ICA_comps.mat'],'bad_comps');
% 
% % Remove components from data
% cfg           = [];
% cfg.component = bad_comps;
% data_ica_clean    = ft_rejectcomponent(cfg, comp1200,data_vis_clean);
data_ica_clean = data_vis_clean;N_clean_trls = size(data_ica_clean.trial,2);
clear data_vis_clean
%% Reconstitute data from FT structure

data_f_clean = [data_ica_clean.trial{1,:}];
clear data_ica_clean
% apply mean field correction
disp("Applying mean field correction")
data_f_clean = S.M*data_f_clean;

%% Source positions

% load AAL locations
VOI_mat_file = [path_meshes,files_VOI];
brain_1mm = ft_read_mri([path_meshes,files_brain]);

cfg = [];
cfg.downsample = 4;
[downsample_brain] = ft_volumedownsample(cfg, brain_1mm);
downsample_brain = ft_convert_units(downsample_brain,'m');
ft_write_mri([path_meshes,files_brain(1:end-4),'_4mm.nii'],downsample_brain,...
    'dataformat','nifti');
VOI = downsample_brain.anatomy > 0;

[sourcepos_vox(:,1),sourcepos_vox(:,2),sourcepos_vox(:,3)] = ind2sub(downsample_brain.dim, find(VOI));

sourcepos = ft_warp_apply(downsample_brain.transform,sourcepos_vox);

if ~cleaning_only
    
    [bf_outs_shell] = run_beamformer('shell',sourcepos,S,0,[],1);
    lead_fields_shell_xyz = bf_outs_shell.LF;
    
    outside_locs = squeeze(sum(lead_fields_shell_xyz(1,:,:)==0)==3);
    lead_fields_shell_xyz(:,:,outside_locs) = [];
    sourcepos(outside_locs,:)=[];
    sourcepos_vox(outside_locs,:)=[];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% convert orientation of sources to polar
    % Load meshes
    
    load([path_meshes,files_meshes],'meshes');
    
    meshes = ft_convert_units(meshes,'m');
    
    X = meshes.pnt; Origin = mean(X,1);
    Ndips = length(sourcepos);
    for n = 1:Ndips
        thispos = sourcepos(n,:);
        [phi,theta1,~] = cart2sph(thispos(1) - Origin(1),thispos(2) - Origin(2) ,thispos(3) - Origin(3));
        theta = pi/2 - theta1;
        Src_Or_theta(n,:) = [cos(theta)*cos(phi) cos(theta)*sin(phi) -sin(theta)];
        Src_Or_phi(n,:) = [-sin(phi) cos(phi) 0];
        Lead_fields(:,1,n) = Src_Or_theta(n,1)*lead_fields_shell_xyz(:,1,n) + Src_Or_theta(n,2)*lead_fields_shell_xyz(:,2,n) + Src_Or_theta(n,3)*lead_fields_shell_xyz(:,3,n);
        Lead_fields(:,2,n) = Src_Or_phi(n,1)*lead_fields_shell_xyz(:,1,n) + Src_Or_phi(n,2)*lead_fields_shell_xyz(:,2,n) + Src_Or_phi(n,3)*lead_fields_shell_xyz(:,3,n);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% make a plot of the geometry...
    figure(1);
    ft_plot_mesh(meshes,'facecolor',[.5 .5 .5],'facealpha',.3,'edgecolor','none')
    hold on
    scatter3(sourcepos(:,1),sourcepos(:,2),sourcepos(:,3),'ro','linewidth',3)
    view([130,30])
    fig = gcf;
    fig.Color = [1,1,1];
    plot3(S.sensor_info.pos(:,1),S.sensor_info.pos(:,2),S.sensor_info.pos(:,3),'o')
    quiver3(S.sensor_info.pos(ch_table.isx==1,1),S.sensor_info.pos(ch_table.isx==1,2),S.sensor_info.pos(ch_table.isx==1,3),...
        S.sensor_info.ors(ch_table.isx==1,1),S.sensor_info.ors(ch_table.isx==1,2),S.sensor_info.ors(ch_table.isx==1,3),'r','linewidth',2)
    quiver3(S.sensor_info.pos(ch_table.isy==1,1),S.sensor_info.pos(ch_table.isy==1,2),S.sensor_info.pos(ch_table.isy==1,3),...
        S.sensor_info.ors(ch_table.isy==1,1),S.sensor_info.ors(ch_table.isy==1,2),S.sensor_info.ors(ch_table.isy==1,3),'g','linewidth',2)
    quiver3(S.sensor_info.pos(ch_table.isz==1,1),S.sensor_info.pos(ch_table.isz==1,2),S.sensor_info.pos(ch_table.isz==1,3),...
        S.sensor_info.ors(ch_table.isz==1,1),S.sensor_info.ors(ch_table.isz==1,2),S.sensor_info.ors(ch_table.isz==1,3),'b','linewidth',2)
    plot3(Origin(1),Origin(2),Origin(3),'bo','linewidth',4)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% take a random lead field and plot it...
    figure(2);
    ft_plot_mesh(meshes,'facecolor',[.5 .5 .5],'facealpha',.3,'edgecolor','none')
    hold on
    ft_plot_topo3d(double(S.sensor_info.pos(ch_table.isz==1,:)),Lead_fields(ch_table.isz==1,2,16))
    alpha(gca,0.5)
    plot3(S.sensor_info.pos(ch_table.isx==1,1),S.sensor_info.pos(ch_table.isx==1,2),S.sensor_info.pos(ch_table.isx==1,3),'go','linewidth',3)
    
    [k2,~] = convhull(sourcepos,'Simplify',true);
    trisurf(k2,sourcepos(:,1),sourcepos(:,2),sourcepos(:,3),'FaceColor','m','FaceAlpha',0.5,'EdgeColor','none')
    scatter3(sourcepos(16,1),sourcepos(16,2),sourcepos(16,3),'r','linewidth',4)
    % quiver3(S.sensor_info.pos(ch_table.isx==1,1),S.sensor_info.pos(ch_table.isx==1,2),S.sensor_info.pos(ch_table.isx==1,3),...
    %     S.sensor_info.ors(ch_table.isx==1,1).*Lead_fields(ch_table.isx==1,2,16) + S.sensor_info.ors(ch_table.isy==1,1).*Lead_fields(ch_table.isy==1,2,16) + S.sensor_info.ors(ch_table.isz==1,1).*Lead_fields(ch_table.isz==1,2,16),...
    %     S.sensor_info.ors(ch_table.isx==1,2).*Lead_fields(ch_table.isx==1,2,16) + S.sensor_info.ors(ch_table.isy==1,2).*Lead_fields(ch_table.isy==1,2,16) + S.sensor_info.ors(ch_table.isz==1,2).*Lead_fields(ch_table.isz==1,2,16),...
    %     S.sensor_info.ors(ch_table.isx==1,3).*Lead_fields(ch_table.isx==1,2,16) + S.sensor_info.ors(ch_table.isy==1,3).*Lead_fields(ch_table.isy==1,2,16) + S.sensor_info.ors(ch_table.isz==1,3).*Lead_fields(ch_table.isz==1,2,16),'r','linewidth',2)
    quiver3(S.sensor_info.pos(ch_table.isz==1,1),S.sensor_info.pos(ch_table.isz==1,2),S.sensor_info.pos(ch_table.isz==1,3),...
        S.sensor_info.ors(ch_table.isz==1,1).*Lead_fields(ch_table.isz==1,2,16),...
        S.sensor_info.ors(ch_table.isz==1,2).*Lead_fields(ch_table.isz==1,2,16),...
        S.sensor_info.ors(ch_table.isz==1,3).*Lead_fields(ch_table.isz==1,2,16),'r','linewidth',2)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clearvars -except path_* files_* sub ses epoch_length data_f_clean cleaning_only fs sourcepos sourcepos_vox Lead_fields N_clean_trls trial_info meshes trig_offset
    %%
    %% filter the OPM data to band of interest
    hp = 30;
    lp = 80;
    %         if ~exist(sprintf('%s%s_%d_%d_Hz_Z.mat',path_AEC,files_AEC,hp,lp),'file')
    [b,a] = butter(4,2*[hp lp]/fs);
    data_f = [filtfilt(b,a,data_f_clean')]';
    data_f_mat = reshape(data_f,size(data_f,1),[],height(trial_info));
    
    
    active_window = round([0.3 1].*fs-trig_offset*fs);active_inds = active_window(1):active_window(2);
    control_window = round([-0.8 -0.1].*fs-trig_offset*fs);control_inds = control_window(1):control_window(2);
    
    %index
    circles_trials = data_f_mat;
    N_circles_trials = size(circles_trials,3);
    
    C_circles = cov(reshape(circles_trials,size(circles_trials,1),prod(size(circles_trials,2,3)))');
    
    Ca_circles = cov(reshape(circles_trials(:,active_inds,:),...
        size(circles_trials(:,active_inds,:),1),prod(size(circles_trials(:,active_inds,:),2,3)))');
    
    Cc_circles = cov(reshape(circles_trials(:,control_inds,:),...
        size(circles_trials(:,control_inds,:),1),prod(size(circles_trials(:,control_inds,:),2,3)))');
    
    %% Beamform
    pseudoT = get_pseudoT(C_circles,Ca_circles,Cc_circles,Lead_fields);
    
    cax_val = [0.5*max(pseudoT), max(pseudoT)];
    [~, max_voxel] = max(pseudoT);
    fg = figure;
    fg.Name = 'pseudoT'
    ft_plot_mesh(meshes,'facecolor',[.5 .5 .5],'facealpha',.1,'edgecolor','none')
    hold on
    scatter3(sourcepos(:,1),sourcepos(:,2),sourcepos(:,3),50,pseudoT,'filled','MarkerFaceAlpha',0.2)
    colormap hot;caxis([cax_val]);cb = colorbar;cb.Label.String = 'Tstat';
    axis equal
    view([-120 20])
    fig = gcf;
    fig.Color = 'w';
    ax = gca;
    ax.FontSize = 14;
    
    % Find voxel with min/max value of Tstat
    dip_loc_circles = sourcepos(max_voxel,:);
    plot3(dip_loc_circles(1),dip_loc_circles(2),dip_loc_circles(3),'gx','markersize',10,'LineWidth',5)
    
    % get peak voxel timecourse

    %%
    [mean_Env_circles] = get_envelope(circles_trials,C_circles,Lead_fields,max_voxel,control_inds);
    figure
    set(gcf,'Color','w');

    trial_time = linspace(0+trig_offset,size(mean_Env_circles,1)./fs+trig_offset,size(mean_Env_circles,1));
    plot(trial_time,mean_Env_circles,'r','Displayname','Circles peak envelope');
    
    %%
    % save timecourses
    %save([path_Tstat,'circles_peak_env.mat'],'mean_Env_circles')
    
    dip_loc_circles_mm = dip_loc_circles.*1000;
    %save([path_Tstat,'circles_peak.txt'],'dip_loc_circles_mm','-ascii','-double')
    %% save Tstat maps

    brain_1mm = ft_read_mri([path_meshes,files_brain]);
    
    cfg = [];
    cfg.downsample = 4;
    [downsample_brain] = ft_volumedownsample(cfg, brain_1mm);
    downsample_brain = ft_convert_units(downsample_brain,'m');
    positions = sourcepos_vox;%ft_warp_apply(inv(downsample_brain.transform),sourcepos_vox);

    image = downsample_brain;
    image.anatomy = image.anatomy.*0;
    image.anatomy(sub2ind(downsample_brain.dim,...
        positions(:,1),positions(:,2),positions(:,3))) = pseudoT;
    ft_write_mri([path_Tstat,files_Tstat,'circles.nii'],image,...
        'dataformat','nifti');
end
end 
function trl = sampleinfo2trl(data)
% borrowed from private FieldTrip funcs
% SAMPLEINFO2TRL constructs the trial definition from the sampleinfo, the time axes
% and optionally from the trialinfo
%
% Use as
%   trl = sampleinfo2trl(data)
%

% get the begin and end sample of each trial
begsample = data.sampleinfo(:,1);
endsample = data.sampleinfo(:,2);

% recreate the offset
offset = zeros(numel(data.trial), 1);
for i=1:numel(data.trial)
    offset(i) = round(data.time{i}(1)*data.fsample);
end

if isfield(data, 'trialinfo') && istable(data.trialinfo)
    trl = table(begsample, endsample, offset);
    trl = horzcat(trl, data.trialinfo);
elseif isfield(data, 'trialinfo') && isnumeric(data.trialinfo)
    trl = [begsample endsample offset data.trialinfo];
else
    trl = [begsample endsample offset];
end
end

function pseudoT = get_pseudoT(C,Ca,Cc,Lead_fields)
% Calculate pseudo T statistic
mu = 0.05;
Cr = C + mu*max(svd(C))*eye(size(C));
Cr_inv = inv(Cr);
pseudoT = zeros(1,size(Lead_fields,3));
infolength = 0;
for n = 1:size(Lead_fields,3)
    this_L = Lead_fields(:,:,n);
    W_v = inv((this_L'*Cr_inv*this_L))*(this_L'*Cr_inv);
    iPower_v = this_L'*Cr_inv*this_L;
    [v,d] = svd(iPower_v);
    [~,id] = min(diag(d));
    lopt = this_L*v(:,id); % turn to nAm amplitude
    w = (lopt'*Cr_inv/(lopt'*Cr_inv*lopt))';
    pseudoT(:,n) = (w'*Ca*w - w'*Cc*w)./(0.5.*(w'*Ca*w + w'*Cc*w));
    
    fprintf(repmat('\b',1,infolength));
    infolength = fprintf('Evaluating pseudo-T (%d/%d)\n',n,size(Lead_fields,3));
end
end

function [mean_Env] = get_envelope(trial_data,C,Lead_fields,min_voxel_index,control_inds)
mu = 0.05;
Cr = C + mu*max(svd(C))*eye(size(C));
Cr_inv = inv(Cr);

this_L = Lead_fields(:,:,min_voxel_index);
W_v = inv((this_L'*Cr_inv*this_L))*(this_L'*Cr_inv);
iPower_v = this_L'*Cr_inv*this_L;
[v,d] = svd(iPower_v);
[~,id] = min(diag(d));
lopt = this_L*v(:,id); % turn to nAm amplitude
w = (lopt'*Cr_inv/(lopt'*Cr_inv*lopt))';
N_trials = size(trial_data,3);
VEs = zeros(size(trial_data,2),N_trials);
Envs = VEs;
for tr_i = 1:N_trials
    VEs(:,tr_i) = w'*trial_data(:,:,tr_i) ./sqrt(w'*w);
    Envs(:,tr_i) = abs(hilbert(VEs(:,tr_i)));
    Envs(:,tr_i) = (Envs(:,tr_i) - mean(Envs(control_inds,tr_i)))./mean(Envs(control_inds,tr_i));
end
mean_Env = mean(Envs,2);
end