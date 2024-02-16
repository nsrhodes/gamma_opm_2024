%% Housekeeping
clear all
close all
clc

%% Fieldtrip set up

addpath('R:\DRS-KidsOPM\Paediatric_OPM_Notts\fieldtrip-20220906')
ft_defaults

%% Select participants

base_dir = 'R:\DRS-KidsOPM\';
load('C:\Users\ppynr2\OneDrive - The University of Nottingham\phd\Gamma\demographics_alldata.mat')
proj_dir_nk = [base_dir 'Paediatric_OPM_Notts\Data\BIDS\derivatives\VEs\'];
proj_dir_na = [base_dir 'Paediatric_OPM_Notts_AdultData\Data\BIDS\derivatives\VEs\'];
proj_dir_sk = [base_dir 'SickKids\gamma_VEs\'];


%% File set-up
script_dir = 'C:\Users\ppynr2\OneDrive - The University of Nottingham\Documents\GitHub\Gamma_circles_development\';
addpath(script_dir)
addpath([script_dir 'functions'])

%% Read in all TFS and VEs

% Organise into data table
data = demographics; TFS = cell(size(data,1),1); VE = cell(size(data,1),1);N_trials = cell(size(data,1),1);
data = addvars(data,TFS);data = addvars(data,VE);data = addvars(data,N_trials);
inds_nk=[];inds_na=[];inds_sk=[];inds_sa=[];
count = 0;
% Define TFS parameters
highpass = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39 41 43 45 47 49 51 53 55 57 59 61 63 65 67 69 71 73 75 77 79 81 83 85 87 89];
lowpass = [5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39 41 43 45 47 49 51 53 55 57 59 61 63 65 67 69 71 73 75 77 79 81 83 85 87 89 91 93];
fre = highpass + ((lowpass - highpass)./2);
trig_offset = -1;
duration = 3;
fs = 1200;
trial_time = linspace(0+trig_offset,size(TFS,2)./fs+trig_offset,size(TFS,2));
cbar_lim = 0.3;
OFF = [-0.8 0];
% Control window
conwin = OFF*fs-trig_offset*fs;
for sub_i = 1:size(demographics,1)
    sub = demographics.subject{sub_i}
    if sub(1)=='0' %nottingham kids
        project_dir = proj_dir_nk;
        inds_nk = [inds_nk sub_i];
    elseif sub(1)=='1' %nottingham adults
        project_dir = proj_dir_na;
        inds_na = [inds_na sub_i];
    elseif sub(1)=='2'%sk kids
        project_dir = proj_dir_sk;
        inds_sk = [inds_sk sub_i];
    elseif sub(1)=='3' %sk adults
        project_dir = proj_dir_sk;
        inds_sa = [inds_sa sub_i];
    end
    % Load in TFS
%     load([project_dir 'sub-' sub '\sub-' sub '_ses-001_task-faces_circles_run-001_TFS_vis.mat'])
%     data.TFS_vis{sub_i} = TFS;
%     load([project_dir 'sub-' sub '\sub-' sub '_ses-001_task-faces_circles_run-001_TFS_cont.mat'])
%     data.TFS_cont{sub_i} = TFS;
    % Load in virtual electrode
    load([project_dir 'sub-' sub '\sub-' sub '_ses-001_task-faces_circles_run-001_VE_AAL.mat'])
    VEs_AAL{sub_i} = VE_AAL;
    data.VE_AAL{sub_i} = VE_AAL;
    VEs_vis{sub_i} = VE_AAL(:,65)';
    VEs_cont{sub_i} = VE_AAL(:,14)';
    inds{sub_i} = ind;
    data.N_trials{sub_i} = length(ind);
    Ntrials = length(ind);

    % Filter data within bands and calculate envelope VIS
    VE_fb = zeros(length(VEs_vis{sub_i}),length(fre));
    for fb = 1:length(highpass)
        filt_VE = nut_filter3(VEs_vis{sub_i}','butter','bp',4,highpass(fb),lowpass(fb),fs,1)';
        VE_fb(:,fb) = abs(hilbert(filt_VE));
    end
    VE_mean = zeros(duration*fs,length(fre));
    for fb = 1:length(highpass)
        VE_fb_trials = [];
        % Chop data
        for i = 1:Ntrials
            VE_fb_trials = cat(1,VE_fb_trials,VE_fb(ind(i)+1:ind(i)+(duration*fs),fb));
        end
        VE_filt = reshape(VE_fb_trials,duration*fs,Ntrials);
        % Average across trials
        VE_mean(:,fb) = mean(VE_filt,2);
    end
    meanrest = mean(VE_mean(conwin(1):conwin(2),:),1);
    meanrestmat = repmat(meanrest,size(VE_mean,1),1);
    TFS = (VE_mean'-meanrestmat')./meanrestmat';
    TFS_all_mat_vis(:,:,sub_i) = TFS;

    % Filter data within bands and calculate envelope CONTROL
    VE_fb = zeros(length(VEs_cont{sub_i}),length(fre));
    for fb = 1:length(highpass)
        filt_VE = nut_filter3(VEs_cont{sub_i}','butter','bp',4,highpass(fb),lowpass(fb),fs,1)';
        VE_fb(:,fb) = abs(hilbert(filt_VE));
    end
    VE_mean = zeros(duration*fs,length(fre));
    for fb = 1:length(highpass)
        VE_fb_trials = [];
        % Chop data
        for i = 1:Ntrials
            VE_fb_trials = cat(1,VE_fb_trials,VE_fb(ind(i)+1:ind(i)+(duration*fs),fb));
        end
        VE_filt = reshape(VE_fb_trials,duration*fs,Ntrials);
        % Average across trials
        VE_mean(:,fb) = mean(VE_filt,2);
    end
    meanrest = mean(VE_mean(conwin(1):conwin(2),:),1);
    meanrestmat = repmat(meanrest,size(VE_mean,1),1);
    TFS = (VE_mean'-meanrestmat')./meanrestmat';
    TFS_all_mat_cont(:,:,sub_i) = TFS;
end
inds_sickkids = [inds_sa,inds_sk];

figure
subplot(1,2,1)
pcolor(trial_time,fre,mean(TFS_all_mat_vis(:,:,[1:76 78:end]),3));shading interp
xlabel('Time (s)');ylabel('Frequency (Hz)')
a = colorbar;%caxis([-cbar_lim cbar_lim])
a.Label.String = 'Relative Change'
axis fill
xlim([-0.8 1.8])
set(gca, 'FontSize',12)
axis square
caxis([-0.15 0.15])
subplot(1,2,2)
pcolor(trial_time,fre,mean(TFS_all_mat_cont(:,:,[1:76 78:end]),3));shading interp
xlabel('Time (s)');ylabel('Frequency (Hz)')
a = colorbar;%caxis([-cbar_lim cbar_lim])
a.Label.String = 'Relative Change'
axis fill
xlim([-0.8 1.8])
set(gca, 'FontSize',12)
axis square
caxis([-0.15 0.15])


%%
addpath C:\Users\ppynr2\OneDrive - The University of Nottingham\Documents\GitHub\Braille_kids\BrainPlots
results_dir = 'C:\Users\ppynr2\OneDrive - The University of Nottingham\Documents\GitHub\Braille_kids\';
cd(results_dir)
addpath('.\BrainPlots\')
addpath('.\gifti-1.8\')
figure
sig_regs = nan(78,1);sig_regs(65)=1;
% PaintBrodmannAreas_chooseview(1:78, 78, 256, [0,79], [], [], []);
PaintBrodmannAreas_chooseview(sig_regs, 78, 256, [0,2], [], [], [])
figure
sig_regs = nan(78,1);sig_regs(14)=1;
% PaintBrodmannAreas_chooseview(1:78, 78, 256, [0,79], [], [], []);
PaintBrodmannAreas_chooseview(sig_regs, 78, 256, [0,2], [], [], [])
