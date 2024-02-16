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

for sub_i = 1:size(demographics,1)
    sub = demographics.subject{sub_i};
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
    load([project_dir 'sub-' sub '\sub-' sub '_ses-001_task-faces_circles_run-001_TFS_vm.mat'])
    data.TFS{sub_i} = TFS;
    TFS_mat(:,:,sub_i) = TFS;
    % Load in virtual electrode
    load([project_dir 'sub-' sub '\sub-' sub '_ses-001_task-faces_circles_run-001_VE_unchopped.mat'])
    VEs{sub_i} = VE_unchopped;
    data.VE{sub_i} = VE_unchopped;
    inds{sub_i} = ind;
    data.N_trials{sub_i} = length(ind);
    Ntrials(sub_i) = length(ind);
    % Get reduced channel count nottingham adult data
    if sub(1)=='1'
        count = count+1;
        load([project_dir 'sub-' sub '\sub-' sub '_ses-001_task-faces_circles_run-001_TFS_SK.mat'])
        TFS_na_skchans(:,:,count) = TFS;
        load([project_dir 'sub-' sub '\sub-' sub '_ses-001_task-faces_circles_run-001_VE_unchopped_SK.mat'])
        VE_na_skchans{count} = VE_unchopped;
        inds_na_skchans{count} = ind;
    end 
end 

%% Demographics figure

binrng = 1:35; % Create Bin Ranges
inds_sickkids = [inds_sa,inds_sk];
demo_sk = demographics(inds_sickkids,:);
demo_uon = demographics([inds_na, inds_nk],:);
counts1 = histcounts(demo_sk.age(demo_sk.sex=='M'), binrng); % Histogram For male
counts2 = histcounts(demo_sk.age(demo_sk.sex=='F'), binrng); % Histogram For female
counts3 = counts1 + counts2;  % Histogram Sum

counts4 = histcounts(demo_uon.age(demo_uon.sex=='M'), binrng);
counts5 = histcounts(demo_uon.age(demo_uon.sex=='F'), binrng); % Histogram For female
counts6 = counts4 + counts5;  % Histogram Sum

% plot
figure(1)
all_m = [counts3 ; counts6];
all_f = [counts1 ; counts4];
bar(binrng(1:end-1), all_m)
hold on
bar(binrng(1:end-1), all_f)
hold off
legend('SK Male', 'Notts Male', 'SK Female', 'Notts Female')
ylim([0 11]);ylabel('Frequency');xlabel('Age (years)')

%% Organise TFS from data table

TFS = cell2mat(data.TFS);
TFS = reshape(TFS,26,[],3600); TFS = permute(TFS,[1 3 2]);
TFS_na = TFS(:,:,inds_na);TFS_nk = TFS(:,:,inds_nk); TFS_sk = TFS(:,:,inds_sk); TFS_sa = TFS(:,:,inds_sa);

%% Plot TFS

% Define TFS parameters
highpass = [1 2 4 6 8 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100 105 110];
lowpass = [4 6 8 10 13 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100 105 110 115 120];
fre = highpass + ((lowpass - highpass)./2);
trig_offset = -1;
fs = 1200;
trial_time = linspace(0+trig_offset,size(TFS,2)./fs+trig_offset,size(TFS,2));
cbar_lim = 0.3;

% Mean TFS
TFS_mean_na = mean(TFS_na,3);
TFS_mean_nk = mean(TFS_nk,3);
TFS_mean_sa = mean(TFS_sa,3);
TFS_mean_sk = mean(TFS_sk,3);
TFS_mean_na_skchans = mean(TFS_na_skchans,3);

figure(2)
subplot(2,2,1)
pcolor(trial_time,fre,TFS_mean_na_skchans);shading interp
xlabel('Time (s)');ylabel('Frequency (Hz)')
a = colorbar;caxis([-cbar_lim cbar_lim])
a.Label.String = 'Relative Change'
axis fill
xlim([-0.8 1.8])
set(gca, 'FontSize',12)
title('Nottingham Adults')

subplot(2,2,2)
pcolor(trial_time,fre,TFS_mean_sa);shading interp
xlabel('Time (s)');ylabel('Frequency (Hz)')
a = colorbar;caxis([-cbar_lim cbar_lim])
a.Label.String = 'Relative Change'
axis fill
xlim([-0.8 1.8])
set(gca, 'FontSize',12)
title('SickKids Adults')

subplot(2,2,3)
pcolor(trial_time,fre,TFS_mean_nk);shading interp
xlabel('Time (s)');ylabel('Frequency (Hz)')
a = colorbar;caxis([-cbar_lim cbar_lim])
a.Label.String = 'Relative Change'
axis fill
xlim([-0.8 1.8])
set(gca, 'FontSize',12)
title('Nottingham Kids')

subplot(2,2,4)
pcolor(trial_time,fre,TFS_mean_sk);shading interp
xlabel('Time (s)');ylabel('Frequency (Hz)')
a = colorbar;caxis([-cbar_lim cbar_lim])
a.Label.String = 'Relative Change'
axis fill
xlim([-0.8 1.8])
set(gca, 'FontSize',12)
title('SickKids Kids')

%% Get envelope from 30-80 Hz

control_inds = 0.2*1200:0.9*1200;
fs = 1200;
hp = 30; lp = 80;
[b,a] = butter(4,2*[hp lp]/fs);
for sub_i = 1:size(demographics,1)
    sub = demographics.subject{sub_i};
    N_trials = data.N_trials{sub_i};
    VE = data.VE{sub_i};
    
    data_f = [filtfilt(b,a,VE')];
    VE_mat = reshape(data_f,[],N_trials);
    for tr_i = 1:N_trials
    Envs(:,tr_i) = abs(hilbert(VE_mat(:,tr_i)));
    Envs(:,tr_i) = (Envs(:,tr_i) - mean(Envs(control_inds,tr_i)))./mean(Envs(control_inds,tr_i));
    Envs_on = mean(Envs(1.3*fs:2*fs,:),1);
    Envs_off = mean(Envs(control_inds,:),1);
    end
    mean_Env(:,sub_i) = mean(Envs,2);
    [h(sub_i),stat_p(sub_i)] = ttest(Envs_on,Envs_off);
end

figure;
plot_area(trial_time, mean(mean_Env(:,inds_na),2), std(mean_Env(:,inds_na))./sqrt(26), 'r',[1 0.4 0.6])
plot_area(trial_time, mean(mean_Env(:,inds_sa),2), std(mean_Env(:,inds_sa))./sqrt(26), 'b',[0.7 0.8 1])

% plot(trial_time,mean(mean_Env(:,inds_na),2)); hold on;
% plot(trial_time,mean(mean_Env(:,inds_sa),2));

for sub_skchans = 1:size(VE_na_skchans,2)
    N_trials = length(inds_na_skchans{sub_skchans});
    VE = VE_na_skchans{sub_skchans};
    data_f = [filtfilt(b,a,VE')];
    VE_mat = reshape(data_f,[],N_trials);
    for tr_i = 1:N_trials
    Envs(:,tr_i) = abs(hilbert(VE_mat(:,tr_i)));
    Envs(:,tr_i) = (Envs(:,tr_i) - mean(Envs(control_inds,tr_i)))./mean(Envs(control_inds,tr_i));
    end
    mean_Env_skchans(:,sub_skchans) = mean(Envs,2);
end

plot_area(trial_time, mean(mean_Env_skchans,2), std(mean_Env_skchans)./sqrt(26),'y',[0.9290 0.6940 0.1250])
%plot(trial_time,mean(mean_Env_skchans,2));
xlim([-0.8 1.8]);xlabel('Time (s)');ylabel('Amplitude')

% statistical testing on envelopes
mean_Env_na = mean_Env(:,inds_na);
mean_Env_sa = mean_Env(:,inds_sa);
mean_Env_na_skchans = mean_Env_skchans;

for i = 1:3600
    [h,p,ci,stats] = ttest(mean_Env_na_skchans(i,:),mean_Env_sa(i,:));
    tstat(i) = stats.tstat;
    pval(i) = p;
end

percent_change = mean(mean_Env(1.1*fs:1.9*fs,:),1)

%% Plot TFS from single subject comparison

figure;
subplot(1,2,1)
pcolor(trial_time,fre,TFS(:,:,27));shading interp
xlabel('Time (s)');ylabel('Frequency (Hz)')
a = colorbar;caxis([-cbar_lim cbar_lim])
a.Label.String = 'Relative Change'
axis fill
xlim([-0.8 1.8])
set(gca, 'FontSize',12)
title('Nottingham')
subplot(1,2,2)
pcolor(trial_time,fre,TFS(:,:,77));shading interp
xlabel('Time (s)');ylabel('Frequency (Hz)')
a = colorbar;caxis([-cbar_lim cbar_lim])
a.Label.String = 'Relative Change'
axis fill
xlim([-0.8 1.8])
set(gca, 'FontSize',12)
title('SickKids')
