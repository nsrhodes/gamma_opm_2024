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
    load([project_dir 'sub-' sub '\sub-' sub '_ses-001_task-faces_circles_run-001_TFS_vm.mat'])
    data.TFS{sub_i} = TFS;
    % Load in virtual electrode
    load([project_dir 'sub-' sub '\sub-' sub '_ses-001_task-faces_circles_run-001_VE_unchopped.mat'])
    VEs{sub_i} = VE_unchopped;
    data.VE{sub_i} = VE_unchopped;
    inds{sub_i} = ind;
    data.N_trials{sub_i} = length(ind);
    Ntrials = length(ind);
   % Ntrials = 18;
    Ntrials_all(sub_i) = Ntrials;
    % load in Nchans
%     load([project_dir 'sub-' sub '\sub-' sub '_Nchans.mat']);
%     data.Nchans{sub_i} = Nchans;
%     Nchans_all(sub_i) = Nchans;
    % Filter data within bands and calculate envelope
    VE_fb = zeros(length(VE_unchopped),length(fre));
    for fb = 1:length(highpass)
        filt_VE = nut_filter3(VE_unchopped','butter','bp',4,highpass(fb),lowpass(fb),fs,1)';
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
    TFS_all_mat(:,:,sub_i) = TFS;
    % Get reduced channel count nottingham adult data
    if sub(1)=='1'
        count = count+1;
        load([project_dir 'sub-' sub '\sub-' sub '_ses-001_task-faces_circles_run-001_TFS_SK.mat'])
        TFS_na_skchans(:,:,count) = TFS;
        load([project_dir 'sub-' sub '\sub-' sub '_ses-001_task-faces_circles_run-001_VE_unchopped_SK.mat'])
        VE_na_skchans{count} = VE_unchopped;
        inds_na_skchans{count} = ind;
        % Filter data within bands and calculate envelope
        VE_fb = zeros(length(VE_unchopped),length(fre));
        for fb = 1:length(highpass)
            filt_VE = nut_filter3(VE_unchopped','butter','bp',4,highpass(fb),lowpass(fb),fs,1)';
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
        TFS_all_mat_SKchans(:,:,count) = TFS;
    end
end
inds_sickkids = [inds_sa(2:end),inds_sk]-1;
inds_uon = [inds_na, inds_nk];
%% Make demographics figure

binrng = 1:35; % Create Bin Ranges
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
legend('SK Male', 'UoN Male', 'SK Female', 'UoN Female')
ylim([0 11]);ylabel('Frequency');xlabel('Age (years)')
%% Make separate demographics figures for sites
figure
subplot(1,2,1)
bar(binrng(1:end-1),histcounts(demo_sk.age,binrng));
xlim([0 35]);ylim([0 11]); axis square; grid on
xlabel('Age (years)');ylabel('Frequency')
set(gca, 'FontSize',12)

subplot(1,2,2)
bar(binrng(1:end-1),histcounts(demo_uon.age,binrng))
xlim([0 35]);ylim([0 11]); axis square; grid on
xlabel('Age (years)');ylabel('Frequency')
set(gca, 'FontSize',12)


%% Adult site comparison
cbar_lim = 0.3;
figure;
subplot(2,2,1)
pcolor(trial_time,fre,mean(TFS_all_mat_SKchans,3));shading interp
xlabel('Time (s)');ylabel('Frequency (Hz)')
a = colorbar;caxis([-cbar_lim cbar_lim])
a.Label.String = 'Relative Change'
axis fill
axis square
xlim([-0.8 1.8])

set(gca, 'FontSize',12)
subplot(2,2,2)
pcolor(trial_time,fre,mean(TFS_all_mat(:,:,inds_sa),3));shading interp
xlabel('Time (s)');ylabel('Frequency (Hz)')
a = colorbar;caxis([-cbar_lim cbar_lim])
a.Label.String = 'Relative Change'
axis fill
axis square
xlim([-0.8 1.8])
set(gca, 'FontSize',12)
subplot(2,2,3:4)
% Get envelope from 30-80 Hz

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
plot_area(trial_time, mean(mean_Env(:,inds_na),2), std(mean_Env(:,inds_na))./sqrt(26), 'r',[1 0.4 0.6])
plot_area(trial_time, mean(mean_Env(:,inds_sa),2), std(mean_Env(:,inds_sa))./sqrt(26), 'b',[0.7 0.8 1])
set(gca, 'FontSize',12)
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
ylim([-0.1 0.4])
% statistical testing on envelopes
mean_Env_na = mean_Env(:,inds_na);
mean_Env_sa = mean_Env(:,inds_sa);
mean_Env_na_skchans = mean_Env_skchans;

for i = 1:3600
    [h,p,ci,stats] = ttest(mean_Env_na_skchans(i,:),mean_Env_sa(i,:));
    tstat(i) = stats.tstat;
    pval(i) = p;
end

percent_change = mean(mean_Env(1.3*fs:2*fs,:),1)

%% ttest envelope
percent_change_NA = mean(mean_Env_na(1.3*fs:2*fs,:),1);
percent_change_SK = mean(mean_Env_sa(1.3*fs:2*fs,:),1);
percent_change_NA_SKchans = mean(mean_Env_na_skchans(1.3*fs:2*fs,:),1);

[h,p_nask] = ttest2(percent_change_NA,percent_change_SK)
[h,p_nasksk_skchans] = ttest2(percent_change_NA_SKchans,percent_change_SK)

%% Individual participant TFS
caxlims = [-0.5 0.5];
figure;
subplot(1,2,1);
pcolor(trial_time,fre,TFS_all_mat(:,:,27));shading interp
xlabel('Time (s)');ylabel('Frequency (Hz)')
a = colorbar;caxis(caxlims)
a.Label.String = 'Relative Change';
axis fill
xlim([-0.8 1.8])
set(gca, 'FontSize',12)
axis square 
subplot(1,2,2);
pcolor(trial_time,fre,TFS_all_mat(:,:,77));shading interp
xlabel('Time (s)');ylabel('Frequency (Hz)')
a = colorbar;caxis(caxlims)
a.Label.String = 'Relative Change';
axis fill
xlim([-0.8 1.8])
set(gca, 'FontSize',12)
axis square 
%% 30 - 80 Hz envelope ttest SUB1
sub1change_uon = percent_change(27)
sub1change_sk = percent_change(77)


%% Get age separated indices

addpath C:\Users\ppynr2\OneDrive - The University of Nottingham\Documents\GitHub\Gamma_circles_development\functions\
age = data.age;
age(77) = [];
[age_order, age_inds] = sort(age,'ascend');
inds_2to4 = age_inds(age_order<=4);
inds_5to8 = age_inds(age_order>4&age_order<9);
inds_9to13 = age_inds(age_order>8&age_order<14);
inds_21to24 = age_inds(age_order>20&age_order<25);
inds_25to28 = age_inds(age_order>24&age_order<29);
inds_29to34 = age_inds(age_order>28&age_order<35);

%% Age separated TFS 
TFS_all_mat(:,:,77)=[];
%%
TFS_2to4 = mean(TFS_all_mat(:,:,inds_2to4),3);
TFS_5to8 = mean(TFS_all_mat(:,:,inds_5to8),3);
TFS_9to13 = mean(TFS_all_mat(:,:,inds_9to13),3);
TFS_21to24 = mean(TFS_all_mat(:,:,inds_21to24),3);
TFS_25to28 = mean(TFS_all_mat(:,:,inds_25to28),3);
TFS_29to34 = mean(TFS_all_mat(:,:,inds_29to34),3);
caxlims = [-0.18 0.18];
figure;
subplot(2,3,1);
pcolor(trial_time,fre,TFS_2to4);shading interp
xlabel('Time (s)');ylabel('Frequency (Hz)')
a = colorbar;caxis(caxlims)
a.Label.String = 'Relative Change';
axis fill
xlim([-0.8 1.8])
set(gca, 'FontSize',12)
title('aged 2-4')
axis square 
subplot(2,3,2);
pcolor(trial_time,fre,TFS_5to8);shading interp
xlabel('Time (s)');ylabel('Frequency (Hz)')
a = colorbar;caxis(caxlims)
a.Label.String = 'Relative Change';
axis fill
xlim([-0.8 1.8])
set(gca, 'FontSize',12)
title('aged 5-9')
axis square 
subplot(2,3,3);
pcolor(trial_time,fre,TFS_9to13);shading interp
xlabel('Time (s)');ylabel('Frequency (Hz)')
a = colorbar;caxis(caxlims)
a.Label.String = 'Relative Change'
axis fill
xlim([-0.8 1.8])
set(gca, 'FontSize',12)
title('aged 10-13')
axis square 
subplot(2,3,4);
pcolor(trial_time,fre,TFS_21to24);shading interp
xlabel('Time (s)');ylabel('Frequency (Hz)')
a = colorbar;caxis(caxlims)
a.Label.String = 'Relative Change'
axis fill
xlim([-0.8 1.8])
set(gca, 'FontSize',12)
title('aged 21-24')
axis square 
subplot(2,3,5);
pcolor(trial_time,fre,TFS_25to28);shading interp
xlabel('Time (s)');ylabel('Frequency (Hz)')
a = colorbar;caxis(caxlims)
a.Label.String = 'Relative Change'
axis fill
xlim([-0.8 1.8])
set(gca, 'FontSize',12)
title('aged 25-28')
axis square 
subplot(2,3,6);
pcolor(trial_time,fre,TFS_29to34);shading interp
xlabel('Time (s)');ylabel('Frequency (Hz)')
a = colorbar;caxis(caxlims)
a.Label.String = 'Relative Change'
axis fill
xlim([-0.8 1.8])
set(gca, 'FontSize',12)
title('aged 29-34')
axis square 
%% Plot spectra
fre = highpass + ((lowpass - highpass)./2);

psd_norm = squeeze(mean(TFS_all_mat(2:end,1.3*fs:2*fs,:),2))';
fre = fre(2:end);
cols = lines(6);
figure();
plot_area_colours(fre, mean(psd_norm(inds_2to4,:),1)', (std(psd_norm(inds_2to4,:))./sqrt(length(inds_2to4)))', cols(1,:),cols(1,:))
%set(gca,'Yscale',['log']);
xlim([5 79])
set(gca,'FontSize',12)
xlabel('Frequency (Hz)'); ylabel('Relative change from baseline (\Delta)'); grid on; hold on
plot_area_colours(fre, mean(psd_norm(inds_5to8,:),1)', (std(psd_norm(inds_5to8,:))./sqrt(length(inds_5to8)))',  cols(2,:),cols(2,:))
plot_area_colours(fre, mean(psd_norm(inds_9to13,:),1)', (std(psd_norm(inds_9to13,:))./sqrt(length(inds_9to13)))',  cols(3,:),cols(3,:))
plot_area_colours(fre, mean(psd_norm(inds_21to24,:),1)', (std(psd_norm(inds_21to24,:))./sqrt(length(inds_21to24)))',  cols(4,:),cols(4,:))
plot_area_colours(fre, mean(psd_norm(inds_25to28,:),1)', (std(psd_norm(inds_25to28,:))./sqrt(length(inds_25to28)))',  cols(5,:),cols(5,:))
plot_area_colours(fre, mean(psd_norm(inds_29to34,:),1)', (std(psd_norm(inds_29to34,:))./sqrt(length(inds_29to34)))',  cols(6,:),cols(6,:))
xline(13,'LineWidth',2);
xline(31,'LineWidth',2);
xline(53,'LineWidth',2);
yline(0,'LineWidth',1);
legend('',['2-4 years (n=' num2str(length(inds_2to4)) ')'],'',['5-8 years (n=' num2str(length(inds_5to8)) ')'],'',['9-13 years (n=' num2str(length(inds_9to13)) ')'],'',['21-24 years (n=' num2str(length(inds_21to24)) ')'],'',['25-28 years (n=' num2str(length(inds_25to28)) ')'],'',['29-34 years (n=' num2str(length(inds_29to34)) ')'],'','','','')
grid off; 

%% Single sub spectra
figure;
plot(fre, psd_norm(1,:));
hold on
plot(fre, psd_norm(77,:));
%% Ttest all freqs of spectra

for i = 1:length(fre)
    [r_all(i),p_all(i)] = corr(age,abs(psd_norm(:,i)));
end
figure
plot(fre,r_all);%hold on; plot(fre,p_all);
xlabel('Centre frequency (Hz)'); ylabel('R');
figure
plot(fre,p_all);hold on; yline(0.05,'g'); yline(0.05./44,'r')
xlabel('Centre frequency (Hz)'); ylabel('p');ylim([-0.1 0.9])
legend('p','p=0.05 (uncorrected)','p=0011 (corrected)')

%% figures of selected frequency slopes

f1 = figure;
fwidth = 300;
fheight = 300;
f1.Position([3,4]) = [fwidth,fheight];
plot(age(inds_sickkids),psd_norm(inds_sickkids,fre==13),'bo')
hold on
plot(age(inds_uon),psd_norm(inds_uon,fre==13),'ro')
xlabel('Age (years)');ylabel('\Delta')
[r13,p13]=corr(age,psd_norm(:,fre==13))
X = [ones(length(age),1) age];
y = psd_norm(:,fre==13);
b = X\y;
refl = refline(b(2),b(1));
axis square
box off
xlim([0 38])
set(gca,'FontSize',12)

f2 = figure;
f2.Position([3,4]) = [fwidth,fheight];
plot(age,psd_norm(:,fre==31),'ko')
xlabel('Age (years)');ylabel('\Delta')
[r31,p31]=corr(age,psd_norm(:,fre==31))
y = psd_norm(:,fre==31);
b = X\y;
refl = refline(b(2),b(1));
axis square
box off
xlim([0 38])
set(gca,'FontSize',12)

f3 = figure;
f3.Position([3,4]) = [fwidth,fheight];
plot(age,psd_norm(:,fre==53),'ko')
xlabel('Age (years)');ylabel('\Delta')
[r53,p53]=corr(age,psd_norm(:,fre==53))
y = psd_norm(:,fre==53);
b = X\y;
refl = refline(b(2),b(1));
axis square
box off
xlim([0 38])
set(gca,'FontSize',12)

%% Get peak frequencies
psd_gamma = psd_norm(:,14:38);
for i = 1:size(psd_gamma,1)
[val,index]=max(psd_gamma(i,:));
peak_gamma_freq(i)=fre(index+13);
peak_gamma_amp(i)=val;
end 
fwidth = 350;
fheight = 300;
%gamma_inds = find(peak_gamma_amp>0.0);
f4 = figure;
f4.Position([3,4]) = [fwidth,fheight];
plot(age,peak_gamma_freq,'ko')
xlabel('Age (years)');ylabel('Peak gamma frequency (Hz)')
[rpeak,ppeak]=corr(age,peak_gamma_freq')
set(gca,'FontSize',12)
y = peak_gamma_freq';
b = X\y;
refl = refline(b(2),b(1));
box off; xlim([0 35])
f5 = figure;
f5.Position([3,4]) = [fwidth,fheight];
plot(age,peak_gamma_amp,'ko')
xlabel('Age (years)');ylabel('Amplitude (A.U.)')
[ramp,pamp]=corr(age,peak_gamma_amp')
set(gca,'FontSize',12)
y = peak_gamma_amp';
b = X\y;
refl = refline(b(2),b(1));
box off; xlim([0 35])
%% examples 
%good: 13, 37, 53
%bad: 50, 65, 101
f6 = figure;
subplot(3,1,1)
plot(fre(14:38),psd_gamma(13,:),'LineWidth',2,'Color',cols(1,:))
xlabel('Frequency (Hz)');ylabel('Amplitude (A.U.)')
set(gca,'FontSize',12)
axis([30 80 -0.1 1.2])
axis square
subplot(3,1,2)
plot(fre(14:38),psd_gamma(37,:),'LineWidth',2,'Color',cols(2,:))
xlabel('Frequency (Hz)');ylabel('Amplitude (A.U.)')
set(gca,'FontSize',12)
axis([30 80 -0.1 1.2])
axis square
subplot(3,1,3)
plot(fre(14:38),psd_gamma(53,:),'LineWidth',2,'Color',cols(3,:))
xlabel('Frequency (Hz)');ylabel('Amplitude (A.U.)')
set(gca,'FontSize',12)
axis([30 80 -0.1 1.2])
axis square

f7 = figure;
subplot(3,1,1)
plot(fre(14:38),psd_gamma(50,:),'LineWidth',2,'Color',cols(4,:))
xlabel('Frequency (Hz)');ylabel('Amplitude (A.U.)')
set(gca,'FontSize',12)
axis([30 80 -0.1 1.2])
axis square
subplot(3,1,2)
plot(fre(14:38),psd_gamma(65,:),'LineWidth',2,'Color',cols(5,:))
xlabel('Frequency (Hz)');ylabel('Amplitude (A.U.)')
set(gca,'FontSize',12)
axis([30 80 -0.1 1.2])
axis square
subplot(3,1,3)
plot(fre(14:38),psd_gamma(101,:),'LineWidth',2,'Color',cols(6,:))
xlabel('Frequency (Hz)');ylabel('Amplitude (A.U.)')
set(gca,'FontSize',12)
axis([30 80 -0.1 1.2])
axis square

%% Channel and trial count

% UoN_Nchans_mean = mean(Nchans_all(1:52));
% UoN_Nchans_sd = std(Nchans_all(1:52));
% SK_Nchans_mean = mean(Nchans_all(53:end));
% SK_Nchans_sd = std(Nchans_all(53:end));

Adults_Ntrials_mean = mean(Ntrials_all([inds_na inds_sa]));
Adults_Ntrials_sd = std(Ntrials_all([inds_na inds_sa]));
Kids_Ntrials_mean = mean(Ntrials_all([inds_nk inds_sk]));
Kids_Ntrials_sd = std(Ntrials_all([inds_nk inds_sk]));
