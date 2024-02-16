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
inds_sickkids = [inds_sa,inds_sk];


% Define TFS parameters
highpass = [1 2 4 6 8 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100 105 110];
lowpass = [4 6 8 10 13 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100 105 110 115 120];
fre = highpass + ((lowpass - highpass)./2);
trig_offset = -1;
fs = 1200;
trial_time = linspace(0+trig_offset,size(TFS,2)./fs+trig_offset,size(TFS,2));
cbar_lim = 0.3;



%% Get VEs

fs = 1200;
hp = 30; lp = 80;
[b,a] = butter(4,2*[hp lp]/fs);
VE_on_all = [];
age = data.age;
age(1) = [];
[age_order, age_inds] = sort(age,'ascend');
inds_2to4 = age_inds(1:23);
inds_5to8 = age_inds(24:37);
inds_9to13 = age_inds(38:49);
inds_21to24 = age_inds(50:69);
inds_25to28 = age_inds(70:87);
inds_29to34 = age_inds(88:end);

for sub_i = 1:size(demographics,1)
    sub = demographics.subject{sub_i};
    N_trials = data.N_trials{sub_i};
    VE = data.VE{sub_i};
    data_f = [filtfilt(b,a,VE')];
    VE_mat = reshape(data_f,[],N_trials);
    VE_on_mat = VE_mat(1.3*fs:2*fs,:);
    VE_on = VE_on_mat(:)';
    VE_on_all = [VE_on_all VE_on];
    [orig, osc, frac] = do_fooof(VE);
    psd(sub_i,:) = orig.powspctrm;
    psd_fooof_osc(sub_i,:) = osc.powspctrm;
    psd_fooof_frac(sub_i,:) = frac.powspctrm;
end

psd_norm = psd;%./mean(psd')';
% display the spectra on a log-log scale
addpath C:\Users\ppynr2\OneDrive - The University of Nottingham\Documents\GitHub\Gamma_circles_development\functions\
figure();
cols = lines(6);
plot_area_colours(orig.freq, mean(psd_norm(inds_2to4,:),1)', (std(psd_norm(inds_2to4,:))./sqrt(length(inds_2to4)))', cols(1,:),cols(1,:))
set(gca,'Yscale',['log']); xlim([0 100])
xlabel('Frequency'); ylabel('Power'); grid on; hold on
plot_area_colours(orig.freq, mean(psd_norm(inds_5to8,:),1)', (std(psd_norm(inds_5to8,:))./sqrt(length(inds_5to8)))',  cols(2,:),cols(2,:))
plot_area_colours(orig.freq, mean(psd_norm(inds_9to13,:),1)', (std(psd_norm(inds_9to13,:))./sqrt(length(inds_9to13)))',  cols(3,:),cols(3,:))
plot_area_colours(orig.freq, mean(psd_norm(inds_21to24,:),1)', (std(psd_norm(inds_21to24,:))./sqrt(length(inds_21to24)))',  cols(4,:),cols(4,:))
plot_area_colours(orig.freq, mean(psd_norm(inds_25to28,:),1)', (std(psd_norm(inds_25to28,:))./sqrt(length(inds_25to28)))',  cols(5,:),cols(5,:))
plot_area_colours(orig.freq, mean(psd_norm(inds_29to34,:),1)', (std(psd_norm(inds_29to34,:))./sqrt(length(inds_29to34)))',  cols(6,:),cols(6,:))

[TFS] = VE_TFS(VE,[0.2 1],0,3,ind,N_trials,trial_time,fs);

%% wavelet transform 
clear TFS;
[wt,wf] = cwt(VE','amor',fs); % morlet wavelet

centf = 3:2:100
half_width = 2.5
for fi = 1:length(centf)
    lo = centf(fi)-half_width;
    hi = centf(fi)+half_width;
    band_freqs = find(wf > lo & wf < hi);
     % Envelope of oscillations
            TFS(fi,:) = mean(abs(wt(band_freqs,:)),1)';
end
TFS_avg = mean(reshape(TFS,length(centf),length(trial_time),N_trials),3,'omitnan');
TFS_change = TFS_avg - mean(TFS_avg(:,0.2*fs:0.9*fs),2);
figure
pcolor(trial_time,centf,TFS_change);shading interp
xlabel('Time (s)');ylabel('Frequency (Hz)')
a = colorbar;%caxis([-cbar_lim cbar_lim])
a.Label.String = 'Relative Change'
axis fill
xlim([-0.8 1.8])
set(gca, 'FontSize',12)
%% Fooof

function [orig, osc, frac] = do_fooof(VE)
fooofdata.trial{1,1} = VE;
fooofdata.time{1,1} = 1/1200:1/1200:length(VE)/1200;
fooofdata.label{1}     = 'chan';
fooofdata.trialinfo(1,1) = 1;
fooofdata.type = 'raw';
% chunk into 3-second segments
cfg               = [];
cfg.length        = 2;
cfg.overlap       = 0.5;
fooofdata              = ft_redefinetrial(cfg, fooofdata);

% compute the fractal and original spectra
cfg               = [];
cfg.foilim        = [1 200];
cfg.pad           = 4;
cfg.tapsmofrq     = 2;
cfg.method        = 'mtmfft';
cfg.output        = 'fooof_aperiodic';
frac = ft_freqanalysis(cfg, fooofdata);
cfg.output        = 'pow';
orig = ft_freqanalysis(cfg, fooofdata);

% subtract the fractal component from the power spectrum
cfg               = [];
cfg.parameter     = 'powspctrm';
cfg.operation     = 'x2-x1';
osc = ft_math(cfg, frac, orig);
end