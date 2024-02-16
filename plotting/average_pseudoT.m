clear all
close all
clc
%% Fieldtrip set up 

addpath('R:\DRS-KidsOPM\Paediatric_OPM_Notts\fieldtrip-20220906')
ft_defaults

%% Select participants 

base_dir = 'R:\DRS-KidsOPM\';
load('C:\Users\ppynr2\OneDrive - The University of Nottingham\phd\Gamma\demographics_alldata.mat')
proj_dir_nk = [base_dir 'Paediatric_OPM_Notts\Data\BIDS\derivatives\Tstats\'];
proj_dir_na = [base_dir 'Paediatric_OPM_Notts_AdultData\Data\BIDS\derivatives\Tstats\'];
proj_dir_sk = [base_dir 'SickKids\Tstats\'];


%% File set-up
script_dir = 'C:\Users\ppynr2\OneDrive - The University of Nottingham\Documents\GitHub\Gamma_circles_development\';
addpath(script_dir)
addpath([script_dir 'functions'])

%% Load MNI mesh

load('C:\Users\ppynr2\OneDrive - The University of Nottingham\Documents\GitHub\Gamma_circles_development\plotting\mni_meshes.mat');
 % plot meshes
 figure(1);ft_plot_mesh(meshes,'facecolor',[.5 .5 .5],'facealpha',.3,'edgecolor','none')
 hold on
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
    % Load in pseudoT stat
    Tstat = ft_read_mri([project_dir 'sub-' sub '\sub-' sub '_ses-001_task-faces_circles_run-001_pseudoT_circles_mm_MNI_4mm.nii.gz']);
    all_Tstats(:,:,:,sub_i) = Tstat.anatomy;
  
end 

%% Get UoN SK adult comparison
mean_Tstat_na = Tstat; mean_Tstat_na.anatomy=mean(all_Tstats(:,:,:,inds_na),4);
mean_Tstat_sa = Tstat; mean_Tstat_sa.anatomy=mean(all_Tstats(:,:,:,inds_sa),4);

ft_write_mri('mean_Tstat_na.nii',mean_Tstat_na,'dataformat','nifti');
ft_write_mri('mean_Tstat_sa.nii',mean_Tstat_sa,'dataformat','nifti');

%% Get age_separated
age = data.age;
age(77) = [];
[age_order, age_inds] = sort(age,'ascend');
inds_2to4 = age_inds(age_order<=4);
inds_5to8 = age_inds(age_order>4&age_order<9);
inds_9to13 = age_inds(age_order>8&age_order<14);
inds_21to24 = age_inds(age_order>20&age_order<25);
inds_25to28 = age_inds(age_order>24&age_order<29);
inds_29to34 = age_inds(age_order>28&age_order<35);

mean_Tstat_2to4 = Tstat; mean_Tstat_2to4.anatomy=mean(all_Tstats(:,:,:,inds_2to4),4);
mean_Tstat_5to8 = Tstat; mean_Tstat_5to8.anatomy=mean(all_Tstats(:,:,:,inds_5to8),4);
mean_Tstat_9to13 = Tstat; mean_Tstat_9to13.anatomy=mean(all_Tstats(:,:,:,inds_9to13),4);
mean_Tstat_21to24 = Tstat; mean_Tstat_21to24.anatomy=mean(all_Tstats(:,:,:,inds_21to24),4);
mean_Tstat_25to28 = Tstat; mean_Tstat_25to28.anatomy=mean(all_Tstats(:,:,:,inds_25to28),4);
mean_Tstat_29to34 = Tstat; mean_Tstat_29to34.anatomy=mean(all_Tstats(:,:,:,inds_29to34),4);

mean_Tstat = Tstat;
mean_Tstat.anatomy = mean(all_Tstats,4);
%ft_write_mri('mean_Tstat.nii',mean_Tstat,'dataformat','nifti');

ft_write_mri('mean_Tstat_2to4.nii',mean_Tstat_2to4,'dataformat','nifti');
ft_write_mri('mean_Tstat_5to8.nii',mean_Tstat_5to8,'dataformat','nifti');
ft_write_mri('mean_Tstat_9to13.nii',mean_Tstat_9to13,'dataformat','nifti');
ft_write_mri('mean_Tstat_21to24.nii',mean_Tstat_21to24,'dataformat','nifti');
ft_write_mri('mean_Tstat_25to28.nii',mean_Tstat_25to28,'dataformat','nifti');
ft_write_mri('mean_Tstat_29to34.nii',mean_Tstat_29to34,'dataformat','nifti');
