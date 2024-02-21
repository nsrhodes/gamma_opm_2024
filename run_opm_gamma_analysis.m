%% Housekeeping

clear all
close all
clc

%% select analysis: 0 = data quality, 1 = 4mm tstat, 2 = visual mask, 3 = alpha, 4 = vis and control, 5 = convert tstat units

opt_analysis = 5; 

%% select dataset: 0 = UoN kids, 1 = UoN adults, 2 = SK kids, 3 = SK adults
dataset = 3;

if dataset == 0
    load('C:\Users\ppynr2\OneDrive - The University of Nottingham\phd\Gamma\demographics_alldata.mat')
    script_dir = 'C:\Users\ppynr2\OneDrive - The University of Nottingham\Documents\GitHub\OPM_gamma_UoN_SK\';
    addpath('R:\DRS-KidsOPM\Paediatric_OPM_Notts\fieldtrip-20220906')
    project_dir = 'R:\DRS-KidsOPM\Paediatric_OPM_Notts\';
    power_line = 50; % 50 Hz powerline in UoN
    subs = [1:4 6:27];
elseif dataset == 1
    load('C:\Users\ppynr2\OneDrive - The University of Nottingham\phd\Gamma\demographics_alldata.mat')
    script_dir = 'C:\Users\ppynr2\OneDrive - The University of Nottingham\Documents\GitHub\OPM_gamma_UoN_SK\';
    addpath('R:\DRS-KidsOPM\Paediatric_OPM_Notts\fieldtrip-20220906')
     project_dir = 'R:\DRS-KidsOPM\Paediatric_OPM_Notts_AdultData\';
    power_line = 50; % 50 Hz powerline in UoN
    subs = [1:26];
elseif dataset == 2
    load('/d/mjt/9/users/natalierhodes/Gamma_Kids_SK/demographics_alldata.mat')
    script_dir = '/d/mjt/9/users/natalierhodes/OPM_gamma_UoN_SK/';
    addpath('/d/mjt/s4/toolboxes/fieldtrip/fieldtrip-20220214')
    project_dir = '/d/mjt/9/users/natalierhodes/Gamma_Kids_SK/';
    power_line = 60; % 60 Hz powerline in SK
    subs = [1:24];
elseif dataset == 3
    load('/d/mjt/9/users/natalierhodes/Gamma_Kids_SK/demographics_alldata.mat')
    script_dir = '/d/mjt/9/users/natalierhodes/OPM_gamma_UoN_SK/';
    addpath('/d/mjt/s4/toolboxes/fieldtrip/fieldtrip-20220214')
    project_dir = '/d/mjt/9/projects/OPM/SKvsNott_OPM_Comps/';
    power_line = 60; % 60 Hz powerline in SK
    subs = [1:26];
end
ft_defaults

%% run analysis

addpath(script_dir)
addpath([script_dir 'functions'])
for ss = subs
    s1 = sprintf('%2d',ss);s1(s1 == ' ') = '0'
    sub = strcat(num2str(dataset),s1);
    ses_i = 1;
    ses = sprintf('%3d',ses_i);ses(ses == ' ') = '0';
    age = demographics.age(demographics.subject==sub);
    if opt_analysis == 0
        data_quality_func(sub,ses,project_dir,power_line,age);
    elseif opt_analysis == 1
        opm_gamma_4mm(sub,ses,project_dir,power_line,age);
    elseif opt_analysis == 2
        opm_gamma_visualmask(sub,ses,project_dir,power_line,age);
    elseif opt_analysis == 3
        opm_gamma_4mm_alpha(sub,ses,project_dir,power_line,age);
    elseif opt_analysis == 4
        opm_gamma_vis_control(sub,ses,project_dir,power_line,age);
    elseif opt_analysis == 5
        convert_tstat_units(sub,ses,project_dir,'alpha') %tstat_type = 'alpha' or 'circles' for gamma
    end 
end 
