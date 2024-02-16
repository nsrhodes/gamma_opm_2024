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
colors = ([255:-7:1;zeros(1,37);1:7:255])./255;
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
    % Load in peak in MNI
    fid = fileread([project_dir, 'sub-' sub '\circles_peak_vm_mni.txt']);
    dip_loc_mm = str2num(fid(43:end));
    dip_loc = dip_loc_mm./1000;
    
    figure(1)
    plot3(dip_loc(1),dip_loc(2),dip_loc(3),'.','markerSize',20, 'Color',colors(:,demographics.age(sub_i))')
    all_dip_locs(sub_i,:) = dip_loc;
%     load([project_dir,'sub-',sub '\dip_loc_vm.mat'])
%     dip_loc_circles_mm = dip_loc_circles.*1000;
%     save([project_dir,'sub-', sub,'\circles_peak_vm.txt'],'dip_loc_circles_mm','-ascii','-double')
     if sub(1)=='1'
         count = count+1;
        fid = fileread([project_dir, 'sub-' sub '\circles_peak_vm__SKchans_mni.txt']);
         dip_loc_mm = str2num(fid(43:end));
         dip_loc = dip_loc_mm./1000;
         all_dip_locs_SKchans(count,:) = dip_loc;
    end
end

%% Plot individual peaks
figure(3)
ft_plot_mesh(meshes,'facecolor',[.5 .5 .5],'facealpha',.3,'edgecolor','none')
hold on
plot3(all_dip_locs(1,1),all_dip_locs(1,2),all_dip_locs(1,3),'r.','MarkerSize',30);
plot3(all_dip_locs(77,1),all_dip_locs(77,2),all_dip_locs(77,3),'b.','MarkerSize',30);

dist_indiv = sqrt((all_dip_locs(1,1)-all_dip_locs(77,1))^2+(all_dip_locs(1,2)-all_dip_locs(77,2))^2+(all_dip_locs(1,3)-all_dip_locs(77,3))^2);
%% Plot system comparison ellipsoid of peak locations overlaid
figure(4)
ft_plot_mesh(meshes,'facecolor',[.5 .5 .5],'facealpha',.3,'edgecolor','none')
hold on

uon_dip_locs = [mean(all_dip_locs(inds_na,:)); std(all_dip_locs(inds_na,:))];
[uon_x,uon_y,uon_z] = ellipsoid(uon_dip_locs(1,1),uon_dip_locs(1,2),uon_dip_locs(1,3),...
    uon_dip_locs(2,1),uon_dip_locs(2,2),uon_dip_locs(2,3));
uon_dip_locs_skchans = [mean(all_dip_locs_SKchans); std(all_dip_locs_SKchans)];
[uonskchans_x,uonskchans_y,uonskchans_z] = ellipsoid(uon_dip_locs_skchans(1,1),uon_dip_locs_skchans(1,2),uon_dip_locs_skchans(1,3),...
    uon_dip_locs_skchans(2,1),uon_dip_locs_skchans(2,2),uon_dip_locs_skchans(2,3));
sk_dip_locs = [mean(all_dip_locs(inds_sa,:)); std(all_dip_locs(inds_sa,:))];
[sk_x,sk_y,sk_z] = ellipsoid(sk_dip_locs(1,1),sk_dip_locs(1,2),sk_dip_locs(1,3),...
    sk_dip_locs(2,1),sk_dip_locs(2,2),sk_dip_locs(2,3));

plot3(uon_dip_locs(1,1),uon_dip_locs(1,2),uon_dip_locs(1,3),'r.')
h7 = surf(uon_x,uon_y,uon_z,'EdgeColor','none','FaceColor','r','FaceAlpha',0.3);
plot3(uon_dip_locs_skchans(1,1),uon_dip_locs_skchans(1,2),uon_dip_locs_skchans(1,3),'g.')
h8 = surf(uonskchans_x,uonskchans_y,uonskchans_z,'FaceColor','g','EdgeColor','none','FaceAlpha',0.3);
plot3(sk_dip_locs(1,1),sk_dip_locs(1,2),sk_dip_locs(1,3),'.','Color',[0.8500 0.3250 0.0980])
h9 = surf(sk_x,sk_y,sk_z,'EdgeColor','none','FaceColor','b','FaceAlpha',0.3);

[r,p] = ttest2(all_dip_locs(inds_na,:),all_dip_locs(inds_sa,:))
%% Plot age group ellipsoid of peak locations overlaid
figure(2);ft_plot_mesh(meshes,'facecolor',[.5 .5 .5],'facealpha',.3,'edgecolor','none')
hold on

age = data.age;
age(77) = [];
[age_order, age_inds] = sort(age,'ascend');
inds_2to4 = age_inds(age_order<=4);
inds_5to8 = age_inds(age_order>4&age_order<9);
inds_9to13 = age_inds(age_order>8&age_order<14);
inds_21to24 = age_inds(age_order>20&age_order<25);
inds_25to28 = age_inds(age_order>24&age_order<29);
inds_29to34 = age_inds(age_order>28&age_order<35);

all_dip_locs(77,:)=[];

dip_locs_2to4 = all_dip_locs(inds_2to4,:);
dip_locs_5to8 = all_dip_locs(inds_5to8,:);
dip_locs_9to13 = all_dip_locs(inds_9to13,:);
dip_locs_21to24 = all_dip_locs(inds_21to24,:);
dip_locs_25to28 = all_dip_locs(inds_25to28,:);
dip_locs_29to34 = all_dip_locs(inds_29to34,:);

%% Get mean and error for ellipsoid spheres for each age

age2to4 = [mean(all_dip_locs(inds_2to4,:)); std(all_dip_locs(inds_2to4,:))];
[age2to4_x,age2to4_y,age2to4_z] = ellipsoid(age2to4(1,1),age2to4(1,2),age2to4(1,3),...
    age2to4(2,1),age2to4(2,2),age2to4(2,3));
figure(2)
plot3(age2to4(1,1),age2to4(1,2),age2to4(1,3),'k.')
h1 = surf(age2to4_x,age2to4_y,age2to4_z,'EdgeColor','none','FaceColor','b','FaceAlpha',0.3);

age5to8 = [mean(all_dip_locs(inds_5to8,:)); std(all_dip_locs(inds_5to8,:))];
[age5to8_x,age5to8_y,age5to8_z] = ellipsoid(age5to8(1,1),age5to8(1,2),age5to8(1,3),...
    age5to8(2,1),age5to8(2,2),age5to8(2,3));
figure(2)
plot3(age5to8(1,1),age5to8(1,2),age5to8(1,3),'k.')
h2 = surf(age5to8_x,age5to8_y,age5to8_z,'EdgeColor','none','FaceColor','r','FaceAlpha',0.3);

age9to13 = [mean(all_dip_locs(inds_9to13,:)); std(all_dip_locs(inds_9to13,:))];
[age9to13_x,age9to13_y,age9to13_z] = ellipsoid(age9to13(1,1),age9to13(1,2),age9to13(1,3),...
    age9to13(2,1),age9to13(2,2),age9to13(2,3));
figure(2)
plot3(age9to13(1,1),age9to13(1,2),age9to13(1,3),'k.')
h3 = surf(age9to13_x,age9to13_y,age9to13_z,'EdgeColor','none','FaceColor','g','FaceAlpha',0.3);

age21to24 = [mean(all_dip_locs(inds_21to24,:)); std(all_dip_locs(inds_21to24,:))];
[age21to24_x,age21to24_y,age21to24_z] = ellipsoid(age21to24(1,1),age21to24(1,2),age21to24(1,3),...
    age21to24(2,1),age21to24(2,2),age21to24(2,3));
figure(2)
plot3(age21to24(1,1),age21to24(1,2),age21to24(1,3),'k.')
h4 = surf(age21to24_x,age21to24_y,age21to24_z,'EdgeColor','none','FaceColor','m','FaceAlpha',0.3);

age25to28 = [mean(all_dip_locs(inds_25to28,:)); std(all_dip_locs(inds_25to28,:))];
[age25to28_x,age25to28_y,age25to28_z] = ellipsoid(age25to28(1,1),age25to28(1,2),age25to28(1,3),...
    age25to28(2,1),age25to28(2,2),age25to28(2,3));
figure(2)
plot3(age25to28(1,1),age25to28(1,2),age25to28(1,3),'k.')
h5 = surf(age25to28_x,age25to28_y,age25to28_z,'EdgeColor','none','FaceColor','y','FaceAlpha',0.3);

age29to34 = [mean(all_dip_locs(inds_29to34,:)); std(all_dip_locs(inds_29to34,:))];
[age29to34_x,age29to34_y,age29to34_z] = ellipsoid(age29to34(1,1),age29to34(1,2),age29to34(1,3),...
    age29to34(2,1),age29to34(2,2),age29to34(2,3));
figure(2)
plot3(age29to34(1,1),age29to34(1,2),age29to34(1,3),'k.')
h6 = surf(age29to34_x,age29to34_y,age29to34_z,'EdgeColor','none','FaceColor','c','FaceAlpha',0.3);
legend('','','','age 2-4','','age 5-8','','age 9-13','','age 21-24','','age 25-28','','age 29-34')
%% stat test on dipole location lateralisation

[r_x,p_x] = corr(age,all_dip_locs(:,1))
[r_y,p_y] = corr(age,all_dip_locs(:,2))
[r_z,p_z] = corr(age,all_dip_locs(:,3))


