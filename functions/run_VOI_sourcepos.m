%% save out VOI sourcepos
for ss = 1:27
clearvars -except ss
close all
clc
addpath('R:\DRS-KidsOPM\Paediatric_OPM_Notts\fieldtrip-20220906')
ft_defaults
%%

project_dir = 'R:\DRS-KidsOPM\Paediatric_OPM_Notts\';
datadir = [project_dir,'Data',filesep,'BIDS',filesep];
run = 'run-001';
s1 = sprintf('%2d',ss);s1(s1 == ' ') = '0'
    sub = strcat('0',s1);
    ses_i = 1;
    ses = sprintf('%3d',ses_i);ses(ses == ' ') = '0';

path_meshes = [datadir,'derivatives',filesep,'sourcespace',filesep,'sub-',sub,filesep];
files_VOI = ['sub-',sub,'_AAL_VOI.mat'];
files_AAL_regions = ['sub-',sub,'_AAL_regions.nii.gz'];


% load AAL locations
VOI_mat_file = [path_meshes,files_VOI];
AAL_regions = ft_read_mri([path_meshes,files_AAL_regions]);
AAL_regions = ft_convert_units(AAL_regions,'m');
    
    VOI = (AAL_regions.anatomy(:,:,:,26) > 0) | (AAL_regions.anatomy(:,:,:,65)>0);

    SE = strel('sphere',5);
    dilatedVOI = imdilate(VOI,SE) & (AAL_regions.anatomy(:,:,:,79) > 0);
    [sourcepos_vox(:,1),sourcepos_vox(:,2),sourcepos_vox(:,3)] = ind2sub(AAL_regions.dim, find(dilatedVOI));
    
    sourcepos = ft_warp_apply(AAL_regions.transform,sourcepos_vox);
    save([path_meshes,files_VOI(1:end-4) '_sourcepos_vm.mat'],'sourcepos')
end 