function convert_tstat_units(sub,ses,project_dir,tstat_type)

tstat = ft_read_mri([project_dir 'Data' filesep 'BIDS' filesep 'derivatives' filesep 'Tstats' filesep 'sub-' sub filesep 'sub-' sub '_ses-' ses '_task-faces_circles_run-001_pseudoT_' tstat_type '_adj.nii']);
tstat.unit = 'm';
tstat = ft_convert_units(tstat,'mm');
ft_write_mri([project_dir 'Data' filesep 'BIDS' filesep 'derivatives' filesep 'Tstats' filesep 'sub-' sub filesep 'sub-' sub '_ses-' ses '_task-faces_circles_run-001_pseudoT_' tstat_type '_adj_mm.nii'],tstat,'dataformat','nifti')

end
