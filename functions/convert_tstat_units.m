function convert_tstat_units(sub,ses,project_dir)

tstat = ft_read_mri([project_dir 'sub-' sub '/sub-' sub '_' ses '_task-faces_circles_run-001_pseudoT_circles_adj.nii']);
tstat.unit = 'm';
tstat = ft_convert_units(tstat,'mm');
ft_write_mri([project_dir 'sub-' sub '/sub-' sub '_' ses '_task-faces_circles_run-001_pseudoT_circles_adj_mm.nii'],tstat,'dataformat','nifti')

end
