%% MNI Plot per subject

datdir = (['projects/p31578/DATA/WSU/' subid '/'])
load([datdir 'data_montage.mat']);

elec_mni_frv_2 = [];
elec_mni_frv_2 = data.elec; % run if different path and if loaded BPM data, do not need to run rest of block
elec_mni_frv_2.chanpos = montage.coordinates_new;
elec_mni_frv_2.elecpos = elec_mni_frv_2.chanpos;
elec_mni_frv_2.label = 1:length(elec_mni_frv_2.chanpos);
elec_mni_frv_2.label = elec_mni_frv_2.label';

surface_template_l = load('/projects/p31578/smg1656/Recons/ft_templates/surface_pial_left.mat'); % from ft
surface_template_r = load('/projects/p31578/smg1656/Recons/ft_templates/surface_pial_right.mat');
ft_plot_mesh(surface_template_l.mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', 0.3);
ft_plot_mesh(surface_template_r.mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', 0.3);


view([0 90]); lighting gouraud; camlight; material dull; 
set(gcf, 'color', 'w')

ft_plot_sens(elec_mni_frv_2, 'elecshape', 'sphere', 'color', 'red');

