function [mean_Env] = get_envelope(trial_data,C,Lead_fields,min_voxel_index,control_inds)
mu = 0.05;
Cr = C + mu*max(svd(C))*eye(size(C));
Cr_inv = inv(Cr);

this_L = Lead_fields(:,:,min_voxel_index);
W_v = inv((this_L'*Cr_inv*this_L))*(this_L'*Cr_inv);
iPower_v = this_L'*Cr_inv*this_L;
[v,d] = svd(iPower_v);
[~,id] = min(diag(d));
lopt = this_L*v(:,id); % turn to nAm amplitude
w = (lopt'*Cr_inv/(lopt'*Cr_inv*lopt))';
N_trials = size(trial_data,3);
VEs = zeros(size(trial_data,2),N_trials);
Envs = VEs;
for tr_i = 1:N_trials
    VEs(:,tr_i) = w'*trial_data(:,:,tr_i) ./sqrt(w'*w);
    Envs(:,tr_i) = abs(hilbert(VEs(:,tr_i)));
    Envs(:,tr_i) = (Envs(:,tr_i) - mean(Envs(control_inds,tr_i)))./mean(Envs(control_inds,tr_i));
end
mean_Env = mean(Envs,2);
end