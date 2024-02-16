function [TFS] = VE_TFS(VE_unchopped,OFF,trig_offset,duration,ind,Ntrials,trial_time,f)
cbar_lim = [0.5]; % colour bar limit (relative change)
highpass = [1 2 4 6 8 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100 105 110];
lowpass = [4 6 8 10 13 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100 105 110 115 120];
fre = highpass + ((lowpass - highpass)./2);
% Control window
conwin = OFF*f-trig_offset*f;
fx = figure;
fx.Name = 'VE TFS';
% Filter data within bands and calculate envelope
VE_fb = zeros(length(VE_unchopped),length(fre));
for fb = 1:length(highpass)
    filt_VE = nut_filter3(VE_unchopped','butter','bp',4,highpass(fb),lowpass(fb),f,1)';
    VE_fb(:,fb) = abs(hilbert(filt_VE));
end
VE_mean = zeros(duration*f,length(fre));
for fb = 1:length(highpass)
    VE_fb_trials = [];
    % Chop data
    for i = 1:Ntrials
        VE_fb_trials = cat(1,VE_fb_trials,VE_fb(ind(i)+1:ind(i)+(duration*f),fb));
    end
    VE_filt = reshape(VE_fb_trials,duration*f,Ntrials);
    % Average across trials
    VE_mean(:,fb) = mean(VE_filt,2);
end
meanrest = mean(VE_mean(conwin(1):conwin(2),:),1);
meanrestmat = repmat(meanrest,size(VE_mean,1),1);
TFS = (VE_mean'-meanrestmat')./meanrestmat';
figure(fx)
pcolor(trial_time,fre,TFS);shading interp
xlabel('Time (s)');ylabel('Frequency (Hz)')
colorbar;caxis([-cbar_lim cbar_lim])
axis fill
end
