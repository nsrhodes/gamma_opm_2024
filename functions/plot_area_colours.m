function plot_area_colours(Time, cannonall_mat_b, cannonall_mat_b_se, colour, colourvec)
X = [Time, fliplr(Time)];
Y = [cannonall_mat_b' + cannonall_mat_b_se', fliplr(cannonall_mat_b' - cannonall_mat_b_se')];
fill(X,Y,colourvec,'FaceAlpha',0.3,'EdgeColor','none');
hold on
plot(Time,cannonall_mat_b,'Color',colour,'LineWidth',1)
