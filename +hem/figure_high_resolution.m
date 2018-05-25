%plot high resolution images

fig = gcf; %grab current displayed figure
% set(gcf,'InvertHardcopy','off') %toggle between this to adjust background
% color
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 4.5]; %size of image default [0 0 6 4.5]
fig.PaperPositionMode = 'manual';
print('Picture9', '-dtiff', '-r600'); %name of figure, file type and resolution