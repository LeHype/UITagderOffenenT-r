function [] = third(name)


% clc
load('Workspace.mat')
fig_video = MakeDefaultFig(figureData.width,figureData.height,'CallbackListenerFunOn',0,'Screen',figureData.screen,'FigTitle','Video','FigPos',figureData.position);
set(gcf,'DefaultAxesPosition',[0 0 1 1])

if ~exist(['Video' filesep],'dir')
    mkdir('Video')
end

disp('------------------------------')
disp(['Danke für deine Teilnahme ',name]);
disp('Drücke auf Reset um von Vorne zu beginnen')
disp('------------------------------')

movie(frames_opt,1,VideoFPS_actual*VideoSpeed);

writeFramesToVideo(frames_opt,vidObj)
movefile(vidObj.Filename,['Video' filesep])

pause(1.5)
set(fig_video,'WindowState','minimized')

end

