function [] = third(name)


% clc
load('Workspace.mat','fig_height','fig_width','fig_pos','frames','video_speed','fps','screen');
 fig_video = MakeDefaultFig(fig_width,fig_height,'CallbackListenerFunOn',0,'Screen',screen,'FigTitle','Video','FigPos',fig_pos);
 set(gcf,'DefaultAxesPosition',[0 0 1 1])
v = VideoWriter([name,'_optimization']);
v.FrameRate = video_speed*fps;
% video_title = [sprintf('%0.5fs',t_f_opt),'_',name_person,'_at_',datestr(now,'hhMM')];
    
open(v);
writeVideo(v,frames);
close(v);
 disp(['Danke für deine Teilnahme ',name]);
 disp('Drücke auf Reset um von Vorne zu beginnen')

 movie(frames,1,round(video_speed*fps));

end

