function [] = third(name)


% clc
load('Workspace.mat','fig_height','fig_width','fig_pos','frames','video_speed','fps','screen','t_f_opt');
 cd Video/
fig_video = MakeDefaultFig(fig_width,fig_height,'CallbackListenerFunOn',0,'Screen',screen,'FigTitle','Video','FigPos',fig_pos);
 set(gcf,'DefaultAxesPosition',[0 0 1 1])
v = VideoWriter([num2str(t_f_opt),'_s_Fahrtzeit_',name]);
v.FrameRate = video_speed*fps;
% video_title = [sprintf('%0.5fs',t_f_opt),'_',name_person,'_at_',datestr(now,'hhMM')];
    
open(v);
writeVideo(v,frames);
close(v);
disp('------------------------------')
 disp(['Danke für deine Teilnahme ',name]);
 disp('Drücke auf Reset um von Vorne zu beginnen')
disp('------------------------------')
 movie(frames,1,round(video_speed*fps));
cd ..
end

