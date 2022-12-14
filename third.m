function [] = third(name)


% clc
load('Workspace.mat','fig_height','fig_width','fig_pos','frames_opt','video_speed','fps','screen','t_f_opt');
 fig_video = MakeDefaultFig(fig_width,fig_height,'CallbackListenerFunOn',0,'Screen',screen,'FigTitle','Video','FigPos',fig_pos);
 set(gcf,'DefaultAxesPosition',[0 0 1 1])

 if ~exist(['Video' filesep],'dir')
     mkdir('Video')
 end
 

  cd (['Video' filesep])
 v = VideoWriter([num2str(t_f_opt),'_s_Fahrtzeit_',name]);
v.FrameRate = video_speed*fps;
% video_title = [sprintf('%0.5fs',t_f_opt),'_',name_person,'_at_',datestr(now,'hhMM')];
    
open(v);
writeVideo(v,frames_opt);
close(v);
disp('------------------------------')
 disp(['Danke für deine Teilnahme ',name]);
 disp('Drücke auf Reset um von Vorne zu beginnen')
disp('------------------------------')
 movie(frames_opt,1,round(video_speed*fps));
cd ..
end

