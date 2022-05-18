% get frames of F via
%     F(i) = getframe(gcf) ;
%     drawnow

% seems not to work  --> crashes 
% fname = 'help_pic';
% print('-dpng','-r400',fname)     
% read_success = false;
% while read_success==false
%     try
%         I = imread([fname '.png']);
%         read_success = true;
%     catch
%         read_success = false;
%     end
% end
% F(i) = im2frame(I);

function make_video(F,framerate,video_title,directory)
    if nargin<4
        directory = '';
    end   
    % create the video writer with 1 fps
    writerObj = VideoWriter([directory,video_title,'.avi']);
    % set the seconds per image
    writerObj.FrameRate = framerate;
    % open the video writer
    open(writerObj);
    % write the frames to the video
    for i=1:length(F)
        % convert the image to a frame
        frame = F(i) ;    
        writeVideo(writerObj, frame);
    end
    % close the writer object
    close(writerObj);
end