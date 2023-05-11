function writeFramesToVideo(FRAMES,vidObj)

open(vidObj);
for ii=1:length(FRAMES)
    writeVideo(vidObj, FRAMES(ii));
end
close(vidObj);

end