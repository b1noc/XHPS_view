function [] = saveVid(name, dur, aniVec)

fr = length(aniVec)/dur;
wout = VideoWriter(name);
wout.FrameRate = fr;
open(wout);
writeVideo(wout,aniVec);
close(wout);

end
