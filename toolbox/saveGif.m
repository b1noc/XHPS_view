function [] = saveGif(name, frameVec)
	filename = append(name, '.gif')
	for i = 1:length(frameVec)
		im = frame2im(frameVec(i)); 
		[imind,cm] = rgb2ind(im,256); 
		% Write to the GIF File 
		if i == 1 
		  imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
		else 
		  imwrite(imind,cm,filename,'gif','WriteMode','append'); 
		end 
	end
end
