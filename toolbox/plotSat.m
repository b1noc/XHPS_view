function patchData = plotSatellite(path, q_k)
	et = load([path, 'et.txt.']);
	nt = load([path, 'nt.txt.']);

	n = max(et(:,1));
	r_com = [1.5615 0.9720 -0.3310];

	nodes = [];
	x = []; y = []; z = [];

	faces = [et(:,6) et(:,7) et(:,8) et(:,9)];
	vertices = nt;
	patch('Vertices', nt, 'Faces', faces)

    for k = 1:n
        % X,Y,Z coordinates of 1. node from node table
        nx =  nt(et(k,6),1);
        ny =  nt(et(k,6),2);
        nz =  nt(et(k,6),3);

		N1 = [nx; ny; nz] - r_com';

        % X,Y,Z coordinates of 2. node from node table
        nx =  nt(et(k,7),1);
        ny =  nt(et(k,7),2);
        nz =  nt(et(k,7),3);

		N2 = [nx; ny; nz] - r_com';


        % X,Y,Z coordinates of 3. node from node table
        nx =  nt(et(k,8),1);
        ny =  nt(et(k,8),2);
        nz =  nt(et(k,8),3);

		N3 = [nx; ny; nz] - r_com';

        % X,Y,Z coordinates of 4. node from node table
        nx =  nt(et(k,9),1);
        ny =  nt(et(k,9),2);
        nz =  nt(et(k,9),3);
        
		N4 = [nx; ny; nz] - r_com';

		nodes = [N1,N2,N3,N4]';
		x = [x nodes(:,1)]
		y = [y nodes(:,2)]
		z = [z nodes(:,3)]
    end
	patchData = [x; y; z];
	figure
	plot_sat = patch(x,y,z);
	plot_sat.FaceColor = 'y';
	figure
end


