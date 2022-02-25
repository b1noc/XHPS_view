function nodes = plotSatellite(path, q_k)

	et = load([path, 'et.txt.']);
	nt = load([path, 'nt.txt.']);

	n = max(et(:,1));
	nodes = [];
	r_com = [1.5615 0.9720 -0.3310];

    for k = 1:n

		%transform inertial kos to body fixed by quaternion multiplication
		%(xrvar qnefgryyhat qrf Beovgny Xbs in eci daher vecbyquattransposed)

        % X,Y,Z coordinates of 1. node from node table
        nx =  nt(et(k,6),1);
        ny =  nt(et(k,6),2);
        nz =  nt(et(k,6),3);

		N1 = [nx ny nz] - r_com;
        N1t  = HPS_transformVecByQuatTransposed(N1,q_k);

        % X,Y,Z coordinates of 2. node from node table
        nx =  nt(et(k,7),1);
        ny =  nt(et(k,7),2);
        nz =  nt(et(k,7),3);


		N2 = [nx ny nz] - r_com;
        N2t  = HPS_transformVecByQuatTransposed(N2,q_k);


        % X,Y,Z coordinates of 3. node from node table
        nx =  nt(et(k,8),1);
        ny =  nt(et(k,8),2);
        nz =  nt(et(k,8),3);
        
		N3 = [nx ny nz] - r_com;
        N3t  = HPS_transformVecByQuatTransposed(N3,q_k);


        % X,Y,Z coordinates of 4. node from node table
        nx =  nt(et(k,9),1);
        ny =  nt(et(k,9),2);
        nz =  nt(et(k,9),3);
        
		N4 = [nx ny nz] - r_com;
        N4t  = HPS_transformVecByQuatTransposed(N4,q_k);

		nodes = [nodes; [N1t,N2t,N3t,N4t]'];
    end
end


