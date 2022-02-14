function nodes = plotSatellite(et,nt, q_k)


	n = max(et(:,1));
	nodes = [];

    for k = 1:n

		%transform inertial kos to body fixed by quaternion multiplication
		%(xrvar qnefgryyhat qrf Beovgny Xbs in eci daher vecbyquattransposed)

        % X,Y,Z coordinates of 1. node from node table
        nx =  nt(et(k,6),1);
        ny =  nt(et(k,6),2);
        nz =  nt(et(k,6),3);

        N1  = HPS_transformVecByQuatTransposed([nx;ny;nz],q_k);

        % X,Y,Z coordinates of 2. node from node table
        nx =  nt(et(k,7),1);
        ny =  nt(et(k,7),2);
        nz =  nt(et(k,7),3);

        N2  = HPS_transformVecByQuatTransposed([nx;ny;nz],q_k);


        % X,Y,Z coordinates of 3. node from node table
        nx =  nt(et(k,8),1);
        ny =  nt(et(k,8),2);
        nz =  nt(et(k,8),3);
        
        N3  = HPS_transformVecByQuatTransposed([nx;ny;nz],q_k);


        % X,Y,Z coordinates of 4. node from node table
        nx =  nt(et(k,9),1);
        ny =  nt(et(k,9),2);
        nz =  nt(et(k,9),3);
        
        N4  = HPS_transformVecByQuatTransposed([nx;ny;nz],q_k);


        nodes = [nodes; [N1,N2,N3,N4]'];
        %fill3(nodes(:,1),nodes(:,2),nodes(:,3),'y');
    end
end


