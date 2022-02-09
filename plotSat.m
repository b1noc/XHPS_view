function nodes = plotSatellite(et,nt)


	n = max(et(:,1));
	nodes = [];

    for k = 1:n
        % X,Y,Z coordinates of 1. node from node table
        n1x =  nt(et(k,6),1);
        n1y =  nt(et(k,6),2);
        n1z =  nt(et(k,6),3);
        N1  = [n1x;n1y;n1z];

        % X,Y,Z coordinates of 2. node from node table
        n2x =  nt(et(k,7),1);
        n2y =  nt(et(k,7),2);
        n2z =  nt(et(k,7),3);
        N2  = [n2x;n2y;n2z];

        % X,Y,Z coordinates of 3. node from node table
        n3x =  nt(et(k,8),1);
        n3y =  nt(et(k,8),2);
        n3z =  nt(et(k,8),3);
        N3  = [n3x;n3y;n3z];

        % X,Y,Z coordinates of 4. node from node table
        n4x =  nt(et(k,9),1);
        n4y =  nt(et(k,9),2);
        n4z =  nt(et(k,9),3);
        N4  = [n4x;n4y;n4z];

        nodes = [nodes; [N1,N2,N3,N4]'];
        %fill3(nodes(:,1),nodes(:,2),nodes(:,3),'y');
    end
end


