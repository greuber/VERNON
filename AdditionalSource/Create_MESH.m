function [ NUM,MESH ] = Create_MESH( NUM,PAR,random,MESH )

clear j1 j2 i

% initialize x and z
x = [0:NUM.NUMERICS.dx:PAR.W];
z = [0:NUM.NUMERICS.dz:PAR.H];
if random == 1
    x = x + 0.1*[rand(size(x))-0.5]*NUM.NUMERICS.dx;             % Add random distribution additionally for variable spacing
    z = z + 0.1*[rand(size(z))-0.5]*NUM.NUMERICS.dz;             % Add random distribution additionally for variable spacing
end

MESH.GCOORD = zeros(2,NUM.NUMERICS.no_nodes_z);
for i=1:1:NUM.NUMERICS.no_nodes_z
    j1 = ((i-1)*NUM.NUMERICS.no_nodes_x)+1;
    j2 = ((i)*NUM.NUMERICS.no_nodes_x);
    MESH.GCOORD(1,j1:j2) = x;
    MESH.GCOORD(2,j1:j2) = z(i);
end


% Choose perturbation
switch PAR.Perturb
    case 'middlesmooth'
        % make sinus in the middle of the model with smoothing to the boarders
        st                              = PAR.A0/sin(pi/2);
        i                               = length(MESH.GCOORD(1,:)):-NUM.NUMERICS.no_nodes_x:NUM.NUMERICS.no_nodes_x;
        factor                          = zeros(1,length(i));
        
        factor(1,1:(length(i)-1)/2)     = linspace(1e-10,st,(length(i)-1)/2);
        factor(1,(length(i)+1)/2:end)   = linspace(st,1e-10,(length(i)+1)/2);
        
        j = length(MESH.GCOORD(1,:))+1-NUM.NUMERICS.no_nodes_x;
        k = 1;
        for i = length(MESH.GCOORD(1,:)):-NUM.NUMERICS.no_nodes_x:NUM.NUMERICS.no_nodes_x
            k_topo          = linspace((0),(1*pi),NUM.NUMERICS.no_nodes_x);
            %     k_topo          = k_topo + 0.2*[rand(size(k_topo))-0.5];
            z_topo          = [sin(k_topo)];
            z_topo          = factor(1,k).*z_topo;
            MESH.GCOORD(2,j:i)   = MESH.GCOORD(2,j:i) + [z_topo];
            j               = j-NUM.NUMERICS.no_nodes_x;
            k               = k+1;
            
        end
        
    case 'middleconst'
        % make sinus in the middle of the model with constant sinus to the top
        st                              = PAR.A0;
        i                               = length(MESH.GCOORD(1,:)):-NUM.NUMERICS.no_nodes_x:NUM.NUMERICS.no_nodes_x;
        factor                          = zeros(1,length(i));
        
        ind = find(MESH.GCOORD(2,:) < floor(PAR.H/2)+NUM.NUMERICS.dz-1 & MESH.GCOORD(2,:) > floor(PAR.H/2)-NUM.NUMERICS.dz-1);
        ind = floor(min(ind)/NUM.NUMERICS.no_nodes_x);
        factor(1,1:ind)     = st;
        factor(1,ind:end)   = 1;
        
        j = length(MESH.GCOORD(1,:))+1-NUM.NUMERICS.no_nodes_x;
        k = 1;
        for i = length(MESH.GCOORD(1,:)):-NUM.NUMERICS.no_nodes_x:NUM.NUMERICS.no_nodes_x
            k_topo          = linspace((0),(1*pi),NUM.NUMERICS.no_nodes_x);
            %     k_topo          = k_topo + 0.2*[rand(size(k_topo))-0.5];
            z_topo          = [sin(k_topo)];
            z_topo          = factor(1,k).*z_topo;
            MESH.GCOORD(2,j:i)   = MESH.GCOORD(2,j:i) + [z_topo];
            j               = j-NUM.NUMERICS.no_nodes_x;
            k               = k+1;
            
        end
        
    case 'middle'
        % make sinus only in the middle of the model
        % make sinus only in the middle of the model
        st                              = PAR.A0/sin(pi/2);
        k_topo          = linspace((0),(pi),NUM.NUMERICS.no_nodes_x);
        z_topo          = [sin(k_topo)];
        z_topo          = st.*z_topo;
        MESH.GCOORD(2,((round((NUM.NUMERICS.no_nodes_z-1)/2))*NUM.NUMERICS.no_nodes_x)+1:(round((NUM.NUMERICS.no_nodes_z)/2))*NUM.NUMERICS.no_nodes_x)   = MESH.GCOORD(2,((round((NUM.NUMERICS.no_nodes_z-1)/2))*NUM.NUMERICS.no_nodes_x)+1:(round((NUM.NUMERICS.no_nodes_z)/2))*NUM.NUMERICS.no_nodes_x) + [z_topo];
        
        
    case 'middleinterface'
        % to have more than one interface: PAR.H_interface = [hi:NUM.NUMERICS.dz:hi2];
        % make sinus only at a specified interface in PAR.H_interface
        for i = 1:length(PAR.H_interface)
            bound = NUM.NUMERICS.dz - 3*NUM.NUMERICS.dz/4;
            st              = PAR.A0/sin(pi/2);
            k_topo          = linspace((0),(pi/2),NUM.NUMERICS.no_nodes_x);
            z_topo          = [sin(k_topo)];
            z_topo          = st.*z_topo;
            ind = find(MESH.GCOORD(2,:)<PAR.H_interface(i)+bound & MESH.GCOORD(2,:)>PAR.H_interface(i)-bound);
            if length(ind) ~= length(z_topo)
                while length(ind) > length(z_topo)
                    bound = bound/2;
                    ind = find(MESH.GCOORD(2,:)<PAR.H_interface(i)+bound & MESH.GCOORD(2,:)>PAR.H_interface(i)-bound);
                end
                while length(ind) < length(z_topo)
                    bound = bound*2;
                    ind = find(MESH.GCOORD(2,:)<PAR.H_interface(i)+bound & MESH.GCOORD(2,:)>PAR.H_interface(i)-bound);
                end
            end
            MESH.GCOORD(2,ind)   = MESH.GCOORD(2,ind) + [z_topo];
        end
        
        
    case 'random'
        % make random disribution for the interface parameter
        % PAR.H_interface
        if ~isfield(PAR,'H_interface')
            PAR.H_interface = PAR.H/2;
        end
        bound = NUM.NUMERICS.dz - NUM.NUMERICS.dz/2;
        ind = find(MESH.GCOORD(2,:)<PAR.H_interface+bound & MESH.GCOORD(2,:)>PAR.H_interface-bound);
        MESH.GCOORD(2,ind) = MESH.GCOORD(2,ind) + bound*[rand(size(MESH.GCOORD(2,ind)))]/5;
        
    case 'none'
end


end

