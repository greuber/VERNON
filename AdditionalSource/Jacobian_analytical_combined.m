function [ NUM ] = Jacobian_analytical_combined( NUM,PAR,MESH )

% once used to combine get_globals_jacobian_analytical and
% get_globals_picard but this just results in Picard method!!!

%% initilize globals
NUM.Solve.J = zeros(NUM.no_nodes*(NUM.ndof-1)+NUM.no_nodes_linear,NUM.no_nodes*(NUM.ndof-1)+NUM.no_nodes_linear);
NUM.Solve.L = zeros(NUM.no_nodes*(NUM.ndof-1)+NUM.no_nodes_linear,NUM.no_nodes*(NUM.ndof-1)+NUM.no_nodes_linear);
NUM.Solve.FG = zeros(NUM.no_nodes*(NUM.ndof-1)+NUM.no_nodes_linear,1);
NUM.Solve.B  = zeros(3,NUM.no_nodes_ele*2);
NUM.Solve.F  = zeros(NUM.no_nodes_ele*(NUM.ndof-1)+NUM.no_nodes_ele_linear,1);
NUM.Solve.f_res = zeros(NUM.no_nodes*(NUM.ndof-1)+NUM.no_nodes_linear,1);
NUM.mu_2d = zeros(1,NUM.no_nodes)   ;                                       % viscosity for this elemet for plotting later



for i = 1:NUM.no_elems_global
    
    % initialize local matrices
    NUM.Solve.KM_L = zeros(NUM.no_nodes_ele*(NUM.ndof-1),NUM.no_nodes_ele*(NUM.ndof-1));
    NUM.Solve.KM_Jac = zeros(NUM.no_nodes_ele*(NUM.ndof-1),NUM.no_nodes_ele*(NUM.ndof-1));
    NUM.Solve.F  = zeros(NUM.no_nodes_ele*(NUM.ndof-1)+NUM.no_nodes_ele_linear,1);
    NUM.Solve.GM = zeros(NUM.no_nodes_ele*(NUM.ndof-1),NUM.no_nodes_ele_linear);
    NUM.Solve.LM = zeros(NUM.no_nodes_ele*(NUM.ndof-1)+NUM.no_nodes_ele_linear,NUM.no_nodes_ele*(NUM.ndof-1)+NUM.no_nodes_ele_linear);
    
    
    
    LGCOORD = (MESH.GCOORD(:,NUM.Number.number_quad(:,i)))';      % local global coordinates for this element
    
    for j =1:NUM.no_intp
        
        [ NI,dNI ] = Shape_Functions( NUM.no_nodes_ele,j,MESH.INTP.COORD );
        NUM.Solve.N_matrix = [NI(1,1), 0, NI(1,2), 0, NI(1,3), 0, NI(1,4), 0, NI(1,5), 0, NI(1,6), 0, NI(1,7), 0, NI(1,8), 0, NI(1,9),0;
            0, NI(1,1), 0, NI(1,2), 0, NI(1,3), 0, NI(1,4), 0, NI(1,5), 0, NI(1,6), 0, NI(1,7), 0, NI(1,8), 0, NI(1,9)];
        
        
        NUM.Solve.N   = NI;
        NUM.Solve.dNds = dNI;
        NUM.Solve.Jac    = NUM.Solve.dNds*LGCOORD;
        
        % Inverse of Jacobian
        NUM.Solve.invJ = inv(NUM.Solve.Jac);
        NUM.Solve.detJ = det(NUM.Solve.Jac);
        NUM.Solve.dNdX = NUM.Solve.invJ*NUM.Solve.dNds;
        
        
        NUM.Solve.B(1,1:2:end-1) = NUM.Solve.dNdX(1,:);
        NUM.Solve.B(2,2:2:end)   = NUM.Solve.dNdX(2,:);
        NUM.Solve.B(3,1:2:end-1) = NUM.Solve.dNdX(2,:);
        NUM.Solve.B(3,2:2:end)   = NUM.Solve.dNdX(1,:);
        
        
        
        
        [ NUM] = ComputeViscosity(PAR,i,NUM);
        NUM.mu_2d(NUM.Number.number_quad(j,i)) = NUM.mu;
        
        NUM.Solve.P = (1/2)*[ 2 ,        0               ,      0,0
            0                       , 2 ,      0,0
            0                       ,        0             , 1,1];
        
        % Compute D matrix for global L
        NUM.Solve.D_L = 2*NUM.mu*NUM.Solve.P*NUM.Solve.P';
        
        % Compute D matrix for global Jacobian
        I = eye(4,4);
        beta = (1/2)*(1/PAR.powerlaw(NUM.Number.number_quad(9,i)) - 1);
        str_local_dash = NUM.Solve.P'*NUM.str_local;
        nu = (1/NUM.str_invariant).*str_local_dash;
        NUM.Solve.D_Jac = 2*NUM.mu*NUM.Solve.P*(I + beta*nu*nu')*NUM.Solve.P';
        
        % compute KM for L and Jacobian
        NUM.Solve.KM_L = NUM.Solve.KM_L + NUM.Solve.B'*NUM.Solve.D_L*NUM.Solve.B * MESH.INTP.weight(j)*NUM.Solve.detJ;
        NUM.Solve.KM_Jac = NUM.Solve.KM_Jac + NUM.Solve.B'*NUM.Solve.D_L*NUM.Solve.B * MESH.INTP.weight(j)*NUM.Solve.detJ;
        
        NUM.Solve.F(1:end-4,1) = NUM.Solve.F(1:end-4,1) + NUM.Solve.N_matrix'*PAR.rho(NUM.Number.number_quad(9,i))' * PAR.g * MESH.INTP.weight(j) * NUM.Solve.detJ;
        % Compute linear shapw functions
        NUM.no_nodes_ele = 4;
        [ NI,dNI ] = Shape_Functions( NUM.no_nodes_ele,j,MESH.INTP.COORD );
        NUM.Solve.N_p    = NI;
        NUM.Solve.GM = NUM.Solve.GM - (NUM.Solve.B'*NUM.Solve.m*NUM.Solve.N_p* MESH.INTP.weight(j)*NUM.Solve.detJ);
        NUM.no_nodes_ele = 9;
        
        
    end
    
    NUM.Solve.LM_L(1:18,1:18) = NUM.Solve.KM_L;
    NUM.Solve.LM_L(1:18,19:22) = NUM.Solve.GM;
    NUM.Solve.LM_L(19:22,1:18) = NUM.Solve.GM';
    NUM.Solve.LM_Jac(1:18,1:18) = NUM.Solve.KM_Jac;
    NUM.Solve.LM_Jac(1:18,19:22) = NUM.Solve.GM;
    NUM.Solve.LM_Jac(19:22,1:18) = NUM.Solve.GM';
    
    
    % add element stiffness matrix to global stiffness matrix
    NUM.Solve.L(NUM.Number.number_ele_dof(:,i),NUM.Number.number_ele_dof(:,i)) = NUM.Solve.L(NUM.Number.number_ele_dof(:,i),NUM.Number.number_ele_dof(:,i)) + NUM.Solve.LM_L;
    NUM.Solve.J(NUM.Number.number_ele_dof(:,i),NUM.Number.number_ele_dof(:,i)) = NUM.Solve.J(NUM.Number.number_ele_dof(:,i),NUM.Number.number_ele_dof(:,i)) + NUM.Solve.LM_Jac;
    NUM.Solve.FG(NUM.Number.number_ele_dof(:,i),1) = NUM.Solve.FG(NUM.Number.number_ele_dof(:,i),1) + NUM.Solve.F;
    
end

% Add boundaries to the system
for i = 1:1:length(NUM.Boundary.bcdof)
    
    NUM.Solve.L(NUM.Boundary.bcdof(i),:) = 0;
    NUM.Solve.L(NUM.Boundary.bcdof(i),NUM.Boundary.bcdof(i)) = 1;
    NUM.Solve.J(NUM.Boundary.bcdof(i),:) = 0;
    NUM.Solve.J(NUM.Boundary.bcdof(i),NUM.Boundary.bcdof(i)) = 1;
    NUM.Solve.FG(NUM.Boundary.bcdof(i)) = NUM.Boundary.bcval(i);
    
end


end



