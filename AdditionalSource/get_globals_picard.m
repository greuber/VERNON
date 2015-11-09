function [ NUM,MESH ] = get_globals_picard( NUM,PAR,MESH,CHAR)
% Evaluates the Picard matrix and the right-hand-side vector

NUM.Solve.FG = zeros(NUM.NUMERICS.no_nodes*(NUM.NUMERICS.ndof-1)+NUM.NUMERICS.no_nodes_linear,1);
vec = zeros((NUM.NUMERICS.no_nodes_ele*(NUM.NUMERICS.ndof-1)+NUM.NUMERICS.no_nodes_ele_linear)^2,NUM.NUMERICS.no_elems_global);

% Some ifs
if NUM.Plasticity.Plasticity
    NUM.Plasticity.Plastic = zeros(NUM.NUMERICS.no_intp,NUM.NUMERICS.no_elems_global)   ;
end

B = zeros(3,NUM.NUMERICS.no_nodes_ele*2);
m  = [1 1 0]';
P  = (1/2)*[ 2 , 0 , 0 , 0
            0 , 2 , 0 , 0
            0 , 0 , 1 , 1];
N_matrix = zeros(2,NUM.NUMERICS.no_nodes_ele*2);

for i = 1:NUM.NUMERICS.no_elems_global
    
    % initialize local matrices
    KM = zeros(NUM.NUMERICS.no_nodes_ele*(NUM.NUMERICS.ndof-1),NUM.NUMERICS.no_nodes_ele*(NUM.NUMERICS.ndof-1));
    F  = zeros(NUM.NUMERICS.no_nodes_ele*(NUM.NUMERICS.ndof-1)+NUM.NUMERICS.no_nodes_ele_linear,1);
    GM = zeros(NUM.NUMERICS.no_nodes_ele*(NUM.NUMERICS.ndof-1),NUM.NUMERICS.no_nodes_ele_linear);
    LM = zeros(NUM.NUMERICS.no_nodes_ele*(NUM.NUMERICS.ndof-1)+NUM.NUMERICS.no_nodes_ele_linear,NUM.NUMERICS.no_nodes_ele*(NUM.NUMERICS.ndof-1)+NUM.NUMERICS.no_nodes_ele_linear);
    
    
    
    LGCOORD = (MESH.GCOORD(:,NUM.Number.number_quad(:,i)))';      % local global coordinates for this element
    
    for j =1:NUM.NUMERICS.no_intp
        
        NI     = NUM.Shape.NI_vector{j};
        N_p    = NUM.Shape.pNI_vector{j};
                
        N_matrix(1,1:2:end) = NI(1,:);
        N_matrix(2,2:2:end) = NI(1,:);
        
        dNds = NUM.Shape.dNI_vector{j};
        Jac    = dNds*LGCOORD;
        
        % Inverse of Jacobian
        invJ = inv(Jac);
        detJ = det(Jac);
        dNdX = invJ*dNds;
        
        B(1,1:2:end) = dNdX(1,:);
        B(2,2:2:end) = dNdX(2,:);
        B(3,1:2:end) = dNdX(2,:);
        B(3,2:2:end) = dNdX(1,:);
        
        [ NUM,MESH] = ComputeViscosity(i,NUM,j,MESH,CHAR,PAR,B);
        
        D               = 2*NUM.Viscosity.mu*(P*P');
        KM              = KM + B'*D*B * MESH.INTP.weight(j)*detJ;
        GM              = GM - (B'*m*N_p* MESH.INTP.weight(j)*detJ);
        F(1:end-NUM.NUMERICS.no_nodes_ele_linear,1)    = F(1:end-NUM.NUMERICS.no_nodes_ele_linear,1) + N_matrix'*MESH.CompVar.rho(NUM.Number.number_quad(j,i))' * PAR.g * MESH.INTP.weight(j) * detJ;

    end
    
    LM(1:NUM.NUMERICS.no_nodes_ele*2,1:NUM.NUMERICS.no_nodes_ele*2)  = KM;
    LM(1:NUM.NUMERICS.no_nodes_ele*2,NUM.NUMERICS.no_nodes_ele*2+1:NUM.NUMERICS.no_nodes_ele*2+NUM.NUMERICS.no_nodes_ele_linear) = GM;
    LM(NUM.NUMERICS.no_nodes_ele*2+1:NUM.NUMERICS.no_nodes_ele*2+NUM.NUMERICS.no_nodes_ele_linear,1:NUM.NUMERICS.no_nodes_ele*2) = GM';
    
    vec(:,i) = LM(:);
    
    
    % add global right hand side vector
    NUM.Solve.FG(NUM.Number.number_ele_dof(:,i),1) = NUM.Solve.FG(NUM.Number.number_ele_dof(:,i),1) + F;
    
end
 
indx_j = repmat(1:(NUM.NUMERICS.no_nodes_ele*(NUM.NUMERICS.ndof-1)+NUM.NUMERICS.no_nodes_ele_linear),(NUM.NUMERICS.no_nodes_ele*(NUM.NUMERICS.ndof-1)+NUM.NUMERICS.no_nodes_ele_linear),1); 
indx_i = indx_j';
Ai = NUM.Number.number_ele_dof(indx_i(:),:);
Aj = NUM.Number.number_ele_dof(indx_j(:),:);
 
NUM.Solve.L = sparse(Ai,Aj,vec);

clear Ai Aj vec

% Set penalty on the Picard matrix
NUM.Solve.L(NUM.NUMERICS.no_nodes*(NUM.NUMERICS.ndof-1)+1:end,NUM.NUMERICS.no_nodes*(NUM.NUMERICS.ndof-1)+1:end) =  NUM.Solve.L(NUM.NUMERICS.no_nodes*(NUM.NUMERICS.ndof-1)+1:end,NUM.NUMERICS.no_nodes*(NUM.NUMERICS.ndof-1)+1:end)  +  ((-1/(1e10*100)).*speye(NUM.NUMERICS.no_nodes_linear,NUM.NUMERICS.no_nodes_linear));

% Add boundaries to the system
for i = 1:1:length(NUM.Boundary.bcdof)
    
    NUM.Solve.L(NUM.Boundary.bcdof(i),:) = 0;
    NUM.Solve.L(NUM.Boundary.bcdof(i),NUM.Boundary.bcdof(i)) = 1;
    NUM.Solve.FG(NUM.Boundary.bcdof(i)) = NUM.Boundary.bcval(i);
    
end

end

