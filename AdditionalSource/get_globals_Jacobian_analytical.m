function [ NUM,MESH ] = get_globals_Jacobian_analytical( NUM,PAR,MESH ,CHAR)

%% initilize globals
vec = zeros((NUM.NUMERICS.no_nodes_ele*(NUM.NUMERICS.ndof-1)+NUM.NUMERICS.no_nodes_ele_linear)^2,NUM.NUMERICS.no_elems_global);
NUM.Solve.J = zeros(NUM.NUMERICS.no_nodes_ele*(NUM.NUMERICS.ndof-1)+NUM.NUMERICS.no_nodes_ele_linear,NUM.NUMERICS.no_nodes_ele*(NUM.NUMERICS.ndof-1)+NUM.NUMERICS.no_nodes_ele_linear);

% Some ifs
if NUM.Plasticity.Plasticity
    NUM.Plasticity.Plastic = zeros(NUM.NUMERICS.no_intp,NUM.NUMERICS.no_elems_global)   ;
end

I = eye(4,4);
B = zeros(3,NUM.NUMERICS.no_nodes_ele*2);
m = [1 1 0]';
P = (1/2)*[ 2 , 0 , 0 , 0
            0 , 2 , 0 , 0
            0 , 0 , 1 , 1];


for i = 1:NUM.NUMERICS.no_elems_global
    
    % initialize local matrices
    KM = zeros(NUM.NUMERICS.no_nodes_ele*(NUM.NUMERICS.ndof-1),NUM.NUMERICS.no_nodes_ele*(NUM.NUMERICS.ndof-1));
    GM = zeros(NUM.NUMERICS.no_nodes_ele*(NUM.NUMERICS.ndof-1),NUM.NUMERICS.no_nodes_ele_linear);
    VP = zeros(NUM.NUMERICS.no_nodes_ele*(NUM.NUMERICS.ndof-1),NUM.NUMERICS.no_nodes_ele_linear);
    LM = zeros(NUM.NUMERICS.no_nodes_ele*(NUM.NUMERICS.ndof-1)+NUM.NUMERICS.no_nodes_ele_linear,NUM.NUMERICS.no_nodes_ele*(NUM.NUMERICS.ndof-1)+NUM.NUMERICS.no_nodes_ele_linear);
    
    LGCOORD = (MESH.GCOORD(:,NUM.Number.number_quad(:,i)))';      % local global coordinates for this element
    
    for j =1:NUM.NUMERICS.no_intp
        
        N_p         = NUM.Shape.pNI_vector{j};
        
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
        
        if NUM.Plasticity.Plasticity
            sin_phi     = sin(deg2rad(MESH.CompVar.Phi(NUM.Number.number_quad(j,i))));
        end
        
        strain_tensor = [NUM.Strain.Exx(j,i);NUM.Strain.Ezz(j,i);NUM.Strain.Exz(j,i)*2];
        strain_tensor    = P' * strain_tensor;

        if NUM.Plasticity.Plasticity && NUM.Plasticity.Plastic(j,i) == 1 
            nu              =   strain_tensor./NUM.Strain.str_invariant_2d(j,i);
            D               =   P * (2*NUM.Viscosity.mu_2d(j,i) * (I - (1/2) * (nu*nu'))) * P';
            temp_PV         =   P*nu.*sin_phi-m;
            temp_PV         =   B' * temp_PV * N_p * (MESH.INTP.weight(j) * detJ);
            VP              =   VP + temp_PV; 
        else
            beta            =   (1/2)*(1/MESH.CompVar.powerlaw(NUM.Number.number_quad(j,i)) - 1);
            nu              =   strain_tensor./NUM.Strain.str_invariant_2d(j,i);
            D               =   P * (2 * NUM.Viscosity.mu_2d(j,i) * (I + beta*(nu*nu'))) * P';
            temp_P          =   B'*m*N_p * (MESH.INTP.weight(j) * detJ);
            VP              =   VP - temp_P;
        end
        
        temp     =   B'*D*B*(MESH.INTP.weight(j) * detJ);
        KM       =   KM + temp;
        temp_P   =   B'*m*N_p * (MESH.INTP.weight(j) * detJ);
        GM       =   GM - temp_P;
        
        
    end
    
    LM(1:NUM.NUMERICS.no_nodes_ele*2,1:NUM.NUMERICS.no_nodes_ele*2)  = KM;
    LM(1:NUM.NUMERICS.no_nodes_ele*2,NUM.NUMERICS.no_nodes_ele*2+1:NUM.NUMERICS.no_nodes_ele*2+NUM.NUMERICS.no_nodes_ele_linear) = VP;
    LM(NUM.NUMERICS.no_nodes_ele*2+1:NUM.NUMERICS.no_nodes_ele*2+NUM.NUMERICS.no_nodes_ele_linear,1:NUM.NUMERICS.no_nodes_ele*2) = GM';
    
    vec(:,i) = LM(:);
    
end

indx_j = repmat(1:(NUM.NUMERICS.no_nodes_ele*(NUM.NUMERICS.ndof-1)+NUM.NUMERICS.no_nodes_ele_linear),(NUM.NUMERICS.no_nodes_ele*(NUM.NUMERICS.ndof-1)+NUM.NUMERICS.no_nodes_ele_linear),1); 
indx_i = indx_j';
Ai = NUM.Number.number_ele_dof(indx_i(:),:);
Aj = NUM.Number.number_ele_dof(indx_j(:),:);
 
NUM.Solve.J = sparse(Ai,Aj,vec);

clear Ai Aj vec

% Set penalty on the Jacobian 
NUM.Solve.J(NUM.NUMERICS.no_nodes*(NUM.NUMERICS.ndof-1)+1:end,NUM.NUMERICS.no_nodes*(NUM.NUMERICS.ndof-1)+1:end) =  NUM.Solve.J(NUM.NUMERICS.no_nodes*(NUM.NUMERICS.ndof-1)+1:end,NUM.NUMERICS.no_nodes*(NUM.NUMERICS.ndof-1)+1:end)  +  ((-1/(1e10*100)).*speye(NUM.NUMERICS.no_nodes_linear,NUM.NUMERICS.no_nodes_linear));


% Add boundaries to the system
for i = 1:1:length(NUM.Boundary.bcdof)
    NUM.Solve.J(NUM.Boundary.bcdof(i),:) = 0;
    % NUM.Solve.J(:,NUM.Boundary.bcdof(i)) = 0;
    NUM.Solve.J(NUM.Boundary.bcdof(i),NUM.Boundary.bcdof(i)) = 1; 
end


end

