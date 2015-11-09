function [ NUM,f_quad_val,f_line_val,MESH ] = get_element_res( p,u,PAR ,MESH,NUM,i,CHAR,f_quad_val,f_line_val)
% Evaluates the residual element-wise
fac        = 0;
x          = MESH.GCOORD(:,NUM.Number.number_quad(:,i)) +  u(NUM.Number.number_dof(:,1:NUM.NUMERICS.no_nodes_ele)) * PAR.dt *fac;    %%% CHANGED *0

B = zeros(3,NUM.NUMERICS.no_nodes_ele*2);
m = [1 1 0]';
N_matrix = zeros(2,NUM.NUMERICS.no_nodes_ele*2);

for j =1:NUM.NUMERICS.no_intp
    
    N_p      = NUM.Shape.pNI_vector{j};
    NI       = NUM.Shape.NI_vector{j};
    N_matrix(1,1:2:end) = NI(1,:);
    N_matrix(2,2:2:end) = NI(1,:);
    
    dNds = NUM.Shape.dNI_vector{j};
    Jac  = dNds*x';
    
    % Inverse of Jacobian
    invJ = inv(Jac);
    detJ = det(Jac);
    dNdX = invJ*dNds;
    
    B(1,1:2:end) = dNdX(1,:);
    B(2,2:2:end) = dNdX(2,:);
    B(3,1:2:end) = dNdX(2,:);
    B(3,2:2:end) = dNdX(1,:);
    
    % [ NUM,MESH] = ComputeStress_residual(i,NUM,j,MESH,CHAR,PAR,B,p,u);
    
    strain_tensor = [NUM.Strain.Exx(j,i);NUM.Strain.Ezz(j,i);NUM.Strain.Exz(j,i)*2];
    tau_local   =   [NUM.Stress.Txx(j,i);NUM.Stress.Tzz(j,i);NUM.Stress.Txz(j,i)];
    
    res_P       =   N_p'*(m'*strain_tensor)* MESH.INTP.weight(j) * detJ;
    % f_quad_val  =   f_quad_val + (B'*(tau_local - m*N_p*p) - N_matrix'*MESH.CompVar.rho(NUM.Number.number_quad(j,i))' * PAR.g * (MESH.INTP.weight(j) * detJ));
    f_quad_val  =   f_quad_val + ((B'*(tau_local - m*N_p*p) - N_matrix'*MESH.CompVar.rho(NUM.Number.number_quad(j,i))' * PAR.g) * (MESH.INTP.weight(j) * detJ));
    f_line_val  =   f_line_val - res_P;
    
end
end

