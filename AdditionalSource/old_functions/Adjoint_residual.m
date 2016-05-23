function [ drdx ] = Adjoint_residual( NUM,PAR,MESH,h )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

MESH.CompVar.rho_ini = MESH.CompVar.rho;
for id1 = 1:2
    NUM.Solve.f_res_brutal = zeros(NUM.no_nodes*(NUM.ndof-1)+NUM.no_nodes_linear,1);
    if id1 == 1
        MESH.CompVar.rho(unique(NUM.O_ind)) = MESH.CompVar.rho(unique(NUM.O_ind))+h;
    else
        MESH.CompVar.rho(unique(NUM.O_ind)) = MESH.CompVar.rho_ini(unique(NUM.O_ind));
    end
    
    for i = 1:NUM.O_elem     % compute residual
        NUM.f_quad_val = zeros(18,1);
        NUM.f_line_val = zeros(4,1);
        u = NUM.r(NUM.Number.number_ele_dof(1:18,i),1);
        p = NUM.r(NUM.Number.number_ele_dof(19:22,i),1);
        x = MESH.Old(:,NUM.Number.number_quad(:,i));
        [ NUM ] = get_element_res( p,u,x,PAR ,MESH,NUM,i);
        NUM.Solve.f_res_brutal(NUM.Number.number_ele_dof(1:18,i),1) = NUM.Solve.f_res_brutal(NUM.Number.number_ele_dof(1:18,i),1) + NUM.f_quad_val;
        NUM.Solve.f_res_brutal(NUM.Number.number_ele_dof(19:22,i),1) = NUM.Solve.f_res_brutal(NUM.Number.number_ele_dof(19:22,i),1) + NUM.f_line_val;
    end
    % Add boundaries to the system
    for i = 1:1:length(NUM.Boundary.bcdof)
        NUM.Solve.f_res_brutal(NUM.Boundary.bcdof(i)) = 0;
    end
    
    if id1 == 1
        res_plus = NUM.Solve.f_res_brutal;
    else
        res_normal = NUM.Solve.f_res_brutal;
    end
end
drdx = (res_plus - res_normal)./(h);

end

