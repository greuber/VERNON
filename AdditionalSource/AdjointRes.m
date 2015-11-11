function [ drdx_temp ] = AdjointRes( NUM,PAR,MESH,h,input,index,field,par,CHAR)
% Computes the derivative of the residual versus design variable with a FD
% approximation

input_ini = MESH.CompVar.(field);

for id1 = 1:2
    NUM.Solve.f_res_brutal = zeros(NUM.NUMERICS.no_nodes*(NUM.NUMERICS.ndof-1)+NUM.NUMERICS.no_nodes_linear,1);
    if id1 == 1
        MESH.CompVar.(field)(unique(index)) = MESH.CompVar.(field)(unique(index))+h;
    else
        MESH.CompVar.(field)(unique(index)) = input_ini(unique(index));
    end
    
    for i = 1:NUM.NUMERICS.no_elems_global     % compute residual
        f_quad_val = zeros(NUM.NUMERICS.no_nodes_ele*2,1);
        f_line_val = zeros(NUM.NUMERICS.no_nodes_ele_linear,1);
        u = NUM.Solve.r(NUM.Number.number_ele_dof(1:NUM.NUMERICS.no_nodes_ele*2,i),1);
        p = NUM.Solve.r(NUM.Number.number_ele_dof(NUM.NUMERICS.no_nodes_ele*2+1:NUM.NUMERICS.no_nodes_ele*2+NUM.NUMERICS.no_nodes_ele_linear,i),1);
        [ NUM,f_quad_val,f_line_val,MESH ] = get_element_res( p,u,PAR ,MESH,NUM,i,CHAR,f_quad_val,f_line_val);
        NUM.Solve.f_res_brutal(NUM.Number.number_ele_dof(1:NUM.NUMERICS.no_nodes_ele*2,i),1) = NUM.Solve.f_res_brutal(NUM.Number.number_ele_dof(1:NUM.NUMERICS.no_nodes_ele*2,i),1) + f_quad_val;
        NUM.Solve.f_res_brutal(NUM.Number.number_ele_dof(NUM.NUMERICS.no_nodes_ele*2+1:NUM.NUMERICS.no_nodes_ele*2+NUM.NUMERICS.no_nodes_ele_linear,i),1) = NUM.Solve.f_res_brutal(NUM.Number.number_ele_dof(NUM.NUMERICS.no_nodes_ele*2+1:NUM.NUMERICS.no_nodes_ele*2+NUM.NUMERICS.no_nodes_ele_linear,i),1) + f_line_val;
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
drdx_temp = (res_plus - res_normal)./(h);

end

