function [ NUM,MESH ] = Compute_elemental_residual( NUM,MESH,PAR,CHAR )
% This function computes the residual element wise

NUM.Solve.f_res = zeros(NUM.NUMERICS.no_nodes*(NUM.NUMERICS.ndof-1)+NUM.NUMERICS.no_nodes_linear,1);
% Some ifs
if NUM.Plasticity.Plasticity
    NUM.Plasticity.Plastic = zeros(NUM.NUMERICS.no_intp,NUM.NUMERICS.no_elems_global)   ;
end

% % This would be the old school way to do it
% [ NUM,MESH ]     = get_globals_picard( NUM,PAR,MESH,CHAR);
% NUM.Solve.f_res  = NUM.Solve.L*NUM.Solve.r - NUM.Solve.FG;

for i = 1:NUM.NUMERICS.no_elems_global
    f_quad_val = zeros(NUM.NUMERICS.no_nodes_ele*2,1);
    f_line_val = zeros(NUM.NUMERICS.no_nodes_ele_linear,1);
    u = NUM.Solve.r(NUM.Number.number_ele_dof(1:NUM.NUMERICS.no_nodes_ele*2,i),1);
    p = NUM.Solve.r(NUM.Number.number_ele_dof(NUM.NUMERICS.no_nodes_ele*2+1:NUM.NUMERICS.no_nodes_ele*2+NUM.NUMERICS.no_nodes_ele_linear,i),1);
    [ NUM,f_quad_val,f_line_val,MESH ] = get_element_res( p,u,PAR ,MESH,NUM,i,CHAR,f_quad_val,f_line_val);
    NUM.Solve.f_res(NUM.Number.number_ele_dof(1:NUM.NUMERICS.no_nodes_ele*2,i),1) = NUM.Solve.f_res(NUM.Number.number_ele_dof(1:NUM.NUMERICS.no_nodes_ele*2,i),1) + f_quad_val;
    NUM.Solve.f_res(NUM.Number.number_ele_dof(NUM.NUMERICS.no_nodes_ele*2+1:NUM.NUMERICS.no_nodes_ele*2+NUM.NUMERICS.no_nodes_ele_linear,i),1) = NUM.Solve.f_res(NUM.Number.number_ele_dof(NUM.NUMERICS.no_nodes_ele*2+1:NUM.NUMERICS.no_nodes_ele*2+NUM.NUMERICS.no_nodes_ele_linear,i),1) + f_line_val;
end
for i = 1:1:length(NUM.Boundary.bcdof)
    NUM.Solve.f_res(NUM.Boundary.bcdof(i)) = 0;
end

end

