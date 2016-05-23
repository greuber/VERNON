function [ MESH,PAR,NUM ] = Stokes_solve( MESH,NUM,PAR )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


NUM.r = zeros(NUM.no_nodes*(NUM.ndof-1)+NUM.no_nodes_linear,1);
NUM.Solve.f_res_brutal = zeros(NUM.no_nodes*(NUM.ndof-1)+NUM.no_nodes_linear,1);
for iter = 1:100;   % say we make a maximum of 100 iterations
    %                             NUM.r = NUM.r';
    %                             NUM.Result.Ux_vec = NUM.r(NUM.Number.number_dof(1,1:NUM.no_nodes));
    %                             NUM.Result.Uy_vec = NUM.r(NUM.Number.number_dof(2,1:NUM.no_nodes));
    %                             MESH.GCOORD(1,:) = MESH.Old(1,:) + NUM.Result.Ux_vec * PAR.dt;
    %                             MESH.GCOORD(2,:) = MESH.Old(2,:) + NUM.Result.Uy_vec * PAR.dt;
    %                             NUM.r = NUM.r';
    
    NUM.time_solver_iter = cputime;
    [ NUM ] = get_globals_picard( NUM,PAR,MESH);
    
    % ensure that the stiffness matrices are not zero
    % on the diagonal
    NUM.Solve.L = NUM.Solve.L + realmin*eye(size(NUM.Solve.L));
    
    
    NUM.Solve.f_res = NUM.Solve.L*NUM.r - NUM.Solve.FG;
    [ NUM ] = get_globals_Jacobian_analytical( NUM,PAR,MESH );
    
    % ensure that the stiffness matrices are not zero
    % on the diagonal
    NUM.Solve.J = NUM.Solve.J + realmin*eye(size(NUM.Solve.J));
    
    dr = NUM.Solve.J\(-NUM.Solve.f_res);
    
    NUM.r = NUM.r + dr;
    
    % update newton iteration counter
    NUM.time_solver_iter = cputime - NUM.time_solver_iter ;
    display(sprintf('norm residual newton = %6.9f',norm(NUM.Solve.f_res)));
    
    if norm(NUM.Solve.f_res)<PAR.tol
        break
    end
    
    
    
end

