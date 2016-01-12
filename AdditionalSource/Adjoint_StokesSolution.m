function [ NUM,MESH ] = Adjoint_StokesSolution( MESH,PAR,NUM,CHAR )
% Computes the Stokes solution for the adjoint inversion

atol = 1e-10;
NUM.Solve.number_pic = 1;
alpha = 1;
res_ini              = 1;

NUM.Solve.r = zeros(NUM.NUMERICS.no_nodes*(NUM.NUMERICS.ndof-1)+NUM.NUMERICS.no_nodes_linear,1);

% apply bounds on the solution vector
for i = 1:1:length(NUM.Boundary.bcdof)
    NUM.Solve.r(NUM.Boundary.bcdof(i)) = NUM.Boundary.bcval(i);
end

for iter = 1:30;   % say we make a maximum of 30 iterations
    NUM.time_solver_iter = cputime;
    [ NUM,MESH ]     = get_globals_Jacobian_analytical( NUM,PAR,MESH ,CHAR);
    [ NUM,MESH ]     = Compute_elemental_residual( NUM,MESH,PAR,CHAR );  
    dr               = NUM.Solve.J\(-NUM.Solve.f_res);
    NUM.Solve.r      = NUM.Solve.r + alpha * dr;
    
    NUM.time_solver_iter = cputime - NUM.time_solver_iter ;
    display(sprintf('norm residual newton = %6.3e; alpha = %4.4f',(norm(NUM.Solve.f_res)/(res_ini)),alpha));
    
    if norm(NUM.Solve.f_res) < res_ini * PAR.tol + atol
        break
    end
    
    if NUM.Solve.number_pic == 1
        res_ini = norm(NUM.Solve.f_res);
    end
    
    NUM.Solve.number_pic = NUM.Solve.number_pic+1;
end
end

