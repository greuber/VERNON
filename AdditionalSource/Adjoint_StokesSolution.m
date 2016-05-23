function [ NUM,MESH ] = Adjoint_StokesSolution( MESH,PAR,NUM,CHAR )
% Computes the Stokes solution for the adjoint inversion

atol                 = 1e-15;
NUM.Solve.number_pic = 1;
alpha                = 1;
res_ini              = 1;

NUM.Solve.r = zeros(NUM.NUMERICS.no_nodes*(NUM.NUMERICS.ndof-1)+NUM.NUMERICS.no_nodes_linear,1);

% apply bounds on the solution vector
NUM.Solve.r(NUM.Boundary.bcdof) = NUM.Boundary.bcval;

for iter = 1:2;   % say we make a maximum of 2 iterations
    NUM.time_solver_iter = cputime;
    [ NUM,MESH ]     = get_globals_Jacobian_analytical( NUM,PAR,MESH ,CHAR);
    [ NUM,MESH ]     = Compute_elemental_residual( NUM,MESH,PAR,CHAR );
    dr               = NUM.Solve.J\(-NUM.Solve.f_res);
    
    
    switch NUM.Solve.Linesearch
        case 'Bisection'
            % Bisection line-search
            alpha = 1;
            r_0   = NUM.Solve.r;
            
            while norm(N_f_res_old) < norm(NUM.Solve.f_res) && LS <= 10 && NUM.Solve.number_pic>1
                alpha           = alpha/2;
                NUM.Solve.r     = r_0 + alpha * dr;
                [ NUM ]          = Compute_elemental_residual( NUM,MESH,PAR,CHAR );
                LS              = LS + 1;
            end
            N_f_res_old     = norm(NUM.Solve.f_res);
            LS              = 1;
            
            
        case 'Jeremic'
            % Jeremic line-search
            r_0      = NUM.Solve.r;
            res_0    = NUM.Solve.f_res;
            alpha           = 1;
            omega           = 0.25;
            temp            = NUM.Solve.J*dr;
            dF              = temp;
            F_new           = res_0 + alpha*omega * dF;
            
            NUM.Solve.r     = r_0 + alpha*dr;
            [ NUM ]          = Compute_elemental_residual( NUM,MESH,PAR,CHAR );
            
            while norm(NUM.Solve.f_res) > norm(F_new)
                alpha           = alpha/2;
                dF              = temp;
                F_new           = res_0 + alpha*omega * dF;
                NUM.Solve.r     = r_0 + alpha*dr;
                [ NUM ]          = Compute_elemental_residual( NUM,MESH,PAR,CHAR );
                
                if alpha < 0.01
                    break;
                end
            end
            
        case 'none'
            alpha = 1;
            r_0 = NUM.Solve.r;
    end
    
    NUM.Solve.r      = r_0 + alpha * dr;
    
    NUM.time_solver_iter = cputime - NUM.time_solver_iter ;
    display(sprintf('norm residual newton = %6.3e; alpha = %4.4f',(norm(NUM.Solve.f_res)/(res_ini)),alpha));
    
    if NUM.Solve.number_pic == 1
        res_ini = norm(NUM.Solve.f_res);
    end
    
    if norm(NUM.Solve.f_res)<1e-15
        display(sprintf('Convergence due to machine precision; norm = %6.3e',norm(NUM.Solve.f_res)))
        break
    end
    if norm(NUM.Solve.f_res) < res_ini * PAR.tol + atol
        break
    end
    
    NUM.Solve.number_pic = NUM.Solve.number_pic+1;
end
end

