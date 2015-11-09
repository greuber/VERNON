function [ NUM,MESH] = Solver(NUM,PAR,MESH,CHAR )

NUM.Solve.number_pic = 1;
max_iter             = 10;
N_f_res_old          = realmax;
LS                   = 1;
no_picard            = 3;
NUM.Solve.Solver     = NUM.Solve.Solver_method;
res_ini              = 1;
number_new           = 1;

atol = 1e-10;
tol_picard = 1e-2;   % when to switch to Newton iterations

% apply bounds on the solution vector
for i = 1:1:length(NUM.Boundary.bcdof)
    NUM.Solve.r(NUM.Boundary.bcdof(i)) = NUM.Boundary.bcval(i);
end



% get all shape function combinations at once
[NUM ] = ComputeShapeFunctionsVector( MESH.INTP.COORD,NUM );

while NUM.Solve.number_pic<=max_iter
    
    % in case of the joint system the code has to decide the solver
    if strcmp(NUM.Solve.Solver_method,'NewtonPicard') == 1
        if NUM.Solve.number_pic <= no_picard
            NUM.Solve.Solver = 'Picard';
        elseif norm(NUM.Solve.f_res) < tol_picard
            NUM.Solve.Solver        = 'Newton';
            number_new    = 1;
        else
            NUM.Solve.Solver = 'Newton';
            number_new    = 1;
        end
    end
    
    switch NUM.Solve.Solver
        
        case 'Picard'
            if NUM.Solve.BrutalForce == 1
                [ NUM,MESH ]           = Compute_elemental_residual( NUM,MESH,PAR,CHAR );
                NUM.Solve.f_res_brutal = NUM.Solve.f_res; 
                display(sprintf('norm residual brutal = %6.3e',norm(NUM.Solve.f_res_brutal/res_ini)));
                NUM.Solve.f_res        = NUM.Solve.f_res_brutal;
            end
            
            NUM.time_solver_iter = cputime;
            [ NUM,MESH ]    = get_globals_picard( NUM,PAR,MESH,CHAR);
            NUM.Solve.f_res = NUM.Solve.L*NUM.Solve.r - NUM.Solve.FG;
            NUM.Solve.r     = NUM.Solve.L\NUM.Solve.FG;
            NUM.time_solver_iter = cputime - NUM.time_solver_iter;
            

            if strcmp(NUM.Timestep.method,'EulerImplicit') == 1
                MESH.GCOORD(1,:) = MESH.Old(1,:) + NUM.Solve.r(NUM.Number.number_dof(1,1:NUM.NUMERICS.no_nodes))' * PAR.dt;
                MESH.GCOORD(2,:) = MESH.Old(2,:) + NUM.Solve.r(NUM.Number.number_dof(2,1:NUM.NUMERICS.no_nodes))' * PAR.dt;
            end
            
            display(sprintf('norm residual picard = %6.3e',(norm(NUM.Solve.f_res)/(res_ini))));
            
        case 'Newton'
            if NUM.Solve.BrutalForce == 1
                [ NUM ] = Compute_elemental_residual( NUM,MESH,PAR,CHAR );
                NUM.Solve.f_res_brutal = NUM.Solve.f_res; 
                display(sprintf('norm residual brutal = %6.3e',norm(NUM.Solve.f_res_brutal/res_ini)));
                NUM.Solve.f_res = NUM.Solve.f_res_brutal;
            end
            
            switch NUM.Solve.Jacobian
                case 'analytical'
                    NUM.time_solver_iter = cputime;

                    [ NUM,MESH ]     = get_globals_picard( NUM,PAR,MESH,CHAR);
%                     NUM.Solve.f_res  = NUM.Solve.L*NUM.Solve.r - NUM.Solve.FG;
%                     [ NUM,MESH ]     = get_globals_Jacobian_analytical( NUM,PAR,MESH ,CHAR);
                    [ NUM,MESH ]     = Compute_elemental_residual( NUM,MESH,PAR,CHAR );

                    dr               = NUM.Solve.L\(-NUM.Solve.f_res);                              %%% CHANGED PICARD
                    
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
                    
                    if strcmp(NUM.Timestep.method,'EulerImplicit') == 1
                        MESH.GCOORD(1,:) = MESH.Old(1,:) + NUM.Solve.r(NUM.Number.number_dof(1,1:NUM.NUMERICS.no_nodes))' * PAR.dt;
                        MESH.GCOORD(2,:) = MESH.Old(2,:) + NUM.Solve.r(NUM.Number.number_dof(2,1:NUM.NUMERICS.no_nodes))' * PAR.dt;
                    end
                    
                    display(sprintf('norm residual newton = %6.3e; alpha = %4.4f',(norm(NUM.Solve.f_res)/(res_ini)),alpha));
                    number_new = number_new + 1;                            % iteration count for the newton iterations
                    
                case 'adjoint'
                    grad = ones(size(NUM.NUMERICS.no_nodes*(NUM.NUMERICS.ndof-1)+NUM.NUMERICS.no_nodes_linear,1))*realmax;
                    NUM.Adjoint.adjoint_tol = 1e-15;    % to reset an initial very small tolerance (is computed in the code automatically later)
                    
                    npar        = length(NUM.Adjoint.m);
                    
                    % Start loop for adjoint convergence
                    while norm(grad) > NUM.Adjoint.adjoint_tol
                        
                        if NUM.Adjoint.it_adj == 1   % first iteration is special for the gradient and Hessian
                            
                            % Compute Stokes solution
                            [ NUM,MESH ] = Adjoint_StokesSolution( MESH,PAR,NUM,CHAR );
                            
                            % [NUM.Adjoint.m] = denondimensionalize(NUM,NUM.Adjoint.m,npar,CHAR);
                            [NUM.Adjoint.m] = normalize(NUM,NUM.Adjoint.m,npar);
                            
                            % Compute cost function and derivative
                            fcost = (1/2) * (NUM.Solve.r(NUM.Adjoint.ind_cost) - NUM.Solve.r_ini(NUM.Adjoint.ind_cost))' * (NUM.Solve.r(NUM.Adjoint.ind_cost) - NUM.Solve.r_ini(NUM.Adjoint.ind_cost));
                            dcost_du = zeros(size(NUM.Solve.r'));
                            dcost_du(NUM.Adjoint.ind_cost) = (NUM.Solve.r(NUM.Adjoint.ind_cost)' - NUM.Solve.r_ini(NUM.Adjoint.ind_cost)');
                            
                            [NUM.Adjoint.m] = denormalize(NUM,NUM.Adjoint.m,npar);
                            [NUM.Adjoint.m] = nondimensionalize(NUM,NUM.Adjoint.m,npar,CHAR);
                            for par = 1:length(NUM.Adjoint.fields)   % loop over design variables
                                h_max = 1e-28;
                                h = max([(1e-6*abs(NUM.Adjoint.m(par,1))), h_max]);
                                input = MESH.CompVar.(NUM.Adjoint.fields{par});
                                index = NUM.Adjoint.index{par};
                                field = NUM.Adjoint.fields{par};    % the perturbed field in the MESH.CompVar structure
                                % Compute residual by forward finite differences
                                [ drdx_temp ] = AdjointRes( NUM,PAR,MESH,h,input,index,field ,par,CHAR);
                                drdx(:,par) = drdx_temp;
                            end
                            [NUM.Adjoint.m] = denondimensionalize(NUM,NUM.Adjoint.m,npar,CHAR);
                            [NUM.Adjoint.m] = normalize(NUM,NUM.Adjoint.m,npar);
                            
                            % compute gradient using adjoint method
                            psi = (NUM.Solve.J')\(dcost_du');
                            grad = - psi'*drdx;
                            [grad] = nondimensionalize(NUM,grad',npar,CHAR);
                            grad = grad';
                            [grad] = normalize(NUM,grad',npar);
                            grad = grad';
                            
                            % In the first iteration the Hessian is the
                            % identity matrix
                            H = eye(length(NUM.Adjoint.fields),length(NUM.Adjoint.fields));     % ATTENTION: THIS CAN BE NORMALIZED BY H*(1/norm(grad))
                            
                            % setup the tolerance
                            NUM.Adjoint.adjoint_tol = norm(grad)*1e-3;   % to end when the error is 3 orders smaller than the in iitla one
                            
                            % save step-size
                            beta_step = (1./abs(grad'));   
                            
                            % Compute first search-direction
                            dx = -H*grad';
                            
                        else      % every iteration starting from the second                          
                            % Update variables for the line search (needs
                            % the 'initial' values to compare)
                            fcost_ini = fcost;
                            grad_ini = grad;
                            sol_ini = NUM.Adjoint.m;
                            ij = 2;   % iteration counter of the line-search
                            
                            % beta_low is used to compare cost functions
                            % to convergence so in the frist line-search we
                            % assume it to be better than no change (0)
                            
                            % save step-size
                            beta_step = (1./abs(dx))./NUM.Adjoint.LS_Parameter;
                            
                            for par = 1:npar
                                if (sol_ini(par) + beta_step(par)*dx(par) >= NUM.Adjoint.bounds{par}(2))
                                    beta_step(par) = beta_step(par)/dx(par);
                                end
                            end
                            
                            
                            beta_step_max = 1;
                            beta_step_old = 0;
                            
                            % Starting the line-search
                            while 1>0
                                % Compute fcost with an initial guess scaled
                                % with the former beta_step
                                NUM.Adjoint.m = sol_ini + beta_step.*dx;
                                
                                [NUM.Adjoint.m] = denormalize(NUM,NUM.Adjoint.m,npar);
                                [NUM.Adjoint.m] = nondimensionalize(NUM,NUM.Adjoint.m,npar,CHAR);
                                for par = 1:length(NUM.Adjoint.fields)
                                    MESH.CompVar.(NUM.Adjoint.fields{par})(NUM.Adjoint.index{par}) = NUM.Adjoint.m(par,1);
                                end
                                [NUM.Adjoint.m] = denondimensionalize(NUM,NUM.Adjoint.m,npar,CHAR);
                                [NUM.Adjoint.m] = normalize(NUM,NUM.Adjoint.m,npar);
                                
                                [ NUM,MESH ] = Adjoint_StokesSolution( MESH,PAR,NUM,CHAR );
                                fcost = (1/2) * (NUM.Solve.r(NUM.Adjoint.ind_cost) - NUM.Solve.r_ini(NUM.Adjoint.ind_cost))' * (NUM.Solve.r(NUM.Adjoint.ind_cost) - NUM.Solve.r_ini(NUM.Adjoint.ind_cost));
                                
                                % First condition ('sufficient decrease
                                % condition')
                                if (fcost > fcost_ini + 1e-4 * beta_step .* grad_ini') | (fcost >= fcost_old)
                                    [beta_step_star] = zoom_tot( beta_step_old,beta_step,NUM,MESH,PAR,sol_ini,dx,fcost_ini,grad_ini,dcost_du,fcost,CHAR);
                                    beta_step_star = beta_step_star;
                                    break
                                end
                                
                                % Compute derivative of cost and gradient
                                % with the current step-size
                                dcost_du = zeros(size(NUM.Solve.r'));
                                dcost_du(NUM.Adjoint.ind_cost) = (NUM.Solve.r(NUM.Adjoint.ind_cost)' - NUM.Solve.r_ini(NUM.Adjoint.ind_cost)');
                                
                                [NUM.Adjoint.m] = denormalize(NUM,NUM.Adjoint.m,npar);
                                [NUM.Adjoint.m] = nondimensionalize(NUM,NUM.Adjoint.m,npar,CHAR);
                                for par = 1:length(NUM.Adjoint.fields)   % loop over design variables
                                    h_max = 1e-28;
                                    h = max([(1e-6*abs(NUM.Adjoint.m(par,1))), h_max]);
                                    input = MESH.CompVar.(NUM.Adjoint.fields{par});
                                    index = NUM.Adjoint.index{par};
                                    field = NUM.Adjoint.fields{par};    % the perturbed field in the MESH.CompVar structure
                                    % Compute residual by forward finite differences
                                    [ drdx_temp ] = AdjointRes( NUM,PAR,MESH,h,input,index,field,par ,CHAR);
                                    drdx(:,par) = drdx_temp;
                                end
                                [NUM.Adjoint.m] = denondimensionalize(NUM,NUM.Adjoint.m,npar,CHAR);
                                [NUM.Adjoint.m] = normalize(NUM,NUM.Adjoint.m,npar);
                                
                                % Compute gradient using adjoint method
                                psi = (NUM.Solve.J')\(dcost_du');
                                grad = - psi'*drdx;
                                [grad] = nondimensionalize(NUM,grad',npar,CHAR);
                                grad = grad';
                                [grad] = normalize(NUM,grad',npar);
                                grad = grad';
                                
                                % Second condition ('curvature condition')
                                if abs(grad) <= -0.9*grad_ini   % CHANGED       % BUT GRADIENT IS STILL starting OSCILLATING          % MAYBE COULD BE any(abs(grad) <= abs(-0.9*grad_ini))
                                    beta_step_star = beta_step;
                                    break
                                end
%                                 if grad >= 0        % CHANGED       % BUT GRADIENT IS STILL starting OSCILLATING          % MAYBE COULD BE any(abs(grad) <= abs(-0.9*grad_ini))
%                                     [beta_step_star] = zoom_tot( beta_step,beta_step_old,NUM,MESH,PAR,sol_ini,dx,fcost_ini,grad_ini,dcost_du,fcost,CHAR);
%                                     break
%                                 end
                                
                                % if non of the conditions is fullfilled we
                                % can increase beta_step (can be done cubic
                                % but here just double the value) and
                                % update all old variables
                                
                                %                                 % cubic interpolation
                                %                                 r1 = dcost_du_old + dcost_du - 3*((fcost_old - fcost)/(beta_step_old - beta_step));
                                %                                 r2 = sqrt(r1.^2 - dcost_du_old .* dcost_du);
                                %                                 beta_step_new = mean(beta_step - (beta_step - beta_step_old) * ((dcost_du + r2 -r1)./(dcost_du - dcost_du_old + 2*r2)));
                                %
                                if beta_step >= (beta_step_max-1e-3)
                                    beta_step_star = beta_step_max;
                                    break
                                end
                                
                                %                                 beta_step_new = real(beta_step_new);
                                % beta_step_old = beta_step;
                                beta_step = beta_step + beta_step*1.1; % + (beta_step_max-beta_step)*rand(1);
                                fcost_old = fcost;
                                ij = ij +1
                                
                                beta_step
                                
                            end
                            beta_step = beta_step_star;                % update beta_step
                            
                            % Compute the gradient again if the difference
                            % is below zero (sometimes the gradient is
                            % simply not computed in the line-search)
                            if grad - grad_ini == 0
                                % Compute derivative of cost and gradient
                                % with the current step-size
                                dcost_du = zeros(size(NUM.Solve.r'));
                                dcost_du(NUM.Adjoint.ind_cost) = (NUM.Solve.r(NUM.Adjoint.ind_cost)' - NUM.Solve.r_ini(NUM.Adjoint.ind_cost)');
                                
                                [NUM.Adjoint.m] = denormalize(NUM,NUM.Adjoint.m,npar);
                                [NUM.Adjoint.m] = nondimensionalize(NUM,NUM.Adjoint.m,npar,CHAR);
                                for par = 1:length(NUM.Adjoint.fields)   % loop over design variables
                                    h_max = 1e-28;
                                    h = max([(1e-6*abs(NUM.Adjoint.m(par,1))), h_max]);
                                    input = MESH.CompVar.(NUM.Adjoint.fields{par});
                                    index = NUM.Adjoint.index{par};
                                    field = NUM.Adjoint.fields{par};    % the perturbed field in the MESH.CompVar structure
                                    % Compute residual by forward finite differences
                                    [ drdx_temp ] = AdjointRes( NUM,PAR,MESH,h,input,index,field,par,CHAR );
                                    drdx(:,par) = drdx_temp;
                                end
                                [NUM.Adjoint.m] = denondimensionalize(NUM,NUM.Adjoint.m,npar,CHAR);
                                [NUM.Adjoint.m] = normalize(NUM,NUM.Adjoint.m,npar);
                                
                                % Compute gradient using adjoint method
                                psi = (NUM.Solve.J')\(dcost_du');
                                grad = - psi'*drdx;
                                [grad] = nondimensionalize(NUM,grad',npar,CHAR);
                                grad = grad';
                                [grad] = normalize(NUM,grad',npar);
                                grad = grad';
                            end
                            
                            
                            NUM.Adjoint.m = sol_ini + beta_step.*dx;
                            [NUM.Adjoint.m] = denormalize(NUM,NUM.Adjoint.m,npar);
                            [NUM.Adjoint.m] = nondimensionalize(NUM,NUM.Adjoint.m,npar,CHAR);
                            for par = 1:length(NUM.Adjoint.fields)
                                MESH.CompVar.(NUM.Adjoint.fields{par})(NUM.Adjoint.index{par}) = NUM.Adjoint.m(par,1);
                            end
                            [NUM.Adjoint.m] = denondimensionalize(NUM,NUM.Adjoint.m,npar,CHAR);
                            [NUM.Adjoint.m] = normalize(NUM,NUM.Adjoint.m,npar);
                            
                            % Compute parameters for BFGS
                            dp = NUM.Adjoint.m - sol_ini;
                            vp = grad - grad_ini;
                            
                            %                             if NUM.Adjoint.it_adj == 2
                            %                                 % First initial guess of Hessian (approximates
                            %                             % an eigenvalue of the Hessian (can maybe be
                            %                             % left))
                            %                             H = (vp*dx)/(vp*vp');
                            %                             else
                            H = eye(length(NUM.Adjoint.fields),length(NUM.Adjoint.fields));
                            %                             end
                            
                            r = (dp./(dp.*vp')) - ((H*vp')/(vp*H*vp'));
                            
                            
                            
                            % Approximate inverse Hessian with BFGS
                            % algorithm
                            H = H - ((H*vp'*(H*vp')')/(vp*H*vp')) + ((dp*dp')/(dp'*vp')) + (vp*H*vp'*(r*r'));   % = H^-1
                            
                            % check for positive definiteness of H
                            [~,p] = chol(H);
                            if p>0
                                warning('Hessian no longer positive definite')
                            end
                            
                            % Compute new search-direction
                            dx = -H*grad';
                            
                        end
                        
                        
                        
                        % Plot solution (change in the design variable)
                        X = MESH.GCOORD(1,:);
                        X = X(NUM.Number.number_2d);
                        Z = MESH.GCOORD(2,:);
                        Z = Z(NUM.Number.number_2d);
                        X_ele = MESH.GCOORD(1,:);
                        X_ele = X_ele(NUM.Number.number_2d_linear);
                        Z_ele = MESH.GCOORD(2,:);
                        Z_ele = Z_ele(NUM.Number.number_2d_linear);
                        
                        figure(1),clf
                        for par = 1:length(NUM.Adjoint.fields)
                            design1_2d = MESH.CompVar.(NUM.Adjoint.fields{par})(NUM.Number.number_2d);
                            subplot(length(NUM.Adjoint.fields),1,par)
                            pcolor(X,Z,design1_2d)
                            colorbar
                            %                             if par == 1
                            %                                 caxis([1 2])
                            %                             elseif par ==2
                            %                                 caxis([1 10])
                            %                             elseif par == 3
                            %                                 caxis([5 10])
                            %                             elseif par ==4
                            %                                 caxis([1 2])
                            %                             end
                            title(['design variable ',num2str(par),' after ',num2str(NUM.Adjoint.it_adj),' iteration; Norm cost function = ',num2str(norm(fcost))])
                        end
                        drawnow
                        %                         fname = (['Iteration_jioint_3_',num2str(NUM.Adjoint.it_adj)]);
                        %                         print(fname,'-dpng','-zbuffer','-r300')
                        
                        % MESH.CompVar.powerlaw
                        MESH.CompVar.rho*CHAR.rho
                        
                        %                                                                         rho_sur_temp = mean(MESH.CompVar.rho(unique((NUM.Adjoint.index{3}))));
                        %                                                                         power_sur_temp = mean(MESH.CompVar.powerlaw(unique((NUM.Adjoint.index{1}))));
                        %                                                                         power_block_temp = mean(MESH.CompVar.powerlaw(unique((NUM.Adjoint.index{4}))));
                        %                                                                         mu_temp = mean(MESH.CompVar.mu_ref(unique((NUM.Adjoint.index{2}))));
                        %                                                                         save(['Adjoint_vectorized_4_',num2str(NUM.Adjoint.it_adj),'Rising_Sphere'],'power_sur_temp','mu_temp','rho_sur_temp','power_block_temp')
                        
                        display('------------------------------------------------------ ')
                        display(['NORM OF THE GRADIENT = ',num2str(norm(grad))])
                        display(['COST FUNCTION = ',num2str(fcost)])
                        display(['BETA STEP = ',num2str(beta_step')])
                        display('------------------------------------------------------ ')
                        
                        % update old variables
                        NUM.Adjoint.it_adj = NUM.Adjoint.it_adj + 1;
                        fcost_old = fcost;
                        dcost_du_old = dcost_du;
                        
                        % [ NUM ] = Update_StressesAndStrains(NUM,PAR,MESH );
                        
                    end
                    
                    
                    
                    
                    
                case 'numerical'
                    NUM.time_solver_iter = cputime;
                    vec = zeros((NUM.NUMERICS.no_nodes_ele*(NUM.NUMERICS.ndof-1)+NUM.NUMERICS.no_nodes_ele_linear)^2,NUM.NUMERICS.no_elems_global);
                    NUM.Solve.f_res = zeros(NUM.NUMERICS.no_nodes*(NUM.NUMERICS.ndof-1)+NUM.NUMERICS.no_nodes_linear,1);
                    if NUM.Plasticity.Plasticity
                        NUM.Plasticity.Plastic = zeros(NUM.NUMERICS.no_intp,NUM.NUMERICS.no_elems_global)   ;
                    end
                    for i = 1:NUM.NUMERICS.no_elems_global
                        f_quad_val = zeros(NUM.NUMERICS.no_nodes_ele*2,1);
                        f_line_val = zeros(NUM.NUMERICS.no_nodes_ele_linear,1);
                        J = zeros(NUM.NUMERICS.no_nodes_ele*(NUM.NUMERICS.ndof-1)+NUM.NUMERICS.no_nodes_ele_linear,NUM.NUMERICS.no_nodes_ele*(NUM.NUMERICS.ndof-1)+NUM.NUMERICS.no_nodes_ele_linear);
                        u = NUM.Solve.r(NUM.Number.number_ele_dof(1:NUM.NUMERICS.no_nodes_ele*2,i),1);
                        p = NUM.Solve.r(NUM.Number.number_ele_dof(NUM.NUMERICS.no_nodes_ele*2+1:NUM.NUMERICS.no_nodes_ele*2+NUM.NUMERICS.no_nodes_ele_linear,i),1);
                        [ NUM,f_quad_val,f_line_val,MESH ] = get_element_res( p,u,PAR ,MESH,NUM,i,CHAR,f_quad_val,f_line_val);
                        f_res_old(1:NUM.NUMERICS.no_nodes_ele*2,1)  = f_quad_val;
                        f_res_old(NUM.NUMERICS.no_nodes_ele*2+1:NUM.NUMERICS.no_nodes_ele*2+NUM.NUMERICS.no_nodes_ele_linear,1) = f_line_val;
                        
                        for j = 1:NUM.NUMERICS.no_nodes_ele*2
                            tmp = u(j);
                            u(j) = u(j) + PAR.perturb ;
                            [ NUM,f_quad_val,f_line_val,MESH ] = get_element_res( p,u,PAR ,MESH,NUM,i,CHAR,f_quad_val,f_line_val);
                            f_res_new(1:NUM.NUMERICS.no_nodes_ele*2,1)  = f_quad_val;
                            f_res_new(NUM.NUMERICS.no_nodes_ele*2+1:NUM.NUMERICS.no_nodes_ele*2+NUM.NUMERICS.no_nodes_ele_linear,1) = f_line_val;
                            J(:,j) = (f_res_new - f_res_old) /PAR.perturb;
                            u(j) = tmp ;
                            
                        end
                        
                        for j = 1:NUM.NUMERICS.no_nodes_ele_linear
                            tmp = p(j);
                            p(j) = p(j) + PAR.perturb ;
                            [ NUM,f_quad_val,f_line_val,MESH ] = get_element_res( p,u,PAR ,MESH,NUM,i,CHAR,f_quad_val,f_line_val);
                            f_res_new(1:NUM.NUMERICS.no_nodes_ele*2,1)  = f_quad_val;
                            f_res_new(NUM.NUMERICS.no_nodes_ele*2+1:NUM.NUMERICS.no_nodes_ele*2+NUM.NUMERICS.no_nodes_ele_linear,1) = f_line_val;
                            J(:,j+(NUM.NUMERICS.no_nodes_ele*2)) = (f_res_new - f_res_old) /PAR.perturb;
                            p(j) = tmp;
                            
                        end
                        vec(:,i) = J(:);
                        NUM.Solve.f_res(NUM.Number.number_ele_dof(:,i),1) = NUM.Solve.f_res(NUM.Number.number_ele_dof(:,i),1) + f_res_old;
                    end
                    
                    indx_j = repmat(1:(NUM.NUMERICS.no_nodes_ele*(NUM.NUMERICS.ndof-1)+NUM.NUMERICS.no_nodes_ele_linear),(NUM.NUMERICS.no_nodes_ele*(NUM.NUMERICS.ndof-1)+NUM.NUMERICS.no_nodes_ele_linear),1);
                    indx_i = indx_j';
                    Ai = NUM.Number.number_ele_dof(indx_i(:),:);
                    Aj = NUM.Number.number_ele_dof(indx_j(:),:);
                    
                    J_global = sparse(Ai,Aj,vec);
                    
                    clear Ai Aj vec
                    
                    % Set penalty on the Jacobian
                    J_global(NUM.NUMERICS.no_nodes*(NUM.NUMERICS.ndof-1)+1:end,NUM.NUMERICS.no_nodes*(NUM.NUMERICS.ndof-1)+1:end) =  J_global(NUM.NUMERICS.no_nodes*(NUM.NUMERICS.ndof-1)+1:end,NUM.NUMERICS.no_nodes*(NUM.NUMERICS.ndof-1)+1:end)  +  ((-1/(1e10*100)).*speye(NUM.NUMERICS.no_nodes_linear,NUM.NUMERICS.no_nodes_linear));
                    
                    for i = 1:1:length(NUM.Boundary.bcdof)
                        J_global(NUM.Boundary.bcdof(i),:) = 0;
                        J_global(NUM.Boundary.bcdof(i),NUM.Boundary.bcdof(i)) = 1;
                        NUM.Solve.f_res(NUM.Boundary.bcdof(i)) = 0;
                    end
                    
                    % for i = 1:1:length(NUM.Boundary.bcdof)
                    % NUM.Solve.r(NUM.Boundary.bcdof(i)) = NUM.Boundary.bcval(i);
                    % end
                    
                    dr    = J_global\(-NUM.Solve.f_res);
                    NUM.Solve.r = NUM.Solve.r + dr;
                    
                    if strcmp(NUM.Timestep.method,'EulerImplicit') == 1
                        MESH.GCOORD(1,:) = MESH.Old(1,:) + NUM.Solve.r(NUM.Number.number_dof(1,1:NUM.NUMERICS.no_nodes))' * PAR.dt;
                        MESH.GCOORD(2,:) = MESH.Old(2,:) + NUM.Solve.r(NUM.Number.number_dof(2,1:NUM.NUMERICS.no_nodes))' * PAR.dt;
                    end
                    
                    NUM.time_solver_iter = cputime - NUM.time_solver_iter;
                    display(sprintf('norm residual newton = %6.3e',(norm(NUM.Solve.f_res)/(res_ini))));
                    
            end
    end
    
    if NUM.Solve.number_pic == 1
        res_ini = norm(NUM.Solve.f_res);
    end
    
    if norm(NUM.Solve.f_res) < res_ini * PAR.tol + atol || NUM.Solve.number_pic >=max_iter
        % Stresses and Strains are updated in 'Timestepping'
        break
    end
    
    NUM.Solve.number_pic = NUM.Solve.number_pic+1;
    
end














% FUNCTIONS
% 1. normalize the variables
function [m] = normalize(NUM,m,npar)
for par = 1:npar
    m(par,1) = (m(par,1) - NUM.Adjoint.bounds{par}(1))/(NUM.Adjoint.bounds{par}(2) - NUM.Adjoint.bounds{par}(1));
end

% 2. denormalize the variables
function [m] = denormalize(NUM,m,npar)
for par = 1:npar
    m(par,1) = m(par,1) * (NUM.Adjoint.bounds{par}(2) - NUM.Adjoint.bounds{par}(1)) + NUM.Adjoint.bounds{par}(1);
end

% 3. nondimensionalize the variables
function [m] = nondimensionalize(NUM,m,npar,CHAR)
for par = 1:npar
    if strcmp(NUM.Adjoint.fields{par},'powerlaw') == 1
        m(par,1) = m(par,1);
    elseif strcmp(NUM.Adjoint.fields{par},'str_ref') == 1
        m(par,1) = m(par,1)/(1/CHAR.Time);
    elseif isfield(CHAR,NUM.Adjoint.fields{par}) == 1
        m(par,1) = m(par,1)/CHAR.(NUM.Adjoint.fields{par});
    end
end

% 4. denondimensionalize the variables
function [m] = denondimensionalize(NUM,m,npar,CHAR)
for par = 1:npar
    if strcmp(NUM.Adjoint.fields{par},'powerlaw') == 1
        m(par,1) = m(par,1);
    elseif strcmp(NUM.Adjoint.fields{par},'str_ref') == 1
        m(par,1) = m(par,1)*(1/CHAR.Time);
    elseif isfield(CHAR,NUM.Adjoint.fields{par}) == 1
        m(par,1) = m(par,1)*(CHAR.(NUM.Adjoint.fields{par}));
    end
end
