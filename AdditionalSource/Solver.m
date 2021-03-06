function [ NUM,MESH] = Solver(NUM,PAR,MESH,CHAR )
%% --------------- %% Solver function %% --------------- %%
% Solves the Stokes equations with either analytical Jacobian or Picard
% matrix. Order:
% 1) Compute the chosen global matrix
% 2) Compute the residual element wise based on the current solution
% 3) Compute the update on the solution vector and finally update the
% solution
% 4) If the relative residual is lower then the tolerance update the
% stresses and strains
%%-------------------------------------------------------%%

NUM.Solve.number_pic = 1;
max_iter             = 2;
N_f_res_old          = realmax;
LS                   = 1;
no_picard            = 1;
NUM.Solve.Solver     = NUM.Solve.Solver_method;
res_ini              = 1;
number_new           = 1;

atol = 1e-15;
tol_picard = 5e-1;   % when to switch to Newton iterations

% apply bounds on the solution vector
NUM.Solve.r(NUM.Boundary.bcdof) = NUM.Boundary.bcval;

if isfield(NUM,'Convergence_normalized')
    NUM  = rmfield(NUM,'Convergence_normalized');
end
if isfield(NUM,'Convergence_total')
    NUM  = rmfield(NUM,'Convergence_total');
end

% get all shape function combinations at once
[NUM ] = Compute_ShapeFunctionsVector( MESH.INTP.COORD,NUM );

while NUM.Solve.number_pic<=max_iter
    
    % in case of the joint system the code has to decide the solver
    if strcmp(NUM.Solve.Solver_method,'NewtonPicard') == 1
        if NUM.Solve.number_pic <= no_picard
            NUM.Solve.Solver = 'Picard';
        elseif norm(NUM.Solve.f_res)/(res_ini) < tol_picard
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
            
% % %           Oldschool Picard
%             [ NUM,MESH ]    = get_globals_picard( NUM,PAR,MESH,CHAR);
%             NUM.Solve.f_res = NUM.Solve.L*NUM.Solve.r - NUM.Solve.FG;
%             NUM.Solve.r     = NUM.Solve.L\NUM.Solve.FG;

                    [ NUM,MESH ]     = get_globals_picard( NUM,PAR,MESH,CHAR);
                    [ NUM,MESH ]     = Compute_elemental_residual( NUM,MESH,PAR,CHAR );
                    dr               = NUM.Solve.L\(-NUM.Solve.f_res); 

                    alpha = 1;
                    r_0 = NUM.Solve.r;
                    
                    NUM.Solve.r      = r_0 + alpha * dr;
                    
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
                    [ NUM,MESH ]     = get_globals_Jacobian_analytical( NUM,PAR,MESH ,CHAR);
                    % [ NUM,MESH ]    = get_globals_picard( NUM,PAR,MESH,CHAR);
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
                            % [NUM.Adjoint.m] = normalize(NUM,NUM.Adjoint.m,npar);
                            
%                             % Compute cost function and derivative  %
%                             fcost = (1/2) * ((NUM.Solve.r(NUM.Adjoint.ind_cost)) - NUM.Solve.r_ini(NUM.Adjoint.ind_cost))' * ((NUM.Solve.r(NUM.Adjoint.ind_cost)) - NUM.Solve.r_ini(NUM.Adjoint.ind_cost));
%                             dcost_du = zeros(size(NUM.Solve.r'));
%                             dcost_du(NUM.Adjoint.ind_cost) = ((NUM.Solve.r(NUM.Adjoint.ind_cost))' - NUM.Solve.r_ini(NUM.Adjoint.ind_cost)');
                            
                            % Compute cost function and derivative  %
                            % ATTENTION THIS IS CHANGED
                            fcost = NUM.Solve.r(NUM.Adjoint.ind_cost);
                            dcost_du = zeros(size(NUM.Solve.r'));
                            dcost_du(NUM.Adjoint.ind_cost) = 1;
                            
                            % [NUM.Adjoint.m] = denormalize(NUM,NUM.Adjoint.m,npar);
                            [NUM.Adjoint.m] = nondimensionalize(NUM,NUM.Adjoint.m,npar,CHAR);
                            for par = 1:length(NUM.Adjoint.fields)   % loop over design variables
                                % Compute residual by forward finite differences
                                [ drdx_temp ] = Adjoint_Res( NUM,PAR,MESH,par,CHAR);
                                drdx(:,par) = drdx_temp;
                            end
                            [NUM.Adjoint.m] = denondimensionalize(NUM,NUM.Adjoint.m,npar,CHAR);
                            % [NUM.Adjoint.m] = normalize(NUM,NUM.Adjoint.m,npar);
                            
                            % compute gradient using adjoint method
                            psi = (NUM.Solve.J')\(dcost_du');
                            grad = - psi'*drdx;
                            % [grad] = nondimensionalize(NUM,grad',npar,CHAR);
                            % grad = grad';
                            % [grad] = normalize(NUM,grad',npar);
                            % grad = grad';
                            
                            % In the first iteration the Hessian is the
                            % identity matrix
                            NUM.Adjoint.H = eye(length(NUM.Adjoint.fields),length(NUM.Adjoint.fields));     % ATTENTION: THIS CAN BE NORMALIZED BY NUM.Adjoint.H*(1/norm(grad))
                            
                            % setup the tolerance
                            NUM.Adjoint.adjoint_tol = norm(grad)*1e-9;   % to end when the error is 9 orders smaller than the in iitla one
                            
                            % save step-size
                            beta_step = (1./abs(grad)')./NUM.Adjoint.LS_Parameter;
                            
                            beta_step_star = zeros(size(beta_step));
                            
                            % Compute first search-direction
                            dx = -NUM.Adjoint.H*grad';
                            
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
                            
                            %                             for par = 1:npar
                            %                                 if (sol_ini(par) + beta_step(par)*dx(par) >= NUM.Adjoint.bounds{par}(2))
                            %                                     beta_step(par) = beta_step(par)/dx(par);
                            %                                 end
                            %                             end
                            
                            
                            beta_step_max = 1;
                            beta_step_old = 0;
                            
                            % Starting the line-search
                            while 1>0
                                % Compute fcost with an initial guess scaled
                                % with the former beta_step
                                NUM.Adjoint.m = sol_ini + beta_step.*dx;
                                
                                % [NUM.Adjoint.m] = denormalize(NUM,NUM.Adjoint.m,npar);
                                [NUM.Adjoint.m] = nondimensionalize(NUM,NUM.Adjoint.m,npar,CHAR);
                                for par = 1:length(NUM.Adjoint.fields)
                                    if isfield(MESH.CompVar,NUM.Adjoint.fields{par}) == 1
                                        MESH.CompVar.(NUM.Adjoint.fields{par})(NUM.Adjoint.index{par}) = NUM.Adjoint.m(par,1);
                                    elseif strcmp(NUM.Adjoint.fields{par},'rad') == 1
                                        
                                        MESH.CompVar.rho         = ones(1,NUM.NUMERICS.no_nodes) * NUM.Adjoint.m(1,1);
                                        MESH.CompVar.mu_ref      = ones(1,NUM.NUMERICS.no_nodes) * NUM.Adjoint.m(4,1);
                                        MESH.CompVar.str_ref     = ones(1,NUM.NUMERICS.no_nodes) * (-NUM.ebg)/(1/CHAR.Time);
                                        MESH.CompVar.powerlaw    = ones(1,NUM.NUMERICS.no_nodes) * PAR.n1;
                                        MESH.CompVar.Phase       = ones(1,NUM.NUMERICS.no_nodes);
                                        MESH.CompVar.G           = ones(size(MESH.GCOORD(2,:)))*PAR.G/CHAR.Stress;
                                        rad = NUM.Adjoint.m(par,1);
                                        ind = find((MESH.GCOORD(1,:) - (PAR.W/2)).^2 + (MESH.GCOORD(2,:) - (PAR.H/2)).^2 < rad^2);
                                        MESH.CompVar.rho(ind)         = NUM.Adjoint.m(3,1);
                                        MESH.CompVar.mu_ref(ind)      = PAR.mu_2/CHAR.Viscosity ;
                                        MESH.CompVar.str_ref(ind)     = (-NUM.ebg)/(1/CHAR.Time);
                                        MESH.CompVar.powerlaw(ind)    = PAR.n1;
                                        MESH.CompVar.Phase(ind)       = 2;
                                        MESH.CompVar.G(ind)           = PAR.G/CHAR.Stress;
                                        
                                    elseif isfield(PAR,NUM.Adjoint.fields{par}) == 1    
                                        PAR.(NUM.Adjoint.fields{par})(NUM.Adjoint.index{par}) = NUM.Adjoint.m(par,1);
                                    else
                                        display('Design variable seems to not be defined for the variable update')
                                    end
                                end
                                [NUM.Adjoint.m] = denondimensionalize(NUM,NUM.Adjoint.m,npar,CHAR);
                                % [NUM.Adjoint.m] = normalize(NUM,NUM.Adjoint.m,npar);
                                
                                [ NUM,MESH ] = Adjoint_StokesSolution( MESH,PAR,NUM,CHAR );
                                fcost = (1/2) * (NUM.Solve.r(NUM.Adjoint.ind_cost) - NUM.Solve.r_ini(NUM.Adjoint.ind_cost))' * (NUM.Solve.r(NUM.Adjoint.ind_cost) - NUM.Solve.r_ini(NUM.Adjoint.ind_cost));
                                
                                % First condition ('sufficient decrease
                                % condition')
                                if (fcost > fcost_ini + 1e-4 * beta_step .* grad_ini')  | (fcost >= fcost_old)
                                    [beta_step_star] = Adjoint_zoom( beta_step_old,beta_step,NUM,MESH,PAR,sol_ini,dx,fcost_ini,grad_ini,dcost_du,fcost,CHAR);
                                    beta_step_star = beta_step_star;
                                    break
                                end
                                
                                % Compute derivative of cost and gradient
                                % with the current step-size
                                dcost_du = zeros(size(NUM.Solve.r'));
                                dcost_du(NUM.Adjoint.ind_cost) = (NUM.Solve.r(NUM.Adjoint.ind_cost)' - NUM.Solve.r_ini(NUM.Adjoint.ind_cost)');
                                
                                % [NUM.Adjoint.m] = denormalize(NUM,NUM.Adjoint.m,npar);
                                [NUM.Adjoint.m] = nondimensionalize(NUM,NUM.Adjoint.m,npar,CHAR);
                                for par = 1:length(NUM.Adjoint.fields)   % loop over design variables
                                    % Compute residual by forward finite differences
                                    [ drdx_temp ] = Adjoint_Res( NUM,PAR,MESH,par,CHAR);
                                    drdx(:,par) = drdx_temp;
                                end
                                [NUM.Adjoint.m] = denondimensionalize(NUM,NUM.Adjoint.m,npar,CHAR);
                                % [NUM.Adjoint.m] = normalize(NUM,NUM.Adjoint.m,npar);
                                
                                % Compute gradient using adjoint method
                                psi = (NUM.Solve.J')\(dcost_du');
                                grad = - psi'*drdx;
                                % [grad] = nondimensionalize(NUM,grad',npar,CHAR);
                                % grad = grad';
                                % [grad] = normalize(NUM,grad',npar);
                                % grad = grad';
                                
                                % Second condition ('curvature condition')
                                if abs(grad) <= -0.9*grad_ini           % MAYBE COULD BE any(abs(grad) <= abs(-0.9*grad_ini))
                                    beta_step_star = beta_step;
                                    break
                                end
                                %                                 if grad >= 0        % CHANGED       % BUT GRADIENT IS STILL starting OSCILLATING          % MAYBE COULD BE any(abs(grad) <= abs(-0.9*grad_ini))
                                %                                     [beta_step_star] = Adjoint_zoom( beta_step,beta_step_old,NUM,MESH,PAR,sol_ini,dx,fcost_ini,grad_ini,dcost_du,fcost,CHAR);
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
                                %                                 ind = find(beta_step >= beta_step_max-1e-3);
                                %                                 if ~isempty(ind)
                                %                                     beta_step_star = beta_step;
                                %                                     beta_step_star(ind) = beta_step_max;
                                %                                     break
                                %                                 end
                                
                                %                                 beta_step_new = real(beta_step_new);
                                % beta_step_old = beta_step;
                                beta_step = beta_step + beta_step*1.1; % + (beta_step_max-beta_step)*rand(1);
                                fcost_old = fcost;
                                ij = ij +1
                                
                                % beta_step
                                
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
                                
                                % [NUM.Adjoint.m] = denormalize(NUM,NUM.Adjoint.m,npar);
                                [NUM.Adjoint.m] = nondimensionalize(NUM,NUM.Adjoint.m,npar,CHAR);
                                for par = 1:length(NUM.Adjoint.fields)   % loop over design variables
                                    % Compute residual by forward finite differences
                                    [ drdx_temp ] = Adjoint_Res( NUM,PAR,MESH,par,CHAR);
                                    drdx(:,par) = drdx_temp;
                                end
                                [NUM.Adjoint.m] = denondimensionalize(NUM,NUM.Adjoint.m,npar,CHAR);
                                % [NUM.Adjoint.m] = normalize(NUM,NUM.Adjoint.m,npar);
                                
                                % Compute gradient using adjoint method
                                psi = (NUM.Solve.J')\(dcost_du');
                                grad = - psi'*drdx;
                                % [grad] = nondimensionalize(NUM,grad',npar,CHAR);
                                % grad = grad';
                                % [grad] = normalize(NUM,grad',npar);
                                % grad = grad';
                            end
                            
                            
                            NUM.Adjoint.m = sol_ini + beta_step.*dx;
                            % [NUM.Adjoint.m] = denormalize(NUM,NUM.Adjoint.m,npar);
                            [NUM.Adjoint.m] = nondimensionalize(NUM,NUM.Adjoint.m,npar,CHAR);
                            for par = 1:length(NUM.Adjoint.fields)
                                if isfield(MESH.CompVar,NUM.Adjoint.fields{par}) == 1
                                    MESH.CompVar.(NUM.Adjoint.fields{par})(NUM.Adjoint.index{par}) = NUM.Adjoint.m(par,1);
                                elseif strcmp(NUM.Adjoint.fields{par},'rad') == 1
                                    
                                    MESH.CompVar.rho         = ones(1,NUM.NUMERICS.no_nodes) * NUM.Adjoint.m(1,1);
                                    MESH.CompVar.mu_ref      = ones(1,NUM.NUMERICS.no_nodes) * NUM.Adjoint.m(4,1) ;
                                    MESH.CompVar.str_ref     = ones(1,NUM.NUMERICS.no_nodes) * (-NUM.ebg)/(1/CHAR.Time);
                                    MESH.CompVar.powerlaw    = ones(1,NUM.NUMERICS.no_nodes) * PAR.n1;
                                    MESH.CompVar.Phase       = ones(1,NUM.NUMERICS.no_nodes);
                                    MESH.CompVar.G           = ones(size(MESH.GCOORD(2,:)))*PAR.G/CHAR.Stress;
                                    rad = NUM.Adjoint.m(par,1);
                                    ind = find((MESH.GCOORD(1,:) - (PAR.W/2)).^2 + (MESH.GCOORD(2,:) - (PAR.H/2)).^2 < rad^2);
                                    MESH.CompVar.rho(ind)         = NUM.Adjoint.m(3,1);
                                    MESH.CompVar.mu_ref(ind)      = PAR.mu_2/CHAR.Viscosity ;
                                    MESH.CompVar.str_ref(ind)     = (-NUM.ebg)/(1/CHAR.Time);
                                    MESH.CompVar.powerlaw(ind)    = PAR.n1;
                                    MESH.CompVar.Phase(ind)       = 2;
                                    MESH.CompVar.G(ind)           = PAR.G/CHAR.Stress;
                                    
%                                 elseif strcmp(NUM.Adjoint.fields{par},'Hi') == 1
%                                     
%                                     PAR.H_interface = NUM.Adjoint.m(par,1);
%                                     
%                                     [ NUM,MESH ] = Create_MESH( NUM,PAR,0,MESH );
%                                     
%                                     MESH.CompVar.rho         = ones(size(MESH.GCOORD(2,:)))*0;
%                                     MESH.CompVar.mu_ref      = ones(size(MESH.GCOORD(2,:)))*PAR.mu1;
%                                     MESH.CompVar.str_ref     = ones(size(MESH.GCOORD(2,:)))*abs(NUM.ebg);
%                                     MESH.CompVar.powerlaw    = ones(size(MESH.GCOORD(2,:)))*PAR.n;
%                                     MESH.CompVar.Phase       = ones(size(MESH.GCOORD(2,:)))*1;
%                                     MESH.CompVar.G           = ones(size(MESH.GCOORD(2,:)))*PAR.G;
%                                     
%                                     hi = NUM.Adjoint.m(par,1);
%                                     % Overburden
%                                     ind                 = find(MESH.GCOORD(2,:)<=1 & MESH.GCOORD(2,:)>hi+PAR.A0 & MESH.GCOORD(1,:)>=0 & MESH.GCOORD(1,:)<=PAR.W);
%                                     MESH.CompVar.rho(ind)        = PAR.rho;
%                                     MESH.CompVar.mu_ref(ind)     = PAR.mu2;
%                                     MESH.CompVar.str_ref(ind)    = abs(NUM.ebg);
%                                     MESH.CompVar.powerlaw(ind)   = PAR.n;
%                                     MESH.CompVar.Phase(ind)      = 2;
%                                     MESH.CompVar.G(ind)          = PAR.G;
                                    
                                elseif isfield(PAR,NUM.Adjoint.fields{par}) == 1
                                    PAR.(NUM.Adjoint.fields{par})(NUM.Adjoint.index{par}) = NUM.Adjoint.m(par,1);
                                else
                                    display('Design variable seems to not be defined for the variable update')
                                end
                            end
                            [NUM.Adjoint.m] = denondimensionalize(NUM,NUM.Adjoint.m,npar,CHAR);
                            % [NUM.Adjoint.m] = normalize(NUM,NUM.Adjoint.m,npar);
                            
                            % Compute parameters for BFGS
                            dp = NUM.Adjoint.m - sol_ini;
                            vp = grad - grad_ini;
                            
                            NUM.Adjoint.H = eye(length(NUM.Adjoint.fields),length(NUM.Adjoint.fields));
                            
                            r = (dp./(dp.*vp')) - ((NUM.Adjoint.H*vp')/(vp*NUM.Adjoint.H*vp'));
                            
                            % Approximate inverse Hessian with BFGS
                            % algorithm
                            NUM.Adjoint.H = NUM.Adjoint.H - ((NUM.Adjoint.H*vp'*(NUM.Adjoint.H*vp')')/(vp*NUM.Adjoint.H*vp')) + ((dp*dp')/(dp'*vp')) + (vp*NUM.Adjoint.H*vp'*(r*r'));   % = NUM.Adjoint.H^-1
                            
                            NUM.Adjoint.H(isnan(NUM.Adjoint.H)) = 0;
                            NUM.Adjoint.H(isinf(NUM.Adjoint.H)) = 0;
                            
                            inv(NUM.Adjoint.H)   % print the correct Hessian (as BFGS gives the inverse Hessian)
                            
                            
                            % check for positive definiteness of NUM.Adjoint.H
                            [~,p] = chol(NUM.Adjoint.H);
                            if p>0
                                warning('Hessian no longer positive definite')
                            end
                            
                            % Compute new search-direction
                            dx = -NUM.Adjoint.H*grad';
                            
                            % dx    % print the search direction vector
                            
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
                        
%                         figure(1),clf
%                         for par = 1:length(NUM.Adjoint.fields)
%                             if isfield(MESH.CompVar,(NUM.Adjoint.fields{par})) == 1
%                                 design1_2d = MESH.CompVar.(NUM.Adjoint.fields{par})(NUM.Number.number_2d);
%                                 subplot(length(NUM.Adjoint.fields),1,par)
%                                 pcolor(X,Z,design1_2d)
%                             elseif isfield(PAR,(NUM.Adjoint.fields{par})) == 1
%                                 design1 = PAR.(NUM.Adjoint.fields{par});
%                                 subplot(length(NUM.Adjoint.fields),1,par)
%                                 plot(design1,'ro','MarkerSize',5)
%                             end
%                             colorbar
%                             title(['design variable ',num2str(par),' after ',num2str(NUM.Adjoint.it_adj),' iteration; Norm cost function = ',num2str(norm(fcost))])
%                         end
%                         drawnow
                        
                        % fname = (['Iteration_jioint_3_',num2str(NUM.Adjoint.it_adj)]);
                        % print(fname,'-dpng','-zbuffer','-r300')
                        
                        % MESH.CompVar.powerlaw
                        % MESH.CompVar.rho*CHAR.rho
                        
                        % rho_sur_temp = mean(MESH.CompVar.rho(unique((NUM.Adjoint.index{3}))));
                        % power_sur_temp = mean(MESH.CompVar.powerlaw(unique((NUM.Adjoint.index{1}))));
                        % power_block_temp = mean(MESH.CompVar.powerlaw(unique((NUM.Adjoint.index{4}))));
                        % mu_temp = mean(MESH.CompVar.mu_ref(unique((NUM.Adjoint.index{2}))));
                        % save(['Adjoint_vectorized_4_',num2str(NUM.Adjoint.it_adj),'Rising_Sphere'],'power_sur_temp','mu_temp','rho_sur_temp','power_block_temp')
                        
                        
                        % [NUM.Adjoint.m] = denormalize(NUM,NUM.Adjoint.m,npar);
                        % [grad] = denormalize(NUM,grad',npar);
                        
                        display('------------------------------------------------------ ')
                        display(['NORM OF THE GRADIENT = ',num2str(norm(grad))])
                        display(['COST FUNCTION        = ',num2str(fcost)])
                        display(['BETA STEP            = ',num2str(beta_step')])
                        display(['mu_2 = ',num2str(NUM.Adjoint.m(1))     ])
                        display('------------------------------------------------------ ')
                        
                        if NUM.Adjoint.OnlyGradients == 1
                            NUM.Adjoint.grad       = grad;
                            NUM.Adjoint.norm_grad  = norm(grad);
                            NUM.Adjoint.fcost      = fcost;
                            break
                        end
                        
                        % [NUM.Adjoint.m] = normalize(NUM,NUM.Adjoint.m,npar);
                        
                        % update old variables
                        NUM.Adjoint.it_adj = NUM.Adjoint.it_adj + 1;
                        fcost_old = fcost;
                        dcost_du_old = dcost_du;
                        
                        % % denormalize the gradient for the stopping
                        % % criteria
                        % [grad] = denormalize(NUM,grad',npar);
                        % grad   = grad';
                        
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
                            f_quad_val = zeros(NUM.NUMERICS.no_nodes_ele*2,1);
                            f_line_val = zeros(NUM.NUMERICS.no_nodes_ele_linear,1);
                            tmp = u(j);
                            u(j) = u(j) + PAR.perturb ;
                            [ NUM,f_quad_val,f_line_val,MESH ] = get_element_res( p,u,PAR ,MESH,NUM,i,CHAR,f_quad_val,f_line_val);
                            f_res_new(1:NUM.NUMERICS.no_nodes_ele*2,1)  = f_quad_val;
                            f_res_new(NUM.NUMERICS.no_nodes_ele*2+1:NUM.NUMERICS.no_nodes_ele*2+NUM.NUMERICS.no_nodes_ele_linear,1) = f_line_val;
                            J(:,j) = (f_res_new - f_res_old) /PAR.perturb;
                            u(j) = tmp ;
                        end
                        
                        for j = 1:NUM.NUMERICS.no_nodes_ele_linear
                            f_quad_val = zeros(NUM.NUMERICS.no_nodes_ele*2,1);
                            f_line_val = zeros(NUM.NUMERICS.no_nodes_ele_linear,1);
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
                        J_global(:,NUM.Boundary.bcdof(i)) = 0;
                        J_global(NUM.Boundary.bcdof(i),NUM.Boundary.bcdof(i)) = 1;
                        NUM.Solve.f_res(NUM.Boundary.bcdof(i)) = 0;
                    end
                    
                    dr    = J_global\(-NUM.Solve.f_res);
                    
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
                            temp            = J_global*dr;
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
                    number_new = number_new + 1;    
                    
            end
    end
    
    NUM.Convergence_total(1,NUM.Solve.number_pic) = norm(NUM.Solve.f_res);
    NUM.Convergence_normalized(1,NUM.Solve.number_pic) = norm(NUM.Solve.f_res)/(res_ini);
    
    if NUM.Solve.number_pic == 1
        res_ini = norm(NUM.Solve.f_res);
    end
    
    if norm(NUM.Solve.f_res)<1e-14
        display(sprintf('Convergence due to machine precision; norm = %6.3e',norm(NUM.Solve.f_res)))
        break
    end
    if norm(NUM.Solve.f_res)/(res_ini) < PAR.tol + atol || NUM.Solve.number_pic >=max_iter   % Stop if we reach the stopping criteria || maximum iteration count
        % Stresses and Strains are updated in 'Timestepping'
        break
    end
    
    if norm(NUM.Solve.f_res)/(res_ini) > 1e10
        display('Code did not converge')
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
    elseif strcmp(NUM.Adjoint.fields{par},'mu_ref') == 1 || strcmp(NUM.Adjoint.fields{par},'mu') == 1
        m(par,1) = m(par,1)/CHAR.Viscosity;
    elseif isfield(CHAR,NUM.Adjoint.fields{par}) == 1
        m(par,1) = m(par,1)/CHAR.(NUM.Adjoint.fields{par});
        
    elseif strcmp(NUM.Adjoint.fields{par},'g') == 1
        m(par,1) = m(par,1)/CHAR.Gravity;
    elseif strcmp(NUM.Adjoint.fields{par},'rad') == 1
        m(par,1) = m(par,1)/CHAR.Length;
    elseif strcmp(NUM.Adjoint.fields{par},'Hi') == 1
        m(par,1) = m(par,1)/CHAR.Length;
    end
end

% 4. denondimensionalize the variables
function [m] = denondimensionalize(NUM,m,npar,CHAR)
for par = 1:npar
    if strcmp(NUM.Adjoint.fields{par},'powerlaw') == 1
        m(par,1) = m(par,1);
    elseif strcmp(NUM.Adjoint.fields{par},'str_ref') == 1
        m(par,1) = m(par,1)*(1/CHAR.Time);
    elseif strcmp(NUM.Adjoint.fields{par},'mu_ref') == 1 || strcmp(NUM.Adjoint.fields{par},'mu') == 1
        m(par,1) = m(par,1)*CHAR.Viscosity;
    elseif isfield(CHAR,NUM.Adjoint.fields{par}) == 1
        m(par,1) = m(par,1)*(CHAR.(NUM.Adjoint.fields{par}));
        
    elseif strcmp(NUM.Adjoint.fields{par},'g') == 1
        m(par,1) = m(par,1)*CHAR.Gravity;
    elseif strcmp(NUM.Adjoint.fields{par},'rad') == 1
        m(par,1) = m(par,1)*CHAR.Length;
    elseif strcmp(NUM.Adjoint.fields{par},'Hi') == 1
        m(par,1) = m(par,1)*CHAR.Length;
    end
end
