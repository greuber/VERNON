function [NUM,MESH,PAR,grad,grad_old,fcost,fcost_old,dfcost_du,dfcost_du_old,beta_step,beta_old_step,dn,dn_old,design_old] = Compute_Adjoint(NUM,MESH,PAR,grad,varargin)

if ~isempty(varargin)
    
    beta_step = varargin{1};
    beta_old_step = varargin{2};
    fcost_old= varargin{3};
    dfcost_du_old= varargin{4};
    dn_old= varargin{5};
    design_old= varargin{6};
    grad_old= varargin{7};
end

NUM.r = zeros(NUM.no_nodes*(NUM.ndof-1)+NUM.no_nodes_linear,1);
NUM.Solve.f_res_brutal = zeros(NUM.no_nodes*(NUM.ndof-1)+NUM.no_nodes_linear,1);
for iter = 1:100;   % say we make a maximum of 100 iterations
    
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
    
    % -----------------------------------------------
    % HERE COMES THE ADJOINT METHOD
    % ------------------------------------------------
    % IMPORTANT CODE IS SENSITIVE TO:
    % 1) the initial perturbation on the h_rho (and
    % h_max)
    % 2) initial beta_step
    % 3) factor to get the initial guess on design_new
    % 4) the factor of the residual (per defualt 1e-4
    % taken from literature)
    % ------------------------------------------------
    
    fcost = (1/2) * ((NUM.r - NUM.r_ini)' * (NUM.r - NUM.r_ini));
    % fcost = (NUM.r - NUM.r_ini);
    % fcost = fcost';  % as it has to be a 1xno_design_variables
    % dfcost_du = (-NUM.r_ini')*(-NUM.r_ini);
    dfcost_du = (NUM.r - NUM.r_ini)';   % as d/dx (x'x) = 2*x'
    
    
    display('------------------------------------------------------ ')
    display(['NORM OF THE GRADIENT = ',num2str(norm(grad))])
    display('------------------------------------------------------ ')
    
    h_max = 1e-8;
    h_rho = max([(1e-6*abs(MESH.CompVar.rho)) h_max]);
    %                         h_rho = max([(1e-6*norm(MESH.CompVar.rho)) h_max]);
    
    % Compute residual directly with perturbed design
    % variable
    
    [psize,qsize] = size(NUM.m);
    for id1 = 1:qsize  % loop for all design variables (as powerlaw is spatially constant one evalution of the whole domain is enough we diont have to distuinguish)
        for id2 = 1:2     % needs the perturbation plus and minus h per design variable
            
            if id2 == 1
                %                                     if id1 == 1 % compute powerlaw
                h = h_rho;
                MESH.CompVar.rho = MESH.CompVar.rho+h;
                %                                     else  % compute density
                %                                         h = h_rho;
                %                                         MESH.CompVar.rho = MESH.CompVar.rho+h;
                %                                     end
                
            else
                %                                    if id1 == 1 % compute powerlaw
                h = h_rho;
                MESH.CompVar.rho = MESH.CompVar.rho-2*h;
                %                                     else  % compute density
                %                                         h = h_rho;
                %                                         MESH.CompVar.rho = MESH.CompVar.rho-2*h;
                %                                     end
            end
            
            for i = 1:NUM.no_elems_global     % compute residual
                NUM.f_quad_val = zeros(18,1);
                NUM.f_line_val = zeros(4,1);
                u = NUM.r_ini(NUM.Number.number_ele_dof(1:18,i),1);
                p = NUM.r_ini(NUM.Number.number_ele_dof(19:22,i),1);
                x = MESH.Old(:,NUM.Number.number_quad(:,i));
                [ NUM ] = get_element_res(p,u,x,PAR ,MESH,NUM,i);
                NUM.Solve.f_res_brutal(NUM.Number.number_ele_dof(1:18,i),1) = NUM.Solve.f_res_brutal(NUM.Number.number_ele_dof(1:18,i),1) + NUM.f_quad_val;
                NUM.Solve.f_res_brutal(NUM.Number.number_ele_dof(19:22,i),1) = NUM.Solve.f_res_brutal(NUM.Number.number_ele_dof(19:22,i),1) + NUM.f_line_val;
            end
            % Add boundaries to the system
            for i = 1:1:length(NUM.Boundary.bcdof)
                NUM.Solve.f_res_brutal(NUM.Boundary.bcdof(i)) = 0;
            end
            
            if id2 == 1
                res_plus = NUM.Solve.f_res_brutal;
            else
                res_minus = NUM.Solve.f_res_brutal;
            end
        end
        drdx = (res_plus - res_minus)./(2*h);    % must be a NfxNd matrix
        
        % set the design variable back
        %                            if id1 == 1 % set back powerlaw
        MESH.CompVar.rho = MESH.CompVar.rho + h_rho;
        %                             else % set back rho
        %
        %                                 MESH.CompVar.rho = NUM.m(1:length(MESH.CompVar.powerlaw),2);
        %                             end
        
        
        dfcostdx = 0;   % as its always zero with respect to the design variables (depence implicitly on them)
        
        psi1 = (dfcost_du')\(NUM.Solve.J);   % fcost here is the derivative of the cost function with respect to u(solution)
        grad = dfcostdx - psi1'.*drdx;   % drdx = dresidual/ddesignvariables   %WHY THE HELL EVER BUT THIS IS STILL THE WRONG SIZE OF THE GRADIENT (SHOULD BE TRANSPOSED)
        
        
        %  COMPUTE APPROXIMATE INVERSE HESSIAN WITH A BFGS algorithm
        
        % 1.Hessian is just the identity matrix
        if NUM.it_adj == 1
            
            beta_step = norm(grad)/1e12;   % set the step_size initially to 1
            %                             else
            %                                 MESH.CompVar.rho = NUM.m(1:2:NUM.no_nodes*(NUM.ndof-1),2);
            %                                 MESH.CompVar.rho = MESH.CompVar.rho';
            %                             end
            
            
            H = eye(length(NUM.r),length(NUM.r));
            
            
            
        else
            %  2.Hessian is a first approximation
            % beta_step = norm(grad)/1e12;
            
            %                                     else
            %                                         MESH.CompVar.rho = NUM.m(1:2:NUM.no_nodes*(NUM.ndof-1),2);
            %                                         MESH.CompVar.rho = MESH.CompVar.rho';
            %                                     end
            %
            design_new = design_old(:,id1)*1e-12;
            
            vp = grad - grad_old;
            dp = design_new - design_old(:,id1);
            
            H = diag((vp'*dn_old)\(vp'*vp));
            
            % 3. real Hessian
            r = dp/(dp'*vp) - (H*vp)/(vp'*H*vp);
            H = H - ((H*vp) * (H*vp)')/(vp'*H*vp) + (dp*dp')/(dp'*vp) + (vp'*H*vp*(r*r'));
            
            
        end
        
        % now the perturbation is always computed
        % after solving the system and applied in the
        % next iteration to guarantee dp being nonzero
        
        dn = -H*grad;   % first compute the theoretical dn
        
        
        
        % end
        
        
        
        
        
        
        
        
    end
    
    
    %                         if NUM.it_adj ~= 1
    %                             if  fcost <= fcost_old + 1e-4 * beta_step * grad_old' * dn_old || norm(grad'*dn) >= 0.9 * norm(grad_old' * dn_old)
    %
    %                                 r1 = dfcost_du_old + dfcost_du - 3*((fcost_old - fcost)/(beta_old_step - beta_step));
    %                                 r2 = sqrt(r1.^2 - dfcost_du_old .* dfcost_du);
    %                                 beta_step = beta_step - (beta_step - beta_old_step) * ((dfcost_du + r2 -r1)/(dfcost_du - dfcost_du_old + 2*r2));
    %                             end
    %                         end
    % beta_step = beta_step + (beta_step-2*beta_step)*rand(1,1);   % to guarantee that wolfe is not negative
    
    
    
    NUM.m(:,id1) = NUM.m(:,id1) + beta_step * dn;
    % save old design variable vector
    design_old = NUM.m(:,id1);
    grad_old = grad;
    % if id1 == 1  % update powerlaw
    m1 = NUM.m(2:2:NUM.no_nodes*(NUM.ndof-1),1);
    % MESH.CompVar.rho = m1';
    MESH.CompVar.rho(unique(NUM.O_ind)) = m1(unique(NUM.O_ind));
    % MESH.CompVar.rho = MESH.CompVar.rho';
    
    % Update iteration index
    NUM.it_adj = NUM.it_adj + 1;
    fcost_old = fcost';
    
    beta_old_step = beta_step + (beta_step-beta_step/1000)*rand(1,1);
    
    %                         if NUM.it_adj == 1
    %                             grad_ini = grad;
    %                             NUM.adjoint_tol = grad_ini * 1e-3;   % to end when the error is 4 orders smaller than the in iitla one
    %                         end
    
    
    dn_old = dn;
    dfcost_du_old = dfcost_du;
    
    MESH.CompVar.rho
    
end




