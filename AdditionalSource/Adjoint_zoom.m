function [beta_step_star ] = Adjoint_zoom( beta_low,beta_high,NUM,MESH,PAR,sol_ini,dx,fcost_ini,grad_ini,dcost_du,fcost,CHAR)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
k = 1;
npar        = length(NUM.Adjoint.m);

while k<7
    %     % cubic interpolation
    %     r1 = dcost_du_old + dcost_du - 3*((fcost_old - fcost)/(beta_step_old - beta_step));
    %     r2 = sqrt(r1.^2 - dcost_du_old .* dcost_du);
    %     beta_step_new = mean(beta_step - (beta_step - beta_step_old) * ((dcost_du + r2 -r1)./(dcost_du - dcost_du_old + 2*r2)));
    
    beta_step = (1/2)*(beta_low + beta_high);   % give simple initial guess by decreasing beta  ( additionally get cubic to work)
    beta_step
    
    % Compute cost for beta_low
    NUM.Adjoint.m = sol_ini + beta_low.*dx;
    
    [NUM.Adjoint.m] = denormalize(NUM,NUM.Adjoint.m,npar);
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
    [NUM.Adjoint.m] = normalize(NUM,NUM.Adjoint.m,npar);
    
    [ NUM,MESH ] = Adjoint_StokesSolution( MESH,PAR,NUM,CHAR );
    
    fcost_low = (1/2) * (NUM.Solve.r(NUM.Adjoint.ind_cost) - NUM.Solve.r_ini(NUM.Adjoint.ind_cost))' * (NUM.Solve.r(NUM.Adjoint.ind_cost) - NUM.Solve.r_ini(NUM.Adjoint.ind_cost));
    
    % Compute cost for the interpolated beta
    NUM.Adjoint.m = sol_ini + beta_step.*dx;
    
    [NUM.Adjoint.m] = denormalize(NUM,NUM.Adjoint.m,npar);
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
    [NUM.Adjoint.m] = normalize(NUM,NUM.Adjoint.m,npar);
    
    [ NUM,MESH ] = Adjoint_StokesSolution( MESH,PAR,NUM,CHAR );
    fcost = (1/2) * (NUM.Solve.r(NUM.Adjoint.ind_cost) - NUM.Solve.r_ini(NUM.Adjoint.ind_cost))' * (NUM.Solve.r(NUM.Adjoint.ind_cost) - NUM.Solve.r_ini(NUM.Adjoint.ind_cost));
    
    
    if (fcost > fcost_ini + 1e-4 * beta_step .* grad_ini') | (fcost >= fcost_low)       % minus -13 as equal is hard to get in matlab!
        beta_high = beta_step;
    else
        % Compute gradient
        dcost_du(NUM.Adjoint.ind_cost) = (NUM.Solve.r(NUM.Adjoint.ind_cost)' - NUM.Solve.r_ini(NUM.Adjoint.ind_cost)');
        
        [NUM.Adjoint.m] = denormalize(NUM,NUM.Adjoint.m,npar);
        [NUM.Adjoint.m] = nondimensionalize(NUM,NUM.Adjoint.m,npar,CHAR);
        for par = 1:length(NUM.Adjoint.fields)   % loop over design variables
            % Compute residual by forward finite differences
            [ drdx_temp ] = AdjointRes( NUM,PAR,MESH,par,CHAR);
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
        
        if abs(grad) <= -0.9*grad_ini     % same reason for the e-10 as the solution if it fits perfect never fits this condition
            beta_step_star = beta_step;
            break
        end
        if grad*(norm(beta_high-beta_low))>=0
            beta_high = beta_low;
        end
        beta_low = beta_step;
        
    end
    k = k +1;
    beta_step_star = beta_step;
end

if exist('grad','var')
else
    grad = grad_ini;
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
    end
end