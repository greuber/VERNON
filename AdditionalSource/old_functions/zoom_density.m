function [ beta_low,beta_step_star,beta_high,grad,fcost,dcost_du ] = zoom( beta_low,beta_high,NUM,MESH,PAR,dcost_du_old,dcost_du,fcost,fcost_old,beta_step,beta_step_old,sol_ini,dx,ind_X,fcost_ini,grad_ini)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
k = 1;
while k<10
%     % cubic interpolation
%     r1 = dcost_du_old + dcost_du - 3*((fcost_old - fcost)/(beta_step_old - beta_step));
%     r2 = sqrt(r1.^2 - dcost_du_old .* dcost_du);
%     beta_step_new = mean(beta_step - (beta_step - beta_step_old) * ((dcost_du + r2 -r1)./(dcost_du - dcost_du_old + 2*r2)));

beta_step = beta_step/10;   % give simple initial guess by decreasing beta  ( additionally get cubic to work)
    
    % Compute cost for beta_low
    MESH.CompVar.rho(unique(NUM.O_ind)) = sol_ini(unique(NUM.O_ind)) + beta_low*dx;
    [ NUM,MESH ] = Adjoint_StokesSolution( MESH,PAR,NUM );
    fcost_low = (1/2) * (NUM.r(ind_X) - NUM.r_ini(ind_X))' * (NUM.r(ind_X) - NUM.r_ini(ind_X));
    
    % Compute cost for the interpolated beta
    MESH.CompVar.rho(unique(NUM.O_ind)) = sol_ini(unique(NUM.O_ind)) + beta_step*dx;
    [ NUM,MESH ] = Adjoint_StokesSolution( MESH,PAR,NUM );
    fcost = (1/2) * (NUM.r(ind_X) - NUM.r_ini(ind_X))' * (NUM.r(ind_X) - NUM.r_ini(ind_X));
    
    if (fcost-1e-13) >= fcost_ini + 1e-4 * beta_step * grad_ini' * dx || fcost >= fcost_low       % minus -13 as equal is hard to get in matlab!
        beta_high = beta_step;
    else
        % Compute gradient
        dcost_du(ind_X) = (NUM.r(ind_X)' - NUM.r_ini(ind_X)');
        h_max = 1e-8;
        h = max([(1e-6*abs(MESH.CompVar.rho)), h_max]);
        [ drdx ] = Adjoint_residual( NUM,PAR,MESH,h );
        psi = (NUM.Solve.J')\(dcost_du');
        grad = - psi'*drdx;
        if abs(grad)-1e-10 <= abs(-0.9*grad_ini)     % same reason for the e-10 as the solution if it fits perfect never fits this condition
            beta_step_star = beta_step;
            break
        end
        if grad*(beta_high-beta_low)>0
            beta_high = beta_low;
        end
        beta_low = beta_step;
    end
%     beta_step_old = beta_step;
%     beta_step = beta_step_new;
    k = k+1;
    beta_step_star = beta_step;
end

if exist('grad','var')
else
    grad = grad_ini;
end

end

