function [ NUM ] = Update_StressesAndStrains(NUM,PAR,MESH )

if NUM.Plasticity.Plasticity
    if isfield(NUM,'Tau_yield')
        % Correct stresses for plastic yielding
        ind_pl = find(NUM.Plastic  == 1);

        alpha = NUM.Tau_yield_2d(ind_pl) ./ NUM.Stress.T2nd_2d(ind_pl);
        NUM.Stress.Txx(ind_pl) = alpha.*NUM.Stress.Txx(ind_pl);
        NUM.Stress.Tzz(ind_pl) = alpha.*NUM.Stress.Tzz(ind_pl);
        NUM.Stress.Txz(ind_pl) = alpha.*NUM.Stress.Txz(ind_pl);
        NUM.Stress.T2nd = sqrt((1/2).*(NUM.Stress.Txx.^2+NUM.Stress.Tzz.^2)+NUM.Stress.Txz.^2);
    end
end

NUM.Stress.Txx_old = NUM.Stress.Txx;
NUM.Stress.Tzz_old = NUM.Stress.Tzz;
NUM.Stress.Txz_old = NUM.Stress.Txz;
NUM.Stress.T2nd_old = NUM.Stress.T2nd;

for i = 1:NUM.NUMERICS.no_elems_global     % Compute rotation angle from vorticity*dt
    LGCOORD = (MESH.GCOORD(:,NUM.Number.number_quad(:,i)))';
    
    for j =1:NUM.NUMERICS.no_intp
        
%         v_local = NUM.Solve.r(NUM.Number.number_ele_dof(1:NUM.NUMERICS.no_nodes_ele*2,i));
%         vx_local = v_local(1:2:end);
%         vz_local = v_local(2:2:end);
%         
%         dNds = NUM.Shape.dNI_vector{j};
%         Jac  = dNds*LGCOORD;
%         invJ = inv(Jac);
%         dNdX = invJ*dNds;
%         
%         RotationAngle = 0;
%         RotationAngle = RotationAngle -  0.5*(dNdX(2,j) .* vx_local(j,1)' - dNdX(1,j) .* vz_local(j,1)')  *PAR.dt;
%         Psi           = RotationAngle;
%         
%         
%         Txx_rot         =   NUM.Stress.Txx_old(j,i).*cos(Psi).^2 + NUM.Stress.Tzz_old(j,i).*sin(Psi).^2 - NUM.Stress.Txz_old(j,i).*sin(2*Psi);
%         Tzz_rot         =   NUM.Stress.Txx_old(j,i).*sin(Psi).^2 + NUM.Stress.Tzz_old(j,i).*cos(Psi).^2 + NUM.Stress.Txz_old(j,i).*sin(2*Psi);
%         Txz_rot         =   (NUM.Stress.Txx_old(j,i)-NUM.Stress.Tzz_old(j,i))/2.*sin(2*Psi)+ NUM.Stress.Txz_old(j,i).*(2*cos(Psi).^2 -1);
%         NUM.Stress.Txx_old(j,i)         =   Txx_rot;
%         NUM.Stress.Tzz_old(j,i)         =   Tzz_rot;
%         NUM.Stress.Txz_old(j,i)         =   Txz_rot;
%         
%         NUM.Stress.T2nd_old(j,i) = sqrt(NUM.Stress.Txx_old(j,i)^2+NUM.Stress.Tzz_old(j,i)^2+(2*NUM.Stress.Txz_old(j,i))^2);
%         
        if NUM.Plasticity.Plasticity
            % Compute plastic strain      
              NUM.Strain.E_pl_2nd(j,i) = (1-(NUM.Viscosity.mu_pl_2d(j,i)/NUM.Viscosity.mu_vis_2d(j,i))) * NUM.Strain.str_invariant_2d(j,i);
              if NUM.Strain.E_pl_2nd(j,i) < 0
                  NUM.Strain.E_pl_2nd(j,i) = 0;
              end
              
            % Compute Parameter for strain weakening
              NUM.Strain.APS(j,i)  = NUM.Strain.APS(j,i) + NUM.Strain.E_pl_2nd(j,i) * PAR.dt;
        end
    
    end
end


end