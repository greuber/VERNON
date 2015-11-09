function [ NUM ] = get_element_res( p,u,x,PAR ,MESH,NUM,i)
 
x = x +  u(NUM.Number.number_dof(:,1:9)) * PAR.dt;
NUM.f_quad_val   = zeros(18,1);
NUM.f_line_val   = zeros(4,1);
    
    for j =1:NUM.no_intp
        
        
        NUM.no_nodes_ele = 9;
        [ NI,dNI ] = Shape_Functions( NUM.no_nodes_ele,j,MESH.INTP.COORD );
        NUM.Solve.N_matrix = [NI(1,1), 0, NI(1,2), 0, NI(1,3), 0, NI(1,4), 0, NI(1,5), 0, NI(1,6), 0, NI(1,7), 0, NI(1,8), 0, NI(1,9),0;
            0, NI(1,1), 0, NI(1,2), 0, NI(1,3), 0, NI(1,4), 0, NI(1,5), 0, NI(1,6), 0, NI(1,7), 0, NI(1,8), 0, NI(1,9)];
        
        
        NUM.Solve.N   = NI;
        NUM.Solve.dNds = dNI;
        NUM.Solve.Jac    = NUM.Solve.dNds*x';
        
        % Inverse of Jacobian
        NUM.Solve.invJ = inv(NUM.Solve.Jac);
        NUM.Solve.detJ = det(NUM.Solve.Jac);
        NUM.Solve.dNdX = NUM.Solve.invJ*NUM.Solve.dNds;
        
        
        NUM.Solve.B(1,1:2:end) = NUM.Solve.dNdX(1,:);
        NUM.Solve.B(2,2:2:end)   = NUM.Solve.dNdX(2,:);
        NUM.Solve.B(3,1:2:end) = NUM.Solve.dNdX(2,:);
        NUM.Solve.B(3,2:2:end)   = NUM.Solve.dNdX(1,:);
        
        
        
        
    str_local   = NUM.Solve.B*u;
    exy = str_local(3,1)/2;
    str_invariant = sqrt(((1/2)*((str_local(1,1)^2+str_local(2,1)^2)))+exy^2);
    
    if str_invariant == 0
        str_invariant = mean(MESH.CompVar.str_ref(NUM.Number.number_quad(j,i)));
    end
    % Here comes the derivative of the function d/dn
   NUM.mu = ((MESH.CompVar.mu_ref(NUM.Number.number_quad(j,i)).*(str_invariant./MESH.CompVar.str_ref(NUM.Number.number_quad(j,i))).^(1./MESH.CompVar.powerlaw(NUM.Number.number_quad(j,i)) - 1)) * log(str_invariant./MESH.CompVar.str_ref(NUM.Number.number_quad(j,i))))/(MESH.CompVar.powerlaw(NUM.Number.number_quad(j,i))^2);
   NUM.mu_2d(NUM.Number.number_quad(j,i)) = NUM.mu;
        
        NUM.Solve.P = (1/2)*[ 2 ,        0               ,      0,0
            0                       , 2 ,      0,0
            0                       ,        0             , 1,1];
        
        NUM.Solve.D = 2*NUM.mu*NUM.Solve.P*NUM.Solve.P';
        
        NUM.no_nodes_ele = 4;
        [ NI,dNI ] = Shape_Functions( NUM.no_nodes_ele,j,MESH.INTP.COORD );
        NUM.Solve.N_p    = NI;
        NUM.no_nodes_ele = 9;
        
        NUM.Solve.Tau = NUM.Solve.D*str_local;
        NUM.Solve.P = NUM.Solve.N_p * p;
        NUM.Solve.Sigma = NUM.Solve.Tau - NUM.Solve.P*NUM.Solve.m;
        NUM.Solve.res_P = NUM.Solve.N_p'*(NUM.Solve.m'*str_local)* MESH.INTP.weight(j) * NUM.Solve.detJ;
        NUM.f_quad_val = NUM.f_quad_val + (NUM.Solve.B'*NUM.Solve.Sigma * MESH.INTP.weight(j) * NUM.Solve.detJ) - (NUM.Solve.N_matrix'*MESH.CompVar.rho(NUM.Number.number_quad(9,i))' * PAR.g * MESH.INTP.weight(j) * NUM.Solve.detJ) ;
        NUM.f_line_val = NUM.f_line_val - NUM.Solve.res_P;
        
        
        
    end
    
    
    


end

