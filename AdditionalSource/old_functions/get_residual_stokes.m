function [ NUM ] = get_globals_picard( NUM,PAR,MESH)
    %% initilize globals
NUM.Solve.L = zeros(NUM.no_nodes*(NUM.ndof-1)+NUM.no_nodes_linear,NUM.no_nodes*(NUM.ndof-1)+NUM.no_nodes_linear);
NUM.Solve.FG = zeros(NUM.no_nodes*(NUM.ndof-1)+NUM.no_nodes_linear,1);
NUM.Solve.B  = zeros(3,NUM.no_nodes_ele*2);
NUM.Solve.F  = zeros(NUM.no_nodes_ele*(NUM.ndof-1)+NUM.no_nodes_ele_linear,1);
NUM.Solve.f_res = zeros(NUM.no_nodes*(NUM.ndof-1)+NUM.no_nodes_linear,1);

    
    
    for i = 1:NUM.no_elems_global
        
        % initialize local matrices
        NUM.Solve.KM = zeros(NUM.no_nodes_ele*(NUM.ndof-1),NUM.no_nodes_ele*(NUM.ndof-1));
        NUM.Solve.F  = zeros(NUM.no_nodes_ele*(NUM.ndof-1)+NUM.no_nodes_ele_linear,1);
        NUM.Solve.GM = zeros(NUM.no_nodes_ele*(NUM.ndof-1),NUM.no_nodes_ele_linear);
        NUM.Solve.LM = zeros(NUM.no_nodes_ele*(NUM.ndof-1)+NUM.no_nodes_ele_linear,NUM.no_nodes_ele*(NUM.ndof-1)+NUM.no_nodes_ele_linear);
        
        
        
        LGCOORD = (MESH.GCOORD(:,NUM.Number.number_quad(:,i)))';      % local global coordinates for this element
        
        for j =1:NUM.no_intp
            
            [ NI,dNI ] = Shape_Functions( NUM.no_nodes_ele,j,MESH.INTP.COORD );
            NUM.Solve.N_matrix = [NI(1,1), 0, NI(1,2), 0, NI(1,3), 0, NI(1,4), 0, NI(1,5), 0, NI(1,6), 0, NI(1,7), 0, NI(1,8), 0, NI(1,9),0;
                                    0, NI(1,1), 0, NI(1,2), 0, NI(1,3), 0, NI(1,4), 0, NI(1,5), 0, NI(1,6), 0, NI(1,7), 0, NI(1,8), 0, NI(1,9)];
            
            
            NUM.Solve.N   = NI;
            NUM.Solve.dNds = dNI;
            NUM.Solve.J    = NUM.Solve.dNds*LGCOORD;
            
            % Inverse of Jacobian
            NUM.Solve.invJ = inv(NUM.Solve.J);
            NUM.Solve.detJ = det(NUM.Solve.J);
            NUM.Solve.dNdX = NUM.Solve.invJ*NUM.Solve.dNds;
            
            
            NUM.Solve.B(1,1:2:end-1) = NUM.Solve.dNdX(1,:);
            NUM.Solve.B(2,2:2:end)   = NUM.Solve.dNdX(2,:);
            NUM.Solve.B(3,1:2:end-1) = NUM.Solve.dNdX(2,:);
            NUM.Solve.B(3,2:2:end)   = NUM.Solve.dNdX(1,:);
            
            
            
            
            [ mu,str_local] = ComputeViscosity(  NUM.Solve.B,PAR.mu_ref,PAR.str_ref,i,PAR.powerlaw,NUM.Number.number_ele_dof ,NUM.r );
            
            NUM.Solve.P = (1/2)*[ 2 ,        0               ,      0,0
                                0                       , 2 ,      0,0
                                0                       ,        0             , 1,1];
            
            NUM.Solve.D = 2*mu*NUM.Solve.P*NUM.Solve.P';
            
            NUM.Solve.KM = NUM.Solve.KM + NUM.Solve.B'*NUM.Solve.D*NUM.Solve.B * MESH.INTP.weight(j)*NUM.Solve.detJ;
            NUM.Solve.F(1:end-4,1) = NUM.Solve.F(1:end-4,1) + NUM.Solve.N_matrix'*PAR.rho(NUM.Number.number_quad(9,i))' * PAR.g * MESH.INTP.weight(j) * NUM.Solve.detJ;
            % Compute linear shap functions
            NUM.no_nodes_ele = 4;
            [ NI,dNI ] = Shape_Functions( NUM.no_nodes_ele,j,MESH.INTP.COORD );
            NUM.Solve.N_p    = NI;
            NUM.Solve.GM = NUM.Solve.GM - (NUM.Solve.B'*NUM.Solve.m*NUM.Solve.N_p* MESH.INTP.weight(j)*NUM.Solve.detJ);
            NUM.no_nodes_ele = 9;
            
            
        end
        
        NUM.Solve.LM(1:18,1:18) = NUM.Solve.KM;
        NUM.Solve.LM(1:18,19:22) = NUM.Solve.GM;
        NUM.Solve.LM(19:22,1:18) = NUM.Solve.GM';
       
        
        % add element stiffness matrix to global stiffness matrix
        NUM.Solve.L(NUM.Number.number_ele_dof(:,i),NUM.Number.number_ele_dof(:,i)) = NUM.Solve.L(NUM.Number.number_ele_dof(:,i),NUM.Number.number_ele_dof(:,i)) + NUM.Solve.LM;
        NUM.Solve.FG(NUM.Number.number_ele_dof(:,i),1) = NUM.Solve.FG(NUM.Number.number_ele_dof(:,i),1) + NUM.Solve.F;
        
    end
    
    % Add boundaries to the system
    for i = 1:1:length(NUM.Boundary.bcdof)
    
        NUM.Solve.L(NUM.Boundary.bcdof(i),:) = 0;
        NUM.Solve.L(NUM.Boundary.bcdof(i),NUM.Boundary.bcdof(i)) = 1;
        NUM.Solve.FG(NUM.Boundary.bcdof(i)) = NUM.Boundary.bcval(i);
    
    end


end

