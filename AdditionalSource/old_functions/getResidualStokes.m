function          [ f_quad_val,f_line_val ] = getResidualStokes( NUM,PAR,MESH,str_local,str_invariant,j,i);

        NUM.Solve.Tau = NUM.Solve.D*str_local;
        NUM.Solve.P = NUM.Solve.N_p*NUM.r(NUM.Number.number_ele_dof(19:22,i));
        NUM.Solve.Sigma = NUM.Solve.Tau - NUM.Solve.P*NUM.Solve.m;
        NUM.Solve.res_P = NUM.Solve.N_p'*(NUM.Solve.m'*str_local)* MESH.INTP.weight(j) * NUM.Solve.detJ;
        f_quad_val = f_quad_val + (NUM.Solve.B'*NUM.Solve.Sigma * MESH.INTP.weight(j) * NUM.Solve.detJ) - (NUM.Solve.N_matrix'*PAR.rho(number_quad(9,i))' * PAR.g * MESH.INTP.weight(j) * NUM.Solve.detJ) ;
        f_line_val = f_line_val -NUM.Solve.res_P;


end

