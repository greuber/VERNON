function [ NUM] = Solver(NUM,PAR,MESH )

while NUM.Solve.number_pic>0
    
    switch NUM.Solve.Solver
        case 'Picard'
    [ NUM ] = get_globals_picard( NUM,PAR,MESH);
    
        case 'Jacobian'
            [ output_args ] = get_globals_Jacobian( input_args )
    
    NUM.Solve.f_res = NUM.Solve.L*NUM.r - NUM.Solve.FG;
    NUM.r = NUM.Solve.L\NUM.Solve.FG;
    
    if norm(NUM.Solve.f_res)<PAR.tol
        break
    end
    
    NUM.Solve.number_pic = NUM.Solve.number_pic+1;
    
end



end

