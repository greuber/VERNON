function [ NI,dNI ] = Shape_Functions( no_nodes_ele,j,COORD )

% j is the index of the current intprop

%-----------------------------------
% 1D
%-----------------------------------
switch no_nodes_ele
    case  2
    
    N1 = (1/2) * (1+COORD(1,j));
    N2 = (1/2) * (1-COORD(1,j));
    
    dN1 = (1/2);
    dN2 = -(1/2);
    
    
    NI  = [N1,N2];
    dNI = [dN1,dN2];
    
    case  3
    
    N1 = (1/2)*COORD(1,j) * ((COORD(1,j)-1));
    N2 = 1-(COORD(1,j)^2)                   ;
    N3 = (1/2)*COORD(1,j) * ((COORD(1,j)+1));
    
    dN1 = COORD(1,j)-(1/2);
    dN2 = -2*COORD(1,j);
    dN3 = COORD(1,j) +(1/2);
    
    
    NI  = [N1,N2,N3];
    dNI = [dN1,dN2,dN3];
    
%-----------------------------------
% 2D
%----------------------------------- 
    case 1
        N1 = 1;
        
        dN1x = 0;
        dN1y = 0;
        
        
        NI  = [N1];
        dNI = [dN1x;
               dN1y];
        
        
    case 4
    
    N1 = (1/4) * (1-COORD(1,j)) * (1-COORD(2,j));
    N2 = (1/4) * (1-COORD(1,j)) * (1+COORD(2,j));
    N3 = (1/4) * (1+COORD(1,j)) * (1+COORD(2,j));
    N4 = (1/4) * (1+COORD(1,j)) * (1-COORD(2,j));
    
    dN1x = (-1/4) * (1-COORD(2,j));
    dN2x = (-1/4) * (1+COORD(2,j));
    dN3x = (1/4)  * (1+COORD(2,j));
    dN4x = (1/4)  * (1-COORD(2,j));
    dN1y = (-1/4) * (1-COORD(1,j));
    dN2y = (1/4)  * (1-COORD(1,j));
    dN3y = (1/4)  * (1+COORD(1,j));
    dN4y = (-1/4) * (1+COORD(1,j));
    
    
    NI  = [N1,N2,N3,N4];
    dNI = [dN1x,dN2x,dN3x,dN4x;
           dN1y,dN2y,dN3y,dN4y];
    
 
    
    case  9
    
    R1 = 0.5*(COORD(1,j))      * ((COORD(1,j)) - 1.0);
    R2 = -((COORD(1,j)) + 1.0) * ((COORD(1,j) - 1.0));
    R3 = 0.5*(COORD(1,j))      * ((COORD(1,j) + 1.0));
    S1 = 0.5*(COORD(2,j))      * ((COORD(2,j)) - 1.0);
    S2 = -((COORD(2,j)) + 1.0) * ((COORD(2,j)) - 1.0);
    S3 = 0.5*(COORD(2,j))      * ((COORD(2,j)) + 1.0);
    
    N1 = R1*S1;
    N2 = R1*S2;
    N3 = R1*S3;
    N4 = R2*S3;
    N5 = R3*S3;
    N6 = R3*S2;
    N7 = R3*S1;
    N8 = R2*S1;
    N9 = R2*S2;
    
    dR1 = (COORD(1,j)) - 0.5;
    dR2 = -2.0*(COORD(1,j));
    dR3 = (COORD(1,j)) + 0.5;
    dS1 = (COORD(2,j)) - 0.5;
    dS2 = -2.0*(COORD(2,j));
    dS3 = (COORD(2,j)) + 0.5;
    
    dN1x = dR1*S1;
    dN2x = dR1*S2;
    dN3x = dR1*S3;
    dN4x = dR2*S3;
    dN5x = dR3*S3;
    dN6x = dR3*S2;
    dN7x = dR3*S1;
    dN8x = dR2*S1;
    dN9x = dR2*S2;
    
    dN1y = R1*dS1;
    dN2y = R1*dS2;
    dN3y = R1*dS3;
    dN4y = R2*dS3;
    dN5y = R3*dS3;
    dN6y = R3*dS2;
    dN7y = R3*dS1;
    dN8y = R2*dS1;
    dN9y = R2*dS2;
    
    NI  = [N1,N2,N3,N4,N5,N6,N7,N8,N9];
    dNI = [dN1x,dN2x,dN3x,dN4x,dN5x,dN6x,dN7x,dN8x,dN9x;
           dN1y,dN2y,dN3y,dN4y,dN5y,dN6y,dN7y,dN8y,dN9y];

       %-----------------------------------
       % 3D
       %-----------------------------------
    case 27
    
    R1  = 0.5*COORD(1,j)*(COORD(1,j) - 1.0);
    R2  = -(COORD(1,j) + 1.0)*(COORD(1,j) - 1.0);
    R3  = 0.5*COORD(1,j)*(COORD(1,j) + 1.0);
    S1  = 0.5*COORD(2,j)*(COORD(2,j) - 1.0);
    S2  = -(COORD(2,j) + 1.0)*(COORD(2,j) - 1.0);
    S3  = 0.5*COORD(2,j)*(COORD(2,j) + 1.0);
    T1  = 0.5*COORD(3,j)*(COORD(3,j) - 1.0);
    T2  = -(COORD(3,j) + 1.0)*(COORD(3,j) - 1.0);
    T3  = 0.5*COORD(3,j)*(COORD(3,j) + 1.0);
    
    N1  = R1*S1*T1;
    N2  = R3*S1*T1;
    N3  = R3*S3*T1;
    N4  = R1*S3*T1;
    N5  = R1*S1*T3;
    N6  = R3*S1*T3;
    N7  = R3*S3*T3;
    N8  = R1*S3*T3;
    N9  = R2*S2*T2;
    N10 = R1*S1*T2;
    N11 = R1*S3*T2;
    N12 = R1*S2*T1;
    N13 = R1*S2*T3;
    N14 = R2*S1*T1;
    N15 = R2*S1*T3;
    N16 = R2*S3*T1;
    N17 = R2*S3*T3;
    N18 = R3*S1*T2;
    N19 = R3*S3*T2;
    N20 = R3*S2*T1;
    N21 = R3*S2*T3;
    N22 = R1*S2*T2;
    N23 = R3*S2*T2;
    N24 = R2*S1*T2;
    N25 = R2*S3*T2;
    N26 = R2*S2*T1;
    N27 = R2*S2*T3;
    
    dR1 = COORD(1,j) - 0.5;
    dR2 = -2.0*COORD(1,j);
    dR3 = COORD(1,j) + 0.5;
    dS1 = COORD(2,j) - 0.5;
    dS2 = -2.0*COORD(2,j);
    dS3 = COORD(2,j) + 0.5;
    dT1 = COORD(3,j) - 0.5;
    dT2 = -2.0*COORD(3,j);
    dT3 = COORD(3,j) + 0.5;
    
    dNI(1,:) = [ dR1*S1*T1, dR3*S1*T1, dR3*S3*T1, dR1*S3*T1, dR1*S1*T3, dR3*S1*T3, dR3*S3*T3, dR1*S3*T3, dR2*S2*T2, dR1*S1*T2, dR1*S3*T2, dR1*S2*T1, dR1*S2*T3, dR2*S1*T1, dR2*S1*T3, dR2*S3*T1, dR2*S3*T3, dR3*S1*T2, dR3*S3*T2, dR3*S2*T1, dR3*S2*T3, dR1*S2*T2, dR3*S2*T2, dR2*S1*T2, dR2*S3*T2, dR2*S2*T1, dR2*S2*T3 ];
    dNI(2,:) = [ R1*dS1*T1, R3*dS1*T1, R3*dS3*T1, R1*dS3*T1, R1*dS1*T3, R3*dS1*T3, R3*dS3*T3, R1*dS3*T3, R2*dS2*T2, R1*dS1*T2, R1*dS3*T2, R1*dS2*T1, R1*dS2*T3, R2*dS1*T1, R2*dS1*T3, R2*dS3*T1, R2*dS3*T3, R3*dS1*T2, R3*dS3*T2, R3*dS2*T1, R3*dS2*T3, R1*dS2*T2, R3*dS2*T2, R2*dS1*T2, R2*dS3*T2, R2*dS2*T1, R2*dS2*T3 ];
    dNI(3,:) = [ R1*S1*dT1, R3*S1*dT1, R3*S3*dT1, R1*S3*dT1, R1*S1*dT3, R3*S1*dT3, R3*S3*dT3, R1*S3*dT3, R2*S2*dT2, R1*S1*dT2, R1*S3*dT2, R1*S2*dT1, R1*S2*dT3, R2*S1*dT1, R2*S1*dT3, R2*S3*dT1, R2*S3*dT3, R3*S1*dT2, R3*S3*dT2, R3*S2*dT1, R3*S2*dT3, R1*S2*dT2, R3*S2*dT2, R2*S1*dT2, R2*S3*dT2, R2*S2*dT1, R2*S2*dT3 ];
    NI       = [N1,N2,N3,N4,N5,N6,N7,N8,N9,N10,N11,N12,N13,N14,N15,N16,N17,N18,N19,N20,N21,N22,N23,N24,N25,N26,N27];
    
    
    case  8
    
    N1 = (1/8) * (1-COORD(1,j)) * (1-COORD(2,j)) * (1-COORD(3,j));
    N2 = (1/8) * (1+COORD(1,j)) * (1-COORD(2,j)) * (1-COORD(3,j));
    N3 = (1/8) * (1+COORD(1,j)) * (1+COORD(2,j)) * (1-COORD(3,j));
    N4 = (1/8) * (1-COORD(1,j)) * (1+COORD(2,j)) * (1-COORD(3,j));
    N5 = (1/8) * (1-COORD(1,j)) * (1-COORD(2,j)) * (1+COORD(3,j));
    N6 = (1/8) * (1+COORD(1,j)) * (1-COORD(2,j)) * (1+COORD(3,j));
    N7 = (1/8) * (1+COORD(1,j)) * (1+COORD(2,j)) * (1+COORD(3,j));
    N8 = (1/8) * (1-COORD(1,j)) * (1+COORD(2,j)) * (1+COORD(3,j));
    
    
    dN1x = (-1/8) * (1-COORD(2,j)) * (1-COORD(3,j));
    dN2x = (1/8)  * (1-COORD(2,j)) * (1-COORD(3,j));
    dN3x = (-1/8) * (1+COORD(2,j)) * (1-COORD(3,j));
    dN4x = (1/8)  * (1+COORD(2,j)) * (1-COORD(3,j));
    dN5x = (1/8)  * (1-COORD(2,j)) * (1+COORD(3,j));
    dN6x = (-1/8) * (1-COORD(2,j)) * (1+COORD(3,j));
    dN7x = (1/8)  * (1+COORD(2,j)) * (1+COORD(3,j));
    dN8x = (-1/8) * (1+COORD(2,j)) * (1+COORD(3,j));
    
    dN1y = (-1/8) * (1-COORD(1,j)) * (1-COORD(3,j));
    dN2y = (1/8)  * (1+COORD(1,j)) * (1-COORD(3,j));
    dN3y = (-1/8) * (1+COORD(1,j)) * (1-COORD(3,j));
    dN4y = (1/8)  * (1-COORD(1,j)) * (1-COORD(3,j));
    dN5y = (1/8)  * (1-COORD(1,j)) * (1+COORD(3,j));
    dN6y = (-1/8) * (1+COORD(1,j)) * (1+COORD(3,j));
    dN7y = (1/8)  * (1+COORD(1,j)) * (1+COORD(3,j));
    dN8y = (-1/8) * (1-COORD(1,j)) * (1+COORD(3,j));
    
    dN1z = (-1/8) * (1-COORD(1,j)) * (1-COORD(2,j));
    dN2z = (1/8)  * (1+COORD(1,j)) * (1-COORD(2,j));
    dN3z = (-1/8) * (1+COORD(1,j)) * (1+COORD(2,j));
    dN4z = (1/8)  * (1-COORD(1,j)) * (1+COORD(2,j));
    dN5z = (1/8)  * (1-COORD(1,j)) * (1-COORD(2,j));
    dN6z = (-1/8) * (1+COORD(1,j)) * (1-COORD(2,j));
    dN7z = (1/8)  * (1+COORD(1,j)) * (1+COORD(2,j));
    dN8z = (-1/8) * (1-COORD(1,j)) * (1+COORD(2,j));
    
    
    NI  = [N1,N2,N3,N4,N5,N6,N7,N8];
    dNI = [dN1x,dN2x,dN3x,dN4x,dN5x,dN6x,dN7x,dN8x;
           dN1y,dN2y,dN3y,dN4y,dN5y,dN6y,dN7y,dN8y
           dN1z,dN2z,dN3z,dN4z,dN5z,dN6z,dN7z,dN8z];
    
end

