function [ weight,COORD ] = get_integration_points(no_intp)

if no_intp == 1;
    COORD(1,:) = [0];
    weight     = [2];

elseif no_intp == 2;
    COORD(1,:) = [-sqrt(1/3) , sqrt(1/3)];
    weight     = [1          ,     1    ];
    
elseif no_intp == 3;
    COORD(1,:) = [-(sqrt(3/5)) , 0 , sqrt(3/5)];
    weight     = [5/9          ,8/9,     5/9  ];
    
elseif no_intp == 4;
    COORD(1,:) = [-sqrt(1/3)  , -sqrt(1/3) , sqrt(1/3) ,  sqrt(1/3)];
    COORD(2,:) = [-sqrt(1/3)  ,  sqrt(1/3) , sqrt(1/3) , -sqrt(1/3)];
    weight     = [1           ,           1,          1,      1    ];
    
elseif no_intp == 9;
    COORD(1,:) = [-sqrt(3/5)  , -sqrt(3/5) , -sqrt(3/5) , 0        , sqrt(3/5) , sqrt(3/5), sqrt(3/5)  , 0          , 0   ];
    COORD(2,:) = [-sqrt(3/5)  , 0          , sqrt(3/5)  , sqrt(3/5), sqrt(3/5) , 0        , -sqrt(3/5) , -sqrt(3/5) , 0   ];
    weight     = [25/81       ,40/81       ,25/81       ,40/81     ,25/81      ,40/81     ,25/81       ,40/81       ,64/81];

elseif no_intp == 27;
    COORD(1,:) = [ -sqrt(3/5) , sqrt(3/5)  , sqrt(3/5)  , -sqrt(3/5) , -sqrt(3/5) , sqrt(3/5)  , sqrt(3/5) , -sqrt(3/5) ,           0 , -sqrt(3/5)  , -sqrt(3/5)  , -sqrt(3/5) , -sqrt(3/5) ,           0 ,           0 ,           0 ,           0 , sqrt(3/5)   , sqrt(3/5)   , sqrt(3/5)  , sqrt(3/5) , -sqrt(3/5)  , sqrt(3/5)   ,           0 ,           0 ,           0 ,           0 ];  
    COORD(2,:) = [ -sqrt(3/5) , -sqrt(3/5) ,  sqrt(3/5) ,  sqrt(3/5) , -sqrt(3/5) , -sqrt(3/5) , sqrt(3/5) , sqrt(3/5)  ,           0 , -sqrt(3/5)  ,  sqrt(3/5)  ,          0 ,         0  ,  -sqrt(3/5) , -sqrt(3/5)  ,  sqrt(3/5)  ,  sqrt(3/5)  , -sqrt(3/5)  ,   sqrt(3/5) ,          0 ,       0   ,           0 ,           0 , -sqrt(3/5)  ,  sqrt(3/5)  ,           0 ,           0 ]; 
    COORD(3,:) = [ -sqrt(3/5) , -sqrt(3/5) , -sqrt(3/5) , -sqrt(3/5) , sqrt(3/5)  , sqrt(3/5)  , sqrt(3/5) , sqrt(3/5)  ,           0 ,           0 ,           0 , -sqrt(3/5) , sqrt(3/5)  , -sqrt(3/5)  , sqrt(3/5)   , -sqrt(3/5)  , sqrt(3/5)   ,           0 ,           0 , -sqrt(3/5) , sqrt(3/5) ,           0 ,           0 ,           0 ,           0 , -sqrt(3/5)  , sqrt(3/5)   ];
    weight     = [125/729     , 125/729    , 125/729    , 125/729    , 125/729    , 125/729    , 125/729   , 125/729    , 512/729     , 200/729     , 200/729     , 200/729    , 200/729    , 200/729     , 200/729     , 200/729     , 200/729     , 200/729     , 200/729     , 200/729    , 200/729   , 320/729     , 320/729     , 320/729     , 320/729     , 320/729     ,   320/729 , ];
    
else
    error('Not implemented yet')
    
    
end

