function Y = MODEL_main (x, func_in)

Q = func_in.Q;
HP = func_in.HP;

nt = round(x(4)); % define turbine types
On = round(x(5)); % define turbine numbers

switch On
    
    case 1 % 1 turbine
        [AR,  AC,  WEP_1] = calc_hpp1 ( x , HP , Q , nt, 1 );
        
    case 2 % 2 turbine
        [AR,  AC,  WEP_1] = calc_hpp2 ( x , HP , Q , nt, 2 );
        
    case 3 % 3 turbine
        [AR,  AC,  WEP_1] = calc_hpp3 ( x , HP , Q , nt, 3 );
        
end

    Y(1) = AC ;
    %
    Y(2) = -AR ;
    
    Y(3) = -WEP_1 ;
