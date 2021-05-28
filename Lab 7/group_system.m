function y = group_system(x,beta)

    N = length(x);
    y = zeros(N, 1);
    
    for n = 1:N
        if n-1 > 0
            y(n) = (1.5+beta).*y(n-1);
        end
        if n-2 > 0
            y(n) = y(n) + 1.5*beta.*y(n-2);
        end
        y(n) = y(n) + x(n);
    end

end


    

