function [H] = Hadamard(spin)
    dim = 2 * spin + 1;
    
    H = zeros(dim,dim);
    
    for n = 1:dim
        j = n-1;
        for m = 1:dim
            k = m-1;
            
            tmp = exp(1i * 2 * pi/dim * j * k);
            
            if (abs(real(tmp)) < 1e-12)
                tmp = 0 + 1i*imag(tmp);
            end
            if (abs(imag(tmp)) < 1e-12)
                tmp = real(tmp) + 1i*0;
            end
            
            H(n,m) = tmp;
        end
    end
    
    H = 1/sqrt(dim) * H;
end