function rho = deph(t,T2,rho0,gamma)
   
    rho = rho0;
        
    for m = 1:length(rho0)
        for mp = 1:length(rho0)
            if abs(m-mp)>1e-3
                rho(m,mp) = rho0(m,mp)*exp(t/T2*gamma(m,mp));
            end
        end
    end
    
end