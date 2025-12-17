function [v] = random_superposition(n)
    angle = randn(n,1);
    
    v = 1/sqrt(n) * exp(1i * angle);
end
