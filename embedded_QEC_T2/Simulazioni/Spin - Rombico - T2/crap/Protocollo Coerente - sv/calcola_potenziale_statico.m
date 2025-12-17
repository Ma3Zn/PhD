function V = calcola_potenziale_statico(sim, t)
%   Funzione che calcola il potenziale statico dovuto al campo trasverso da
%   misurare

%%% Per eseguire le simulazioni in interaction picture
%   Calcoliamo l'operatore da sch ad int
    U = diag(exp(1i * mod(diag(sim.sys.H0) * t, 2*pi)));

% % % %%% Per eseguire le simulazini in schrodinger picture
% % %     U = eye(sim.sys.dim);

%   Potenziale statico
    V = sim.Bx * (U * sim.sys.Mx * U');
%     V = sim.Bx * (U * sim.imp{4}.Mx * U');
end