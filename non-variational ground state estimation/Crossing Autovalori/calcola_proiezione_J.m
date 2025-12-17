function J_effettivi = calcola_proiezione_J(n, J, k)
%   Funzoine per il calcolo della proiezione dei J da simulare tra spin 1 e
%   spin 1/2 sul sistema composto da soli spin 1/2

%   Coefficiente di proiezione da spin 1 a coppia di spin 1/2
    cp = 2;

%   Allochiamo lo spazio per i J_effettivi
    J_effettivi = zeros(1, 3/2*n);

%   Verifichiamo che n sia pari
    if (mod(n,2)==1)
        return;
    end

%   Indice di supporto
    l = 1;

%   Cicliamo su tutti i J
    for i = 1:3:3/2 * n
%       Siccome la catena è simmetrica supponiamo di partire con uno spin
%       1/2
        J_effettivi(i) = cp * J(l);
        
%       Il prossimo J_effettivo è quello della coppia che simula lo spin 1
        J_effettivi(i+1) = -k * (J(l) + J(l+1))/2;

%       Altro spin 1/2 della catena
        J_effettivi(i+2) = cp * J(l+1);

%       Aggiorniamo l'indice l
        l = l + 2;
    end
end