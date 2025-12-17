T = T_IDLE;

% ricaviamo il logaritmo dell'unitaria da implemntare
% Rz è diagonale nella base delle codewords e dunque è banale calcolarne il
% logaritmo in tale base
% log_Rz = diag(log(diag(Rz)));
log_Rz = logm(op);

% ricaviamo il termine noto dell'equazione
ratei = 1/T2 * T;

sop_E = 2*genera_superoperatore_sx(E_0) * genera_superoperatore_dx(E_0)...
        - genera_superoperatore_sx(E_0*E_0) ...
        - genera_superoperatore_dx(E_0*E_0);
sop_E = ratei * sop_E;

operatore = ( genera_superoperatore_sx(log_Rz) ...
            - genera_superoperatore_dx(log_Rz));
operatore = operatore + sop_E;

% ricaviamo l'operatore di evoluzione
% [V,v] = eig(operatore);
% exp_operatore = V*(diag(exp(diag(v))))/V;
exp_operatore = expm(operatore);

% simuliamo il gate
rho = exp_operatore * rho;