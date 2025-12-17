function [sys, seq_imp] = genera_sequenza_impulsi(sys)
%   Funzione che genera la sequenza di impulsi lungo z da mandare al
%   sistema per renderne degeneri gli opportuni livelli energetici in modo
%   da riuscire ad eseguire una rotazione logica sotto l'effetto del campo
%   statico incognito Bx

%   Recuperiamo la dimensione del sistema
    dim = sys.dim;

%   Calcoliamo il gap esatto (per poter usare dei Bx più grandi)
    H = sys.H0 + sys.Bx * sys.Mx;
    [V,d] = eig(H);
    d = sort(diag(d));
    d = diag(d(end:-1:1));
    sys.H = d;

% % % %   Non capisco se sia necessario o meno
% % %     sys.Mz = V' * sys.Mz * V;

%%  Gate logico da eseguire

%   Recuperiamo il generatore di tale gate
    h = 1i*logm(sys.GL);

%   Riscriviamo il generatore nella base computazionale
    h = sys.CB * h * sys.CB';

%%  Calcoliamo il commutatore di Mx ed Mz
    comm_zx = (sys.mu_b / sys.h_bar) * diag(diag(sys.Sz))*sys.Sx - sys.Sx*diag(diag(sys.Sz));

%%  Parametri impulsi da eseguire
%   Allochiamo una sequenza di impulsi da eseguire in parallelo
    imp_seq = cell(1,4);
    
%   Settiamo il valore massimo per il campo magnetico utilizzabile per
%   l'impulso
    B_max = sys.Bz_max; %T

%%  Primo Impulso   
    lv = [1 2];

    imp_1.lv      = lv;
    imp_1.w       = sys.w(lv(1),lv(2));
    imp_1.angolo  = abs(h(lv(1),lv(2)));
    imp_1.fase    = angle(h(lv(1),lv(2))) - angle(comm_zx(lv(1),lv(2)));
% % %     imp_1.angolo  = pi/2; %%%%%%%%%%%%%%% DBG
% % %     imp_1.fase    = pi/2; %%%%%%%%%%%%%%% DBG

% % % %   Rimozione completa leackage
% % %     imp_1.Mz               = zeros(dim, dim);
% % %     imp_1.Mz(lv(1), lv(1)) = sys.Mz(lv(1), lv(1));
% % %     imp_1.Mz(lv(2), lv(2)) = sys.Mz(lv(2), lv(2));

%   Rimozione parziale leackage
    imp_1.Mz    = zeros(dim, dim);
    imp_1.Mz    = diag(diag(sys.Mz));

% % %     imp_1.Mx               = zeros(dim, dim);
% % %     imp_1.Mx(lv(1), lv(2)) = sys.Mx(lv(1),lv(2));
% % %     imp_1.Mx(lv(2), lv(1)) = sys.Mx(lv(2),lv(1));

% % % %   Calcoliamo il gap esatto (per poter usare dei Bx più grandi)
% % %     H = sys.H0 + sys.Bx * sys.Mx;
% % %     [~,d] = eig(H);
% % %     d = sort(diag(d));
% % %     d = diag(d(end:-1:1));
% % %     sys.H = d;
    imp_1.w = d(lv(1),lv(1)) - d(lv(2),lv(2));
    
% % %     imp_1.w = sys.w(lv(1),lv(2));

%%  Secondo Impulso   
    lv = [1 3];

    imp_2.lv      = lv;
    imp_2.w       = sys.w(lv(1),lv(2));
    imp_2.angolo  = abs(h(lv(1),lv(2)));
    imp_2.fase    = angle(h(lv(1),lv(2))) - angle(comm_zx(lv(1),lv(2)));
% % %     imp_2.angolo  = pi/2; %%%%%%%%%%%%%%% DBG
% % %     imp_2.fase    = pi/2; %%%%%%%%%%%%%%% DBG
    
% % % %   Rimozione completa leackage
% % %     imp_2.Mz               = zeros(dim, dim);
% % %     imp_2.Mz(lv(1), lv(1)) = sys.Mz(lv(1), lv(1));
% % %     imp_2.Mz(lv(2), lv(2)) = sys.Mz(lv(2), lv(2));
 
%   Rimozione parziale leackage
    imp_2.Mz    = zeros(dim, dim);
    imp_2.Mz    = diag(diag(sys.Mz));

% % %     imp_2.Mx               = zeros(dim, dim);
% % %     imp_2.Mx(lv(1), lv(2)) = sys.Mx(lv(1),lv(2));
% % %     imp_2.Mx(lv(2), lv(1)) = sys.Mx(lv(2),lv(1));

%   Calcoliamo il gap esatto (per poter usare dei Bx più grandi)
% % %     H = sys.H0 + sys.Bx * imp_2.Mx;
% % %     [~,d] = eig(H);
% % %     d = sort(diag(d));
% % %     d = diag(d(end:-1:1));
% % %     sys.H = d;
    imp_2.w = d(lv(1),lv(1)) - d(lv(2),lv(2));

% % %     imp_2.w = sys.w(lv(1),lv(2));

%%  Terzo Impulso   
    lv = [2 4];

    imp_3.lv      = lv;
    imp_3.w       = sys.w(lv(1),lv(2));
    imp_3.angolo  = abs(h(lv(1),lv(2)));
    imp_3.fase    = angle(h(lv(1),lv(2))) - angle(comm_zx(lv(1),lv(2)));
% % %     imp_3.angolo  = pi/2; %%%%%%%%%%%%%%% DBG
% % %     imp_3.fase    = pi/2; %%%%%%%%%%%%%%% DBG

% % % %   Rimozione completa leackage
% % %     imp_3.Mz               = zeros(dim, dim);
% % %     imp_3.Mz(lv(1), lv(1)) = sys.Mz(lv(1), lv(1));
% % %     imp_3.Mz(lv(2), lv(2)) = sys.Mz(lv(2), lv(2));

%   Rimozione parziale leackage
    imp_3.Mz    = zeros(dim, dim);
    imp_3.Mz    = diag(diag(sys.Mz));

% % %     imp_3.Mx               = zeros(dim, dim);
% % %     imp_3.Mx(lv(1), lv(2)) = sys.Mx(lv(1),lv(2));
% % %     imp_3.Mx(lv(2), lv(1)) = sys.Mx(lv(2),lv(1));

% % % %   Calcoliamo il gap esatto (per poter usare dei Bx più grandi)
% % %     H = sys.H0 + sys.Bx * imp_3.Mx;
% % %     [~,d] = eig(H);
% % %     d = sort(diag(d));
% % %     d = diag(d(end:-1:1));
% % %     sys.H = d;
    imp_3.w = d(lv(1),lv(1)) - d(lv(2),lv(2));

% % %     imp_3.w = sys.w(lv(1),lv(2));

%%  Quarto Impulso   
    lv = [3 4];

    imp_4.lv      = lv;
    imp_4.w       = sys.w(lv(1),lv(2));
    imp_4.angolo  = abs(h(lv(1),lv(2)));
    imp_4.fase    = angle(h(lv(1),lv(2))) - angle(comm_zx(lv(1),lv(2)));

% % % %   Rimozione completa leackage
% % %     imp_4.Mz               = zeros(dim, dim);
% % %     imp_4.Mz(lv(1), lv(1)) = sys.Mz(lv(1), lv(1));
% % %     imp_4.Mz(lv(2), lv(2)) = sys.Mz(lv(2), lv(2));

%   Rimozione parziale leackage
    imp_4.Mz    = zeros(dim, dim);
    imp_4.Mz    = diag(diag(sys.Mz));

% % %     imp_4.Mx               = zeros(dim, dim);
% % %     imp_4.Mx(lv(1), lv(2)) = sys.Mx(lv(1),lv(2));
% % %     imp_4.Mx(lv(2), lv(1)) = sys.Mx(lv(2),lv(1));

% % % %   Calcoliamo il gap esatto (per poter usare dei Bx più grandi)
% % %     H = sys.H0 + sys.Bx * imp_4.Mx;
% % %     [~,d] = eig(H);
% % %     d = sort(diag(d));
% % %     d = diag(d(end:-1:1));
% % %     sys.H = d;
    imp_4.w = d(lv(1),lv(1)) - d(lv(2),lv(2));

% % %     imp_4.w = sys.w(lv(1),lv(2));

%%  Settaggio ampiezze campi oscillanti
% % %     K_12 = abs(imp_1.angolo * imp_1.w / abs(comm_zx(1,2)) )
% % %     K_13 = abs(imp_2.angolo * imp_2.w / abs(comm_zx(1,3)) )
% % %     K_24 = abs(imp_3.angolo * imp_3.w / abs(comm_zx(2,4)) )
% % %     K_34 = abs(imp_4.angolo * imp_4.w / abs(comm_zx(3,4)) )
% % % 
% % %     B_max = 0.2
% % % 
% % % %   Settiamo le ampiezze dei campi magnetici per ottenere transizioni della
% % % %   stessa durata
% % %     imp_1.Bz = K_12 / K_13 * B_max;
% % %     imp_2.Bz = B_max;
% % %     imp_3.Bz = K_24 / K_13 * B_max;
% % %     imp_4.Bz = K_34 / K_13 * B_max;
% % % 
% % % % % % %   Controlliamo che i rapporti siano effettivamente uguali
% % %     K_12 / imp_1.Bz
% % %     K_13 / imp_2.Bz
% % %     K_24 / imp_3.Bz
% % %     K_34 / imp_4.Bz

%%  Settiamo le durate sfruttando la conoscenza di Bx (DBG)
    
%   Dobbiamo avere B_max > Bx
    B_max = 0.2;
    
    imp_1.Bz = B_max * imp_1.w * imp_1.angolo / abs(comm_zx(imp_1.lv(1),imp_1.lv(2)));
    imp_2.Bz = B_max * imp_2.w * imp_2.angolo / abs(comm_zx(imp_2.lv(1),imp_2.lv(2)));
    imp_3.Bz = B_max * imp_3.w * imp_3.angolo / abs(comm_zx(imp_3.lv(1),imp_3.lv(2)));
    imp_4.Bz = B_max * imp_4.w * imp_4.angolo / abs(comm_zx(imp_4.lv(1),imp_4.lv(2)));

% % % %   Il tempo lo da la rotazione più lenta (per individuare la transizione
% % % %   più lenta)
% % %     2 * K_12 / sys.Bx
% % %     2 * K_13 / sys.Bx
% % %     2 * K_24 / sys.Bx
% % %     2 * K_34 / sys.Bx

% % % %   Durata rotazione più lenta
% % %     T_gate = 2 * K_24 / (sys.Bx * B_max)
% % % 
% % % %   Settiamo in modo opportuno le ampiezze dei relativi campi
% % %     B_max = 0.2;
% % %     imp_1.Bz = B_max;
% % %     imp_2.Bz = B_max;
% % %     imp_3.Bz = B_max;
% % %     imp_4.Bz = B_max;
% % % 
% % % %   
% % %     T = 500;
% % %     imp_1.Bz = imp_1.Bz / imp_1.d_Sz;
% % %     imp_2.Bz = imp_2.Bz / imp_2.d_Sz;
% % %     imp_3.Bz = imp_3.Bz / imp_3.d_Sz;
% % %     imp_4.Bz = imp_4.Bz / imp_4.d_Sz;
% % %    
% % %     [imp_1.d_Mz * 621/T  imp_1.d_Sz * 621/T]
% % %     [imp_2.d_Mz * 658/T  imp_2.d_Sz * 658/T]
% % %     [imp_3.d_Mz * 177/T  imp_3.d_Sz * 177/T]
% % %     [imp_4.d_Mz * 396/T  imp_4.d_Sz * 396/T]

%%  Inizializzazione sequenza di impulsi
    seq_imp{1} = imp_1;
    seq_imp{2} = imp_2;
    seq_imp{3} = imp_3;
    seq_imp{4} = imp_4;
end