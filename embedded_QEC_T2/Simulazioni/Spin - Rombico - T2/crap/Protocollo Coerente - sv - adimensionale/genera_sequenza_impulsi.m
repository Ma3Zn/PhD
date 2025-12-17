function [sys, seq_imp] = genera_sequenza_impulsi(sys)
%   Funzione che genera la sequenza di impulsi lungo z da mandare al
%   sistema per renderne degeneri gli opportuni livelli energetici in modo
%   da riuscire ad eseguire una rotazione logica sotto l'effetto del campo
%   statico incognito Bx

%   Recuperiamo la dimensione del sistema
    dim = sys.dim;

% % % %   Calcoliamo il gap esatto (per poter usare dei Bx più grandi)
% % %     H = sys.H0 + sys.Bx * sys.Mx;
% % %     [V,d] = eig(H);
% % %     d = sort(diag(d));
% % %     d = diag(d(end:-1:1));
% % %     sys.H = d;

% % % %   Non capisco se sia necessario o meno
% % %     sys.Mz = V' * sys.Mz * V;

%%  Gate logico da eseguire

%   Recuperiamo il generatore di tale gate
    h = 1i*logm(sys.GL);

%   Riscriviamo il generatore nella base computazionale
    h = sys.CB * h * sys.CB';

%%  Calcoliamo il commutatore di Mx ed Mz
    comm_zx = diag(diag(sys.Mz))*sys.Mx - sys.Mx*diag(diag(sys.Mz));

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

% % % %   Calcoliamo il gap esatto (per poter usare dei Bx più grandi)
% % %     imp_1.w = d(lv(1),lv(1)) - d(lv(2),lv(2));
    
% % %     imp_1.w = sys.w(lv(1),lv(2));

%%  Secondo Impulso   
    lv = [1 3];

    imp_2.lv      = lv;
    imp_2.w       = sys.w(lv(1),lv(2));
    imp_2.angolo  = abs(h(lv(1),lv(2)));
    imp_2.fase    = angle(h(lv(1),lv(2))) - angle(comm_zx(lv(1),lv(2)));

%   Calcoliamo il gap esatto (per poter usare dei Bx più grandi)
% % %     imp_2.w = d(lv(1),lv(1)) - d(lv(2),lv(2));

% % %     imp_2.w = sys.w(lv(1),lv(2));

%%  Terzo Impulso   
    lv = [2 4];

    imp_3.lv      = lv;
    imp_3.w       = sys.w(lv(1),lv(2));
    imp_3.angolo  = abs(h(lv(1),lv(2)));
    imp_3.fase    = angle(h(lv(1),lv(2))) - angle(comm_zx(lv(1),lv(2)));

% % % %   Calcoliamo il gap esatto (per poter usare dei Bx più grandi)
% % %     imp_3.w = d(lv(1),lv(1)) - d(lv(2),lv(2));

% % %     imp_3.w = sys.w(lv(1),lv(2));

%%  Quarto Impulso   
    lv = [3 4];

    imp_4.lv      = lv;
    imp_4.w       = sys.w(lv(1),lv(2));
    imp_4.angolo  = abs(h(lv(1),lv(2)));
    imp_4.fase    = angle(h(lv(1),lv(2))) - angle(comm_zx(lv(1),lv(2)));

% % % %   Calcoliamo il gap esatto (per poter usare dei Bx più grandi)
% % %     imp_4.w = d(lv(1),lv(1)) - d(lv(2),lv(2));

% % %     imp_4.w = sys.w(lv(1),lv(2));

%%  Settaggio ampiezze campi oscillanti
    K_12 = abs(imp_1.angolo * imp_1.w / abs(comm_zx(1,2)) )
    K_13 = abs(imp_2.angolo * imp_2.w / abs(comm_zx(1,3)) )
    K_24 = abs(imp_3.angolo * imp_3.w / abs(comm_zx(2,4)) )
    K_34 = abs(imp_4.angolo * imp_4.w / abs(comm_zx(3,4)) )

    B_max = 0.02

%   Settiamo le ampiezze dei campi magnetici per ottenere transizioni della
%   stessa durata
    imp_1.Bz = K_12 / K_13 * B_max;
    imp_2.Bz = B_max;
    imp_3.Bz = K_24 / K_13 * B_max;
    imp_4.Bz = K_34 / K_13 * B_max;

% % % %   Controlliamo che i rapporti siano effettivamente uguali
    K_12 / imp_1.Bz
    K_13 / imp_2.Bz
    K_24 / imp_3.Bz
    K_34 / imp_4.Bz

%%  Settiamo le ampiezza come fa Luca (DBG)
    
% % % %   Dobbiamo avere B_max > Bx
% % %     B_max = 0.02;
% % %     
% % %     imp_1.Bz = B_max * imp_1.w * imp_1.angolo / abs(comm_zx(imp_1.lv(1),imp_1.lv(2)))
% % %     imp_2.Bz = B_max * imp_2.w * imp_2.angolo / abs(comm_zx(imp_2.lv(1),imp_2.lv(2)))
% % %     imp_3.Bz = B_max * imp_3.w * imp_3.angolo / abs(comm_zx(imp_3.lv(1),imp_3.lv(2)))
% % %     imp_4.Bz = B_max * imp_4.w * imp_4.angolo / abs(comm_zx(imp_4.lv(1),imp_4.lv(2)))

%%  Inizializzazione sequenza di impulsi
    seq_imp{1} = imp_1;
    seq_imp{2} = imp_2;
    seq_imp{3} = imp_3;
    seq_imp{4} = imp_4;
end