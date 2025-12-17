function [sim, freq] = calcola_frequenza_oscillazione(sys, Bx)
%   Funzione che simula un un'oscillazione logica dovuta al campo esterno
%   Bx e ne ricava la frequenza individuata dal primo minimo
%   dell'oscillazione

%   Calcoliamo tutti i parametri della simulazine da eseguire
    [sys, seq_imp] = genera_sequenza_impulsi(sys);
%     [sys, seq_imp] = genera_sequenza_impulsi_MOD(sys);
    
%   Inizializziamo la simulazione da eseguire
    sim.sys = sys;
    sim.imp = seq_imp;
    sim.Bx  = Bx;

%   Approssimiamo numericamente l'evoluzione del sistema durante gli
%   impulsi oscillanti lungo Z e calcoliamo la frequenza di oscillazione
%   logica
    [sim, freq] = simula_impulsi(sim);

%   Facciamo il plot dell'evoluzione temporale dei vari stati logici
    grafico_evoluzione(sim); 
end