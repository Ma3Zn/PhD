# =====================================================
# main.jl
#
# Script principale che:
#  1) Carica i parametri principali (N, chi, parametri del circuito, K_list),
#  2) Crea lo stato MPS iniziale,
#  3) Applica il circuito (dipendente da un singolo parametro) sull’MPS,
#  4) Costruisce l’Hamiltoniana di Heisenberg (MPO),
#  5) Calcola iterativamente <ψ|H^k|ψ> per k in K_list,
#  6) Stampa/Salva i risultati.
# =====================================================

include("Circuit.jl")
using .Circuit

include("Utils.jl")
using .Utils

 # -------------------------------------------------
 # 1) Definizione dei parametri
 # -------------------------------------------------
 N_qubit = 8            # numero di qubit
 chi = 100              # massima bond dimension
 cutoff = 1e-8          # magnitudo massima valori singolari mantenuti nella SVD
 N_tstep = 5            # Numero di trotter step da utilizzare nell'evoluzione
 K_list = [1, 2, 3]     # potenze dell'Hamiltoniana da misurare per la SbE

 # -------------------------------------------------
 # 2) Creazione dell'MPS iniziale nello stato opportuno (|0101...01>)
 # -------------------------------------------------
 ψ = create_initial_mps(N_qubit)

 # -------------------------------------------------
 # 3) Evoluzione del MPS iniziale tramite la 
 #    digitalizzazione del processo di annealing
 #    utilizzando N_tstep
 # -------------------------------------------------
 apply_trotterized_annealing!(ψ, N_qubit, N_tstep, J, 
                              drive, lam, drive_sp, lam_sp, chi;
                              cutoff)

 # # -------------------------------------------------
 # # 4) Costruzione dell'Hamiltoniana di Heisenberg (MPO)
 # #    con periodic boundary conditions
 # # -------------------------------------------------
 # heisenberg_mpo = build_heisenberg_mpo(N)

 # # -------------------------------------------------
 # # 5) Calcolo di <ψ|H^k|ψ>, per k in K_list
 # #    (tramite moltiplicazioni ripetute di H su ψ)
 # # -------------------------------------------------
 # results = Dict{Int, Float64}()
 # for k in K_list
 #     # Copia di final_mps su cui iterare l'applicazione di H
 #     psi_temp = copy(final_mps)
 #     for i in 1:k
 #         psi_temp = apply_mpo(heisenberg_mpo, psi_temp, chi)
 #     end
 #     # Calcolo dell'overlap con lo stato originale
 #     results[k] = inner(final_mps, psi_temp)
 # end

 # # -------------------------------------------------
 # # 6) Output dei risultati
 # # -------------------------------------------------
 # output_results(results)