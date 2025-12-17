function fai_grafico(errore_medio, T2_vals, tol)

%     errore_medio = pulisci_matrice(errore_medio.', tol);
    
    % Facciamo il grafico dell'errore osservato
    loglog(1./T2_vals, errore_medio, 'LineWidth', 3);
    grid on
    hold on

end