function grafico_evoluzione(sim)
%   Funzione che prende in input una serie di matrici di densità e ne
%   fa il grafico dell'andamento temporale delle popolazioni

%   Recuperiamo la dimensione del sistema
    sz = size(sim.evol_L);
    dim = sz(1);

    ev = zeros(1,sz(3));

    time = sim.time * sim.sys.h_bar;

%   Stabiliamo se plottare tutti gli stati logici o solo quelli con k=0
    passo = 1;
%     passo = dim/2;

% % %     figure()
% % %     hold on
% % %     grid on
% % %     for i = 1:dim
% % %         for j = 1:dim
% % %             if (i~=j)
% % %                 for k = 1:sz(3)
% % %                    ev(k) = atan(real(sim.evol(i,j,k) / imag(sim.evol(i,j,k))));
% % %                 end
% % %             plot([0 sim.time], ev, 'LineWidth', 3);
% % %             end
% % %         end
% % %     end
% % % 
% % %     figure()
% % %     hold on
% % %     grid on
% % %     for i = 1:dim
% % %         for j = 1:dim
% % %             if (i~=j)
% % %                 for k = 1:sz(3)
% % %                    ev(k) = real(sim.evol(i,k));
% % %                 end
% % %             plot([0 time], ev, 'LineWidth', 3);
% % %             end
% % %         end
% % %     end
% % % 
% % %     figure()
% % %     hold on
% % %     grid on
% % %     for i = 1:dim
% % %         for j = 1:dim
% % %             if (i~=j)
% % %                 for k = 1:sz(3)
% % %                    ev(k) = imag(sim.evol(i,k));
% % %                 end
% % %                 plot([0 sim.time], ev, 'LineWidth', 3);
% % %             end
% % %         end
% % %     end

    figure()
    hold on
    grid on
    for i = 1:dim
            for k = 1:sz(3)
               ev(k) = angle(sim.evol(i,k));
            end
        plot([0 time], ev, 'LineWidth', 3);
    end

    figure()
    hold on
    grid on
    for i = 1:passo:dim
        for j = 1:sz(3)
            ev(j) = abs(sim.evol(i,j))^2;
        end
        plot([0 time], ev,'LineWidth',3);
    end

    figure()
    for i = 1:passo:dim
        for j = 1:sz(3)
            ev(j) = sim.evol_L(i,i,j);
        end
        plot([0 time], ev,'LineWidth',3);
        grid on
        hold on
    end
end