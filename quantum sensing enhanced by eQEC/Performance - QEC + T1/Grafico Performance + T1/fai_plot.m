function fai_plot(x, y, spin_vals, flag, label)
%   Funzione che genera il plot delle curve (x,y), (x,y+3*std_dev), (x,
%   y-3*std_dev) e (y,x), (y+3*std_dev,x), (y-3*std_dev, x) in scala
%   logaritmica

%   Dimensioni dati
    sz_y = size(y);
    sz_x = size(x);
    n_spin = length(spin_vals);
    n_Bx   = length(x(1,:));
    dims   = 2*spin_vals + 1;

    for i = 1:n_spin
        if(flag)
            lab_leg{i+1} = strcat(string(dims(i)-1),'/2');
        else
            lab_leg{i} = strcat(string(dims(i)-1),'/2');
        end
    end

%%  Primo Plot

%   Recuperiamo il numero di dati da plottare
    figure()
    hold on
    grid on
    set(gca,'YScale','log','XScale','log','FontSize',15);

    if (flag)
        loglog(x,x, 'LineWidth',3);
        lab_leg{1} = "Bx";
    end


    for i = 1:n_spin
        loglog(reshape(x(i,:),1,n_Bx), reshape(y(i,:),1,n_Bx), '-','LineWidth', 2)
    end

    xlabel(label{1},"FontSize",20, 'Interpreter', 'latex');
    ylabel(label{2},"FontSize",20,'Rotation',90, 'Interpreter', 'latex');

    legend(lab_leg,'Location','west');

    title(strcat('Intervallo QEC $',label{7},'$ ns'), 'Interpreter', 'latex');

%     saveas(gca, strcat('./images/',label{5},'_',label{6},'_',label{3},label{4},'.fig'));
%     saveas(gca, strcat('./images/',label{5},'_',label{6},'_',label{3},label{4},'.svg'));
    
end