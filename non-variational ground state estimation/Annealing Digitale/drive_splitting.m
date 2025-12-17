function f = drive_splitting(t, par)
%   Funzione che genera il drive del termine di splitting
    
% % % %   Funziona bene in combo con un drive lineare
% % %     f = 0;
% % % 
% % %     if ((t * par.lambda) ^ par.n < 1)
% % %         f = 1 - (par.lambda * t) ^ par.n;
% % %     end

%   Funziona bene in combo con un drive tipo sin
    f = tan(1) - tan((par.lambda*t)^par.n);

% % % %   Tipo la tangente ma peggio
% % %     f = cos((par.lambda*t)^par.n*pi/2);
end