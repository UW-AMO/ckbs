% General van der Pol simulator (includes number of components, and
% cross-talk matrix. 
% ------------------------------------------------------------------------
function [x_out] = vanderpol_sim_full(params, xi, n_out, step, ndim, W)
    % assume xi is a 2*ndim x 1 vector, and x_out is 2*ndom x n_out vector
    % W is the cross-talk matrix 
    
    a1 = params.a1; 
    a2 = params.a2;
    a3 = params.a3;
    a4 = params.a4; 
    a5 = params.a5;


    
    % anonymous function for the ODE
    ode        = @(x) [ -a3*x(2) - a4*x(1) + a4*a5 ; a1 * ( 1 - x(2)*x(2) ) * x(2) + a2*x(1) ];
%    ode        = @(x) [ x(2) ; mu * ( 1 - x(1)*x(1) ) * x(2) - x(1) ];

    x_out      = zeros(2*ndim, n_out);
    x_out(:,1) = xi; 
    step_2     = step / 2;
    step_6     = step / 6;
    
    f1 = zeros(2*ndim,1);
    x2 = f1; 
    f2 = f1;
    x3 = f1;
    f3 = f1;
    x4 = f1;
    f4 = f1;
   
    for k = 2 : n_out
        x1 = x_out(:,k-1);
        for j = 1:ndim
            indsj = (j-1)*2+1:j*2;
            f1(indsj) = ode(x1(indsj)) + W(indsj, indsj)*x1(indsj);
            x2(indsj) = x1(indsj) + step_2 * f1(indsj);
            f2(indsj) = ode(x2(indsj)) + W(indsj, indsj)*x2(indsj);
            x3(indsj) = x1(indsj) + step_2 * f2(indsj);
            f3(indsj) = ode(x3(indsj)) + W(indsj, indsj)*x3(indsj);
            x4(indsj) = x1(indsj) + step * f3(indsj);
            f4(indsj) = ode(x4(indsj));
        end
        x_out(:,k) = x1 + step_6 * (f1 + 2 * f2 + 2 * f3 + f4);
    end
    return
end
% $$ $newlinech %$$
% $end
