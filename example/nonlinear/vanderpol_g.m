% extended vanderpol 



function [gk, Gk] = vanderpol_g(k, xk1, params)
    W = params.W; % 2ndim x 2ndim cross talk matrix 
    ndim = size(W, 1)/2; % number of dimensions 
    
    initial = params.vanderpol_g_initial; % should be a vector of size 2n.
    dt      = params.vanderpol_g_dt;
    mu      = params.vanderpol_g_mu;

    
    % make a little block Gk and gk. 

    if k == 1        
        gk = initial;
        Gk = zeros(2*ndim, 2*ndim);
    else
        gk = zeros(2*ndim,1);
        Gk = zeros(2*ndim); 
        for s = 1:ndim
            % build little guy
            fI = 2*(s-1)+1;
            sI = 2*s;
            gkk      = zeros(2, 1);
            gkk(1)   = xk1(fI) + xk1(sI) * dt;
            gkk(2)   = xk1(sI) + (mu * (1 - xk1(fI)*xk1(fI) * xk1(sI) - xk1(fI))) * dt;
            %
            Gkk      = zeros(2, 2);
            Gkk(1,1) = 1;
            Gkk(1,2) = dt;
            Gkk(2,1) = (mu * (- 2 * xk1(fI)) * xk1(sI) - 1) * dt;
            Gkk(2,2) = 1 + mu * (1 - xk1(fI)*xk1(fI)) * dt;
            gk(fI) = gkk(1);
            gk(sI) = gkk(2);
            Gk(fI:sI, fI:sI) = Gkk;
        end
        Gk = Gk + W; % add in the model
    end
    return
end
% $$ $newlinech %$$
% $end
