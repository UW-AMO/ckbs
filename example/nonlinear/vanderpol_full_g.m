% extended vanderpol 


function [gk, Gk] = vanderpol_full_g(k, xk1, params)
    W = params.W; % 2ndim x 2ndim cross talk matrix 
    ndim = size(W, 1)/2; % number of dimensions 
    
    initial = params.vanderpol_g_initial; % should be a vector of size 2n.
    dt      = params.vanderpol_g_dt;
    mu      = params.vanderpol_g_mu;

   a1 = params.a1; 
    a2 = params.a2;
    a3 = params.a3;
    a4 = params.a4; 
    a5 = params.a5;
    
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
            Gkk      = zeros(2, 2);

            gkk(1)   = xk1(fI) -a3* xk1(sI) * dt - a4*(xk1(fI) - a5)*dt;
            Gkk(1,1) = 1 - a4*dt;
            Gkk(1,2) = -a3*dt;

            
            gkk(2)   = xk1(sI) + (a1 * xk1(sI)*(1 - xk1(sI)*xk1(sI))  + a2*xk1(fI)) * dt;
            Gkk(2,1) =  a2*dt;
            Gkk(2,2) = 1 + a1 * (1 - xk1(sI)*xk1(sI)) * dt + (a1 * (- 2 * xk1(sI)) * xk1(sI) - 1)*dt;

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
