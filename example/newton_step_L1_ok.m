% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2011
% Authors: Bradley Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
%          Aleksandr Aravkin:    saravkin at eos dot ubc dot ca
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin newton_step_L1_ok.m$$ $newlinech %$$
% $spell
%       ckbs
%       Bdia
%       blk
%       ds
%       du
%       dy
%       end end
%       hdia
%       hlow
%       mu
%       tmp
%       blktridiag
%       mul
%       speye
%       blkdiag
%       kuhn
% $$
%
% $section ckbs_newton_step_L1 Example and Test$$
%
% $index ckbs_newton_step_L1, example and test$$
% $index newton_step_L1, example and test$$
% $index example, newton_step_L1$$
% $index test, newton_step_L1$$
%
% $head Source Code$$
% $newlinech $$ $codep
function [ok] = newton_step_L1_ok()
    ok = true;
    % --------------------------------------------------------
    % You can change these parameters
    m    = 2;   % number of measurements per time point
    n    = 3;   % number of state vector components per time point
    N    = 4;   % number of time points
    % ---------------------------------------------------------
    %  Define the problem
    rand('seed', 123);
    mu    = .5;
    s     = rand(m, N);
    y     = rand(n, N);
    r     = rand(m, N);
    b     = rand(m, N);
    d     = rand(n, N);
    Hdia  = zeros(n, n, N);
    Hlow  = zeros(n, n, N);
    Bdia  = zeros(m, n, N);
    pPlus = rand(m, N);
    pMinus = rand(m, N);
    for k = N : -1 : 1
        Bdia(:,:, k)    = rand(m, n);
        tmp             = rand(n, n);
        Hlow(:,:, k)    = tmp;
        Hdia(:,:, k)   = (tmp * tmp') + 4 * eye(n);
        %
    end

    % Form full matrices.
    H = ckbs_blktridiag_mul(Hdia, Hlow, speye(n*N));
    B = ckbs_blkdiag_mul(Bdia, speye(n*N));

    % -------------------------------------------------------------------
    [dPPlus, dPMinus, dR, dS, dY] = ckbs_newton_step_L1(mu, s, y, r, b, d, Bdia, Hdia, Hlow, pPlus, pMinus);
    % -------------------------------------------------------------------
    F      =ckbs_kuhn_tucker_L1(mu, s, y, r, b, d, Bdia, Hdia, Hlow, ...
        pPlus, pMinus);

    dF     = [ ...
        speye(m*N) , -speye(m*N), zeros(m*N), zeros(m*N), -B
        zeros(m*N) , diag(s(:)) , zeros(m*N), diag(pMinus(:)),zeros(m*N,n*N)
        zeros(m*N) , zeros(m*N) , speye(m*N), speye(m*N),zeros(m*N,n*N)
        diag(r(:)) , zeros(m*N) , diag(pPlus(:)), zeros(m*N),zeros(m*N,n*N)
        zeros(n*N, m*N), zeros(n*N, m*N), -B'/2, B'/2, H
        ];
    delta        = [ dPPlus(:); dPMinus(:); dR(:); dS(:); dY(:) ];
    ok           = ok & ( max( abs( F + dF * delta ) ) <= 1e-10 );
    return
end
% $$ $newlinech %$$
% $end
