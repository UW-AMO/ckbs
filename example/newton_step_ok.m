% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2011
% Authors: Bradlely Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
%          Aleksandr Aravkin:    saravkin at eos dot ubc dot ca
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin newton_step_ok.m$$ $newlinech %$$
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
%       mul
%       speye
%       kuhn
%       Btu
%       Hy
%       blktridiag
%       blkdiag
% $$
%
% $section ckbs_newton_step Example and Test$$
%
% $index ckbs_newton_step, example and test$$
% $index newton_step, example and test$$
% $index example, newton_step$$
% $index test, newton_step$$
%
% $head Source Code$$
% $newlinech $$ $codep
function [ok] = newton_step_ok()
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
    r     = m * N;
    p     = n * N;
    s     = rand(r, 1) + 1;
    y     = rand(p, 1);
    u     = rand(r, 1) + 1;
    b     = rand(r, 1);
    d     = rand(p, 1);
    H     = zeros(p, p);
    B     = zeros(r, r);
    Hdia  = zeros(n, n, N);
    Hlow  = zeros(n, n, N);
    Bdia  = rand(m, n, N);
    for k = N : -1 : 1
        tmp             = rand(n, n);
        Hlow(:,:, k)    = tmp;
        Hdia(:,:, k)   = (tmp * tmp') + 4 * eye(n);
    end

    % Form full matrices
    B = ckbs_blkdiag_mul(Bdia, speye(p));
    H = ckbs_blktridiag_mul(Hdia, Hlow, speye(p));

    % -------------------------------------------------------------------
    [ds, dy, du] = ckbs_newton_step(mu, s, y, u, b, d, Bdia, Hdia, Hlow);
    % -------------------------------------------------------------------
    F      = ckbs_kuhn_tucker(mu, s, y, u, b, d, Bdia, Hdia, Hlow);


    dF     = [ ...
        eye(r)      , B           , zeros(r, r)
        zeros(p, r) , H           , B'
        diag(u)     , zeros(r, p) , diag(s) ...
        ];
    delta        = [ ds ; dy ; du ];
    ok           = ok & ( max( abs( F + dF * delta ) ) <= 1e-10 );
    return
end
% $$ $newlinech %$$
% $end
