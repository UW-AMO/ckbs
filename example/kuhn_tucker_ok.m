% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2011
% Authors: Bradlely Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
%          Aleksandr Aravkin:    saravkin at eos dot ubc dot ca
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin kuhn_tucker_ok.m$$ $newlinech %$$
% $spell
%       ckbs
%       Bdiag
%       blk
%       blkdiag
%       mul
%       Kuhn
%       mu
%       Bdia
%       Hdia
%       Hlow
%       Btu
%       Hy
%       blktridiag
% $$
%
% $section ckbs_kuhn_tucker Example and Test$$
%
% $index ckbs_kuhn_tucker, example and test$$
% $index kuhn_tucker, example and test$$
% $index example, kuhn_tucker$$
% $index test, kuhn_tucker$$
%
% $head Source Code$$
% $newlinech $$ $codep
function [ok] = kuhn_tucker_ok()
    ok = true;
    % -------------------------------------------------------------
    % You can change these parameters
    m    = 5;
    n    = 4;
    N    = 3;
    % -------------------------------------------------------------
    rand('seed', 123);
    p     = n * N;
    r     = m * N;
    mu    = 1.5;
    s     = rand(r, 1);
    y     = rand(p, 1);
    u     = rand(r, 1);
    b     = rand(r, 1);
    d     = rand(p, 1);
    Bdia  = rand(m, n, N);
    Hdia  = rand(n, n, N);
    Hlow  = rand(n, n, N);

    By = ckbs_blkdiag_mul(Bdia, y);
    Btu = ckbs_blkdiag_mul_t(Bdia, u);
    Hy = ckbs_blktridiag_mul(Hdia, Hlow, y);

    % -------------------------------------
    F = ckbs_kuhn_tucker(mu, s, y, u, b, d, Bdia, Hdia, Hlow);
    % -------------------------------------
    check = [ ...
        s + b + By; ...
        Hy + Btu + d; ...
        s .* u - mu ...
        ];
    ok    = ok & ( max(abs(F - check)) < 1e-10 );
    return
end
% $$ $newlinech %$$
% $end
