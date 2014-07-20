% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2011
% Authors: Bradlely Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
%          Aleksandr Aravkin:    saravkin at eos dot ubc dot ca
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin kuhn_tucker_L1_ok.m$$ $newlinech %$$
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
%       Vec
%       Hy
%       blktridiag
%       sqrt
%       Bt
%       Sm
% $$
%
% $section ckbs_kuhn_tucker_L1 Example and Test$$
%
% $index ckbs_kuhn_tucker_L1, example and test$$
% $index kuhn_tucker_L1, example and test$$
% $index example, kuhn_tucker_L1$$
% $index test, kuhn_tucker_L1$$
%
% $head Source Code$$
% $newlinech $$ $codep
function [ok] = kuhn_tucker_L1_ok()
    ok = true;
    % -------------------------------------------------------------
    % You can change these parameters
    m    = 5;
    n    = 4;
    N    = 3;
    % -------------------------------------------------------------
    rand('seed', 123);
    %r     = m * N;
    mu    = 1.5;
    s     = rand(m, N);
    y     = rand(n, N);
    r     = rand(m, N);
    b     = rand(m, N);
    d     = rand(n, N);
    Bdia  = rand(m, n, N);
    Hdia  = rand(n, n, N);
    Hlow  = rand(n, n, N);
    pPlus = rand(m, N);
    pMinus = rand(m, N);


    % -------------------------------------
    F =ckbs_kuhn_tucker_L1(mu, s, y, r, b, d, Bdia, Hdia, Hlow, pPlus, pMinus);
    % -------------------------------------

    yVec = y(:);
    dVec = d(:);
    pPlusVec   = pPlus(:);
    pMinusVec  = pMinus(:);
    rVec = r(:);
    sVec = s(:);
    bVec = b(:);
    Hy = ckbs_blktridiag_mul(Hdia, Hlow, y(:));
    Bt_SmR = ckbs_blkdiag_mul_t(Bdia, s(:)-r(:));
    By = ckbs_blkdiag_mul(Bdia, y(:));

    check = [ ...
        pPlusVec - pMinusVec - bVec - By; ...
        pMinusVec.*sVec - mu;...
        rVec + sVec - 2*sqrt(2); ...
        pPlusVec.*rVec - mu; ...
        Hy + dVec + 0.5*Bt_SmR; ...
        ];
    ok    = ok & ( max(abs(F - check)) < 1e-10 );
    return
end
% $$ $newlinech %$$
% $end
