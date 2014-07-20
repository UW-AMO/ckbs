% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2011
% Authors: Bradlely Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
%          Aleksanr Aravin:      saravkin at eos dot ubc dot ca
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin blkdiag_mul_ok.m$$ $newlinech %$$
% $spell
%       ckbs
%       Bdiag
%       blk
%       blkdiag
%       mul
% $$
%
% $section blkdiag_mul Example and Test$$
%
% $index ckbs_blkdiag_mul, example and test$$
% $index blkdiag_mul, example and test$$
% $index example, blkdiag_mul$$
% $index test, blkdiag_mul$$
%
% $head Source Code$$
% $newlinech $$ $codep
function [ok] = blkdiag_mul_ok()
    ok = true;
    % -------------------------------------------------------------
    % You can change these parameters
    m    = 2;
    n    = 3;
    N    = 2;
    k    = 3;
    % -------------------------------------------------------------
    % Define the problem
    rand('seed', 123);
    v     = rand(n * N, k);
    Bdiag = zeros(m, n, N);
    B     = zeros(m * N , n * N);
    blk_m = 1 : m;
    blk_n = 1 : n;
    for k = 1 : N
        Bdiag(:, :, k)  = rand(m, n);
        B(blk_m, blk_n) = Bdiag(:, :, k);
        blk_m           = blk_m + m;
        blk_n           = blk_n + n;
    end
    % -------------------------------------
    w     = ckbs_blkdiag_mul(Bdiag, v);
    % -------------------------------------
    check = B * v;
    ok    = ok & ( max(max(abs(w - check))) < 1e-10 );
    return
end
% $$ $newlinech %$$
% $end
