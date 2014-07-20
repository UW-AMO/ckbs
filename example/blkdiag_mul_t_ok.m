% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2011
% Authors: Bradlely Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
%          Aleksandr Aravkin:    saravkin at eos dot ubc dot ca
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin blkdiag_mul_t_ok.m$$ $newlinech %$$
% $spell
%       ckbs
%       Bdiag
%       blk
%       blkdiag
%       mul
% $$
%
% $section blkdiag_mul_t Example and Test$$
%
% $index ckbs_blkdiag_mul_t, example and test$$
% $index blkdiag_mul_t, example and test$$
% $index example, blkdiag_mul_t$$
% $index test, blkdiag_mul_t%$$
%
% $head Source Code$$
% $newlinech $$ $codep
function [ok] = blkdiag_mul_t_ok()
    ok = true;
    % -------------------------------------------------------------
    % You can change these parameters
    m    = 2;
    n    = 3;
    N    = 2;
    p   =  5;
    % -------------------------------------------------------------
    % Define the problem
    rand('seed', 123);
    v     = rand(m * N, p);
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
    w     = ckbs_blkdiag_mul_t(Bdiag, v);
    % -------------------------------------
    check = B' * v;
    ok    = ok & ( max(max(abs(w - check))) < 1e-10 );
    return
end
% $$ $newlinech %$$
% $end
