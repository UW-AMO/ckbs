% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-11
% Authors: Bradlely Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
%          Aleksandr Aravkin:    saravkin at eos dot ubc dot ca
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin blkbidiag_mul_ok.m$$ $newlinech %$$
% $spell
%       ckbs
%       Bdiag
%       Boffdiag
%       blk
%       blkdiag
%       blkbidiag
%       mul
%       bidiag
%       ak
%       Hdia
%       Hlow
%       Hcheck
%       speye
% $$
%
% $section blkbidiag_mul Example and Test$$
%
% $index ckbs_blkbidiag_mul, example and test$$
% $index blkbidiag_mul, example and test$$
% $index example, blkbidiag_mul$$
% $index test, blkbidiag_mul$$
%
% $head Source Code$$
% $newlinech $$ $codep
function [ok] = blkbidiag_mul_ok()
    ok = true;
    % -------------------------------------------------------------
    % You can change these parameters
    m    = 2;
    n    = 3;
    N    = 4;
    p   = 3;
    % -------------------------------------------------------------
    % Define the problem
    rand('seed', 123);
    v     = rand(n * N, 3);

    rand('seed', 123);
    a = rand(n, n, N);
    r = rand(n * N, m);
    Hdia  = rand(n, n, N);
    Hlow  = rand(n, n, N);
    Hcheck = zeros(n*N);

    blk_n = 1 : n;
    for k = 1 : N
        Hcheck(blk_n, blk_n) = Hdia(:,:, k);
        blk_n           = blk_n + n;
    end
    blk_n = 1:n;
    for k = 2:N;
        blk_minus       = blk_n;
        blk_n           = blk_n + n;
        Hcheck(blk_n, blk_minus) = Hlow(:,:, k);
       end
    % -------------------------------------

    H = ckbs_blkbidiag_mul(Hdia, Hlow, eye(n*N));

    % -------------------------------------
    check = H - Hcheck;
    ok    = ok & ( max(max(abs(check))) < 1e-10 );
    return
end
% $$ $newlinech %$$
% $end
