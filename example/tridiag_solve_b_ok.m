% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2011
% Authors: Bradlely Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
%          Aleksander Aravkin: saravkin at eos dot ubc dot ca
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin tridiag_solve_b_ok.m$$ $newlinech %$$
% $spell
%       ckbs
%       ak
%       blk
%       lhs
%       qk
%       tridiag
%       uk
%       blktridiag
%       mul
% $$
%
% $section ckbs_tridiag_solve_b Example and Test$$
%
% $index ckbs_tridiag_solve_b, example and test$$
% $index tridiag_solve_b, example and test$$
% $index example, tridiag_solve_b$$
% $index test, tridiag_solve_b$$
%
% $head Source Code$$
% $newlinech $$ $codep
function [ok] = tridiag_solve_b_ok()
    ok = true;
    m  = 2;
    n  = 3;
    N  = 4;
    % case where uk = 0, qk = I, and ak is random
    rand('seed', 123);
    a = rand(n, n, N);
    r = rand(n * N, m);
    c = zeros(n, n, N);
    for k = 1 : N
        ak         = a(:, :, k);
        b(:, :, k) = 2 * eye(n) + ak * ak';
        c(:, :, k) = ak';
    end
    % ----------------------------------------
    [e, lambda] = ckbs_tridiag_solve_b(b, c, r);
    % ----------------------------------------

    check = ckbs_blktridiag_mul(b, c, e);

    ok = ok & max(abs(check - r)) < 1e-10;

    return
end
% $$ $newlinech %$$
% $end
