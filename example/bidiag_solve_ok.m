% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2011
% Authors: Bradlely Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
%          Aleksander Aravkin: saravkin at eos dot ubc dot ca
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin bidiag_solve_ok.m$$ $newlinech %$$
% $spell
%       ckbs
%       ak
%       blk
%       lhs
%       qk
%       bidiag
%       uk
%       blkbidiag
%       mul
% $$
%
% $section ckbs_bidiag_solve Example and Test$$
%
% $index ckbs_bidiag_solve, example and test$$
% $index bidiag_solve, example and test$$
% $index example, bidiag_solve$$
% $index test, bidiag_solve$$
%
% $head Source Code$$
% $newlinech $$ $codep
function [ok] = bidiag_solve_ok()
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
    [e, lambda] = ckbs_bidiag_solve(b, c, r);
    % ----------------------------------------

    check = ckbs_blkbidiag_mul(b, c, e);

    ok = ok & max(max(abs(check - r))) < 1e-10;

    return
end
% $$ $newlinech %$$
% $end
