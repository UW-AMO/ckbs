% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2011
% Authors: Bradlely Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
%          Aleksander Aravkin: saravkin at eos dot ubc dot ca
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin diag_solve_ok.m$$ $newlinech %$$
% $spell
%       ckbs
%       ak
%       blk
%       lhs
%       qk
%       diag
%       uk
%       blkdiag
%       mul
% $$
%
% $section ckbs_diag_solve Example and Test$$
%
% $index ckbs_diag_solve, example and test$$
% $index diag_solve, example and test$$
% $index example, diag_solve$$
% $index test, diag_solve$$
%
% $head Source Code$$
% $newlinech $$ $codep
function [ok] = diag_solve_ok()
    ok = true;
    m  = 2;
    n  = 3;
    N  = 4;
    % case where uk = 0, qk = I, and ak is random
    rand('seed', 123);
    a = rand(n, n, N);
    r = rand(n * N, m);
    for k = 1 : N
        ak         = a(:, :, k);
        b(:, :, k) = 2 * eye(n) + ak * ak';
    end
    % ----------------------------------------
    [e, lambda] = ckbs_diag_solve(b, r);
    % ----------------------------------------

    check = ckbs_blkdiag_mul(b, e);

    ok = ok & max(max(abs(check - r))) < 1e-10;

    return
end
% $$ $newlinech %$$
% $end
