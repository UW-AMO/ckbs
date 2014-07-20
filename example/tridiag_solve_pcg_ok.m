% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2011
% Authors: Bradlely Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
%          Aleksandr Aravkin:    saravkin at eos dot ubc dot ca
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin tridiag_solve_pcg_ok.m$$ $newlinech %$$
% $spell
%       pcg
%       eps
%       randn
%       relres
%       ckbs
%       ak
%       iter
%       cg
%       blk
%       lhs
%       qk
%       tridiag
%       uk
%       blktridiag
%       mul
% $$
%
% $section ckbs_tridiag_solve_pcg Example and Test$$
%
% $index ckbs_tridiag_solve_pcg, example and test$$
% $index tridiag_solve_pcg, example and test$$
% $index example, tridiag_solve_pcg$$
% $index test, tridiag_solve_pcg$$
%
% $head Source Code$$
% $newlinech $$ $codep
function [ok] = tridiag_solve_pcg_ok()
    ok = true;
    m  = 1;
    n  = 50;
    N  = 40;
    eps = 1e-5;
    % case where uk = 0, qk = I, and ak is random
    %rand('seed', 123);
    a = rand(n, n, N);
    r = randn(n * N, m);
    c = zeros(n, n, N);
    for k = 1 : N
        ak         = a(:, :, k);
        b(:, :, k) = 5 * eye(n) + ak * ak';
        c(:, :, k) = ak';
    end
    % ----------------------------------------
    [e, flag, relres, iter] = ckbs_tridiag_solve_pcg(b, c, r);
    % ----------------------------------------
    %
    check = ckbs_blktridiag_mul(b, c, e);
    max(abs(check - r))
    ok = ok & max(abs(check - r))/max(abs(r)) < eps;

    return
end
% $$ $newlinech %$$
% $end
