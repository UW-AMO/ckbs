% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2011
% Authors: Bradlely Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
%          Aleksandr Aravkin:    saravkin at eos dot ubc dot ca
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin tridiag_solve_ok.m$$ $newlinech %$$
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
% $section ckbs_tridiag_solve Example and Test$$
%
% $index ckbs_tridiag_solve, example and test$$
% $index tridiag_solve, example and test$$
% $index example, tridiag_solve$$
% $index test, tridiag_solve$$
%
% $head Source Code$$
% $newlinech $$ $codep
function [ok] = tridiag_solve_ok()
    ok = true;
    m  = 1;
    n  = 1;
    N  = 5;
    % case where uk = 0, qk = I, and ak is random
    rand('seed', 123);
    a = rand(n, n, N);
    r = rand(n * N, m);
    c = zeros(n, n, N);
    b = zeros(n, n, N);
    for k = 1 : N
        ak         = a(:, :, k);
%        fprintf('ak is %5.3f\n', ak);
        b(:, :, k) = 2 * eye(n) + ak * ak';
        c(:, :, k) = ak';
    end
    
    % eigs(fullMat)
    %fullMat = ckbs_blktridiag_mul(b, c, eye(N*n))
    
   
    % ----------------------------------------
    [ef, lambda] = ckbs_tridiag_solve(b, c, r);
    [eb, lambda] = ckbs_tridiag_solve_b(b, c, r);
    [emf, lambda] = ckbs_tridiag_solve_mf(b, c, r);
    [enew, lambda] = ckbs_tridiag_solve_new(b, c, r);
    % ----------------------------------------
    %
    check_f = ckbs_blktridiag_mul(b, c, ef)-r;
    check_b = ckbs_blktridiag_mul(b, c, eb)-r;
    check_mf = ckbs_blktridiag_mul(b, c, emf)-r;
    check_new = ckbs_blktridiag_mul(b, c, enew)-r;
    
    ok = ok && max(max(abs(check_f))) < 1e-10;
    ok = ok && max(max(abs(check_b))) < 1e-10;
    ok = ok && max(max(abs(check_mf))) < 1e-10;
    ok = ok && max(max(abs(check_new))) < 1e-10;
    
    if(~ok)
       fprintf('Check_f: %5.3f, Check_b: %5.3f, Check_mf: %5.3f, Check_new: %5.3f\n', norm(check_f(:)), norm(check_b(:)), ...
       norm(check_mf(:)), norm(check_new(:)));
       %check_mf
       [ef eb emf enew]
       
    end
    
    return
end
% $$ $newlinech %$$
% $end
