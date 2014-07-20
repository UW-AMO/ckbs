% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2011
% Authors: Bradlely Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
%          Aleksandr Aravkin:    saravkin at eos dot ubc dot ca
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin tridiag_inv_ok.m$$ $newlinech %$$
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
%       inv
%       ifull
%       ind
%       var
% $$
%
% $section ckbs_tridiag_inv Example and Test$$
%
% $index ckbs_tridiag_inv, example and test$$
% $index tridiag_inv, example and test$$
% $index example, tridiag_inv$$
% $index test, tridiag_inv$$
%
% $head Source Code$$
% $newlinech $$ $codep
function [ok] = tridiag_inv_ok()
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
    [var] = ckbs_tridiag_inv(b, c);
    % ----------------------------------------
    %
    full = ckbs_blktridiag_mul(b, c, eye(N*n));
    ifull = inv(full); 
    ifullDiag = zeros(n, n, N);
    
    for k = 1:N
       indStart = (k-1)*n + 1;
       indEnd = k*n;
       ifullDiag(:,:,k) = ifull(indStart:indEnd, indStart:indEnd); 
    end
    
    ok = ok & sum(sum(sum(abs(var - ifullDiag)))) < 1e-10;

    return
end
% $$ $newlinech %$$
% $end
