% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2014
% Authors: Bradley Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
%          Aleksandr Aravkin:    sasha.aravkin at gmail dot com
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin ckbs_tridiag_solve_mf$$ $newlinech %$$
% $spell
%       ckbs
%       blkdiag
%       ckbs_tridiag_solve
%       Tridiagonal
%       Bradley
% $$
%
% $section New 2-parallel Algorithm for Block Tridiagonal Systems$$
%
% $index ckbs_tridiag_solve_new$$
% $index tridiag_solve$$
%
% $index solve, block tridiagonal$$
% $index tridiagonal, block solve$$
%
% $head Syntax$$
% $codei/[/e/, /lambda/] = ckbs_tridiag_solve(/b/, /c/, /r/)/$$
%
% $head Purpose$$
% This routine solves the following linear system of equations for
% $latex e$$:
% $latex \[
%       A * e = r
% \] $$
% where the symmetric block tridiagonal matrix $latex A$$ is defined by
% $latex \[
% A =
% \left( \begin{array}{ccccc}
% b_1    & c_2^\R{T} & 0          & \cdots & 0      \\
% c_2    & b_2       & c_3^\R{T}  &        & \vdots \\
% \vdots &           & \ddots     &        &        \\
% 0      & \cdots    &            & b_N    & c_N
% \end{array} \right)
% \] $$
%
% This routine implements a new parallel 2-filter soltuion, which 
% can be done in 2-parallel, without a post-processing step
% required by the MF 2-filter smoother. 
%
% $head b$$
% The argument $icode b$$ is a three dimensional array,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       b_k = b(:,:,k)
% \]$$
% and $icode b$$ has size $latex n \times n \times N$$.
%
% $head c$$
% The argument $icode c$$ is a three dimensional array,
% for $latex k = 2 , \ldots , N$$
% $latex \[
%       c_k = c(:,:,k)
% \]$$
% and $icode c$$ has size $latex n \times n \times N$$.
%
% $head r$$
% The argument $icode r$$ is an $latex (n * N) \times m$$ matrix.
%
% $head e$$
% The result $icode e$$ is an $latex (n * N) \times m$$ matrix.
%
% $head lambda$$
% The result $icode lambda$$ is a scalar equal to the
% logarithm of the determinant of $latex A$$.
%
% $head Reference$$
% $icode Kalman smoothing and block tridiagonal systems: new connections and numerical stability results$$,
% Aleksandr Y. Aravkin, Bradley M. Bell, James V. Burke, Gianluigi
% Pillonetto,
% axiv 1303.5237, 2013. 
%
%
% $children#
%       example/tridiag_solve_ok.m
% #$$
%
% $head Example$$
% The file $cref tridiag_solve_ok.m$$ contains an example and test of
% $code ckbs_tridiag_solve$$, as well as of the other routines for this problem.
% It returns true if $code ckbs_tridiag_solve$$ passes the test
% and false otherwise.
%
% $end
% ---------------------------------------------------------------------------
function [e, lambda] = ckbs_tridiag_solve_new(b, c, r)
    %
    % Determine the size of the problem
    n = size(b, 1);
    m = size(r, 2);
    N = size(b, 3);
    %
    % Initialize forward and backward structures
    df     = b;
    db     = b;  
    sf     = zeros(n, m, N);
    sb     = zeros(n, m, N);
    blk_n = 1 : n;
    blk_n_b = (n*(N-1) + 1):n*N;
 
        %
    %% Forward and backward half-passes (can be done in parallel) 
    % Step 3 Forward to midpoint
    sk = r(blk_n, :);
    sf(:,:,1) = sk; 
    
    
    midPoint = floor(N/2);
    if(midPoint < 2) 
        error('try a longer time series, buddy')
    end
    
    for k = 2 : 1 : midPoint
        blk_n      = blk_n + n;
        sk1        = sf(:, :, k-1);
        dk1        = df(:, :, k-1);
        sk         = r(blk_n, :);
        ck         = c(:, :, k);
        df(:, :, k) = b(:, :, k) - ck * chol_solve(dk1, ck');
        sf(:, :, k) = sk - ck * chol_solve(dk1, sk1);
    end

    

    % Step 3 backward to midpoint
    sk = r(blk_n_b, :);
    sb(:,:,N) = sk;
    dk = b(:,:,N);
%    db(:,:,N) = b(:,:,N);
    ck = c(:,:,N);
    
     for k = N-1 : -1 : midPoint + 1
        ck1 = ck;
        ck = c(:,:,k);
        dk1 = dk;
        sk1 = sk;
        blk_n_b = blk_n_b - n;
        sk = r(blk_n_b, :);
        bk = b(:,:,k);
        dk = bk - ck1'* chol_solve(dk1, ck1);
        db(:,:,k) = dk;
        sk = sk - ck1'* chol_solve(dk1, sk1);
        sb(:,:,k) = sk;
     end
    
     %% Combining information from forward and backward half passes 
     
     cc = c(:,:,midPoint+1);
     df(:,:,midPoint) = df(:,:,midPoint) - cc'*chol_solve(db(:,:,midPoint+1), cc);
     sf(:,:,midPoint) = sf(:,:,midPoint) - cc'*chol_solve(db(:,:,midPoint+1), sb(:,:,midPoint+1));
     sb(:,:,midPoint+1) = sb(:,:,midPoint+1) - cc*chol_solve(df(:,:,midPoint), sf(:,:,midPoint));
     % seems correct
     lambda = log(det(df(:,:,midPoint)));
     
     
     %% Runnign backward passes to completion 
     
     % forward
     sf(:, :, midPoint) = chol_solve(df(:,:, midPoint),  sf(:,:, midPoint));
     e((midPoint-1)*n+1:midPoint*n,:) = sf(:,:,midPoint);
     for k = midPoint-1 : -1 : 1
        ek1        = sf(:, :, k+1);
        ck1        = c(:, :, k+1);
        dk         = df(:, :, k);
        sf(:, :, k) = chol_solve(dk,  sf(:, :, k) - ck1' * ek1);
        e((k-1)*n+1:k*n,:) = sf(:,:,k);
        lambda = lambda + log(det(dk));
     end
          
    % backward 
    dk = db(:,:,midPoint+1);
    sk = sb(:,:,midPoint+1);
    %'third chol_solve'
    e((midPoint)*n+1:(midPoint+1)*n,:) = chol_solve(dk, sk);

    blk_n = (midPoint)*n+1:(midPoint+1)*n;
    for k = midPoint+2:1:N
        ek1 = e(blk_n, :);
        blk_n = blk_n + n;
        dk = db(:,:,k);
        lambda = lambda + log(det(dk));
        sk = sb(:,:,k);
        ck = c(:,:,k);
        tempk = sk - ck*ek1;
        %    'fourth chol_solve'
        e(blk_n, :) = chol_solve(dk, tempk);
    end

    return
end

function [x] = chol_solve(A, y)
    A_chol = chol(A);
    x      = A_chol \ ( A_chol' \ y );
    return
end
