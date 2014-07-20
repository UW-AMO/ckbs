% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2013
% Authors: Aleksandr Y. Aravkin: saravkin at us dot ibm dot com
%          Bradley Bell:         bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin ckbs_affine_singular$$ $newlinech %$$
% $spell
%       itr
%       obj
%       ckbs
%       dg
%       dh
%       qinv
%       rinv
%       optimality
%       pairwise
%       semidefinite
% $$
%
% $index ckbs_affine_singular$$
%
% $index affine, singular smoother$$
% $index singular, affine smoother$$
% $index smoother, affine singular$$
%
% $section Singular Affine Kalman Bucy Smoother$$
%
% $head Syntax$$
% $codei/[/xOut/,  /info/] = ckbs_affine_singular(...
%         /z/, /g/, /h/, ...
%        /dg/, /dh/, /q/, /r/)/$$
%
% $head Purpose$$
% This routine finds the smoothed estimate for affine process and measurement 
% models when the variance matrices $latex Q$$ and $latex R$$ may be singular.
% 
%
% $head Notation$$
% The singular Kalman-Bucy smoother state is given by 
% $latex \[
% \begin{array}{rcl}
% G^{-1} (w - Q G^{-T}H^T\Phi^{-1}(HG^{-1}w-v))
% \end{array}
% \] $$
% where 
% $latex \[
% \Phi = HG^{-1}QG^{-T}H^T + R ,
% \] $$
% and the matrices $latex R_k$$ and $latex Q_k$$ are
% symmetric positive semidefinite. 
% Note that $latex g_1$$ is the initial state estimate
% and $latex Q_1$$ is the corresponding covariance.
%
%
% $head z$$
% The argument $icode z$$ is a two dimensional array,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       z_k = z(:, k)
% \]$$
% and $icode z$$ has size $latex m \times N$$.
%
% $head g$$
% The argument $icode g$$ is a two dimensional array,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       g_k = g(:, k)
% \]$$
% and $icode g$$ has size $latex n \times N$$.
% The value $latex g_1$$ serves as the initial state estimate.
%
% $head h$$
% The argument $icode h$$ is a two dimensional array,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       h_k = h(:, k)
% \]$$
% and $icode h$$ has size $latex m \times N$$.
%
%
% $head dg$$
% The argument $icode dg$$ is a three dimensional array,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       G_k = dg(:,:,k)
% \]$$
% and $icode dg$$ has size $latex n \times n \times N$$.
% The initial state estimate $latex g_1$$ does not depend on the value of
% $latex x_0$$, hence $latex G_1$$ should be zero.
%
% $head dh$$
% The argument $icode dh$$ is a three dimensional array,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       H_k = dh(:,:,k)
% \]$$
% and $icode dh$$ has size $latex m \times n \times N$$.
%
% $head q$$
% The argument $icode q$$ is a three dimensional array,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       Q_k = q(:,:, k)
% \]$$
% and $icode q$$ has size $latex n \times n \times N$$.
% The value of $latex Q_k$$ is the variance of the initial state
% estimate $latex g_1$$.
%
% $head r$$
% The argument $icode r$$ is a three dimensional array,
% for $latex k = 1 , \ldots , N$$
% $latex \[
%       R_k = r(:,:, k)
% \]$$
% and $icode r$$ has size $latex m \times m \times N$$.
% It is ok to signify a missing data value by setting the corresponding
% row and column of $icode r$$ to infinity. This also enables use
% to interpolate the state vector $latex x_k$$ to points where
% there are no measurements.
%
% $head xOut$$
% The result $icode xOut$$ contains the optimal sequence
% $latex ( x_1 , \ldots , x_N )$$.
% For $latex k = 1 , \ldots , N$$
% $latex \[
%       x_k = xOut(:, k)
% \]$$
% and $icode xOut$$ is a two dimensional array with size $latex n \times N$$.
%
%
% $head info$$
% Contains infinity norm of each of the three KKT equations for the
% constrained reformulation. 
% 
% $pre
%
% $$
%
% $children#
%       example/affine_singular_ok.m
% #$$
%
% $head Example$$
% The file $cref affine_singular_ok.m$$ contains an example and test of
% $code ckbs_affine_singular$$.
% It returns true if $code ckbs_affine_singular$$ passes the test
% and false otherwise.
%
% $end
% ----------------------------------------------------------------------------


function [xOut, info] = ckbs_affine_singular(z, g, h, dg, dh, q, r)

    info = [];
    if nargin ~= 7
        error('ckbs_affine: improper number of input arguments');
    end
    if nargout ~= 2
        error('ckbs_affine: improper number of return values');
    end
    % size of problem
    n     = size(g, 1);
    N     = size(z, 2);
    m     = size(z,   1);
    %ell   = size(b,   1);
    %
    % check sizes
    if  N ~= size(h,2) || N ~= size(dg,3) || N ~= size(dh,3) || N ~= size(q,3) || N ~= size(r,3)
        N
        size(z,2)
        size(h,2)
        size(dg,3)
        size(dh,3)
        size(q,3)
        size(r,3)
        error('ckbs_affine: argument sizes do not agree');
    end
    if n ~= size(dg,2) || n ~= size(dh,2) || ...
            n ~= size(q,1) ||n ~= size(q,2)
        n
        size(g,1)
        size(dg,2)
        size(dh,2)
        size(q,1)
        size(q,2)
        error('ckbs_affine: argument sizes do not agree');
    end
    if m ~= size(h,1) | m ~= size(dh,1) | m ~= size(r,1) | m ~= size(r,2)
        m
        size(h,1)
        size(dh,1)
        size(r,1)
        size(r,2)
        error('ckbs_affine: argument sizes do not agree');
    end
    %
    % Other usefull sizes
    %meas        = ell * N;
  
  
    
    
    packedIdn = zeros(n, n, N);
    packedIdm = zeros(m, m, N);
    pinvQ = zeros(n, n, N);
    pinvR = zeros(m, m, N);
    hrh = zeros(n, n, N);
    for k = 1:N
       packedIdn(:,:,k) = eye(n);
       packedIdm(:,:,k) = eye(m);
     %  pinvQ(:,:,k) = pinv(q(:,:,k));
     %  pinvR(:,:,k) = pinv(r(:,:,k));
      % hrh(:,:,k) = dh(:,:,k)'*pinv(r(:,:,k))*dh(:,:,k);
    end
  %% Preconditioning using Woodbury  
  %  [gghh, ggo] = ckbs_blkbidiag_symm_mul(packedId, -dg, pinvQ);
  %  gghh = gghh + hrh; 
  %  precon = @(x)preconAction(x, dh, gghh, ggo, pinvR);
    
    w = g(:);
    v = z(:);
    
    GinvW =  ckbs_bidiag_solve(packedIdn, -dg, w);
    HGinvWmV =  ckbs_blkdiag_mul(dh, GinvW) - v;
    
    %% Preconditioning using simpler block tridiagonal system. 
    hg = zeros(m, n, N);
    for k=2:N
       hg(:,:,k) = dh(:,:,k)*dg(:,:,k);
    end    
    %[pd, po] = ckbs_blkbidiag_symm_mul(dh, hg, q); % This is HVQV'H', with V 'like' G^{-1}
    %pd = pd + 10*packedIdm; % put something big in instead of R
    %precon = @(x)ckbs_diag_solve(pd, x);
    
    %preconAct = @(x)specialAction(x, dh, dg, q, 10*packedIdm, packedIdn);
    %precon = @(x)pcg(preconAct, x, 1e-2);
    %precon = @(x)ckbs_diag_solve(r, x);

    %%
    [lambdaTwo, flag, relres, iter] = pcg(@(x)specialAction(x, dh, dg, q, r, packedIdn), HGinvWmV, [], n*N);
    lambdaOne = ckbs_blkdiag_mul_t(dh, lambdaTwo);
    lambdaOne = -ckbs_bidiag_solve_t(packedIdn, -dg, lambdaOne);
    
    QlamOne = ckbs_blkdiag_mul(q, lambdaOne);
    xOut = ckbs_bidiag_solve(packedIdn, -dg, w + QlamOne);
    xOut = reshape(xOut, n, N);
    
%     zv = z(:);
%     
%     iHH = zeros(m, m, N);
%     HiRiH = zeros (n, n, N);
%     
%     for k = 1:N
%        packedId(:,:,k) = eye(n);
%        iHH(:,:,k) = inv(dh(:,:,k) * dh(:,:,k)'); % the HH^T matrix 
%        HiRiH(:,:,k) = dh(:,:,k)' * iHH(:,:,k) * r(:,:,k) * iHH(:,:,k) * dh(:,:,k);
%     end
%     
%     iHHz = ckbs_blkdiag_mul(iHH, zv);
%     HiHHz = ckbs_blkdiag_mul_t(dh, iHHz);
%     GHiHHz = ckbs_blkbidiag_mul(packedId, -dg, HiHHz);
%     
%     w = g(:);
% 
%     rhs = GHiHHz - w;
%     
%     [PhiDia, PhiOffDia] = ckbs_blkbidiag_symm_mul(packedId, -dg, HiRiH);
%     PhiDia = PhiDia + q;
%    
%     lambdaOne = -ckbs_tridiag_solve_b(PhiDia, PhiOffDia, rhs);
%    
%     QlambdaOne = ckbs_blkdiag_mul(q, lambdaOne);
%     
%     xOut = ckbs_bidiag_solve(packedId, -dg, w - QlambdaOne);
%     
%     xOut = reshape(xOut, n, N);
%     
    info = [0 0]; % temporary holder
  
    % KKT check
    GtlambdaOne = ckbs_blkbidiag_mul_t(packedIdn, -dg, lambdaOne);
    HGtlambdaOne = ckbs_blkdiag_mul(dh, GtlambdaOne);
    %lambdaTwo = -ckbs_blkdiag_mul(iHH, HGtlambdaOne);
    
    K0 = specialAction(lambdaTwo, dh, dg, q, r, packedIdn) - HGinvWmV;
    K1 = ckbs_blkbidiag_mul(packedIdn, -dg, xOut(:)) - ckbs_blkdiag_mul(q, lambdaOne) - w;
    K2 = ckbs_blkdiag_mul(dh, xOut(:)) - ckbs_blkdiag_mul(r, lambdaTwo) - v;
    K3 = GtlambdaOne + ckbs_blkdiag_mul_t(dh, lambdaTwo);
   
    info = [norm(K0, inf); norm(K1, inf); norm(K2, inf); norm(K3, inf); iter];
    
return
end

function [v] = specialAction(x, dh, dg, q, r, packedId)

n = size(dg, 1);
N = size(dg, 3);

v = ckbs_blkdiag_mul(r, x);
temp = ckbs_blkdiag_mul_t(dh, x);
temp = ckbs_bidiag_solve_t(packedId, -dg, temp);
temp = ckbs_blkdiag_mul(q, temp);
temp = ckbs_bidiag_solve(packedId, -dg, temp);
temp = ckbs_blkdiag_mul(dh, temp);

v = v + temp;
end

function [v] = preconAction(x, dh, gghh, ggo, pinvr)

temp = ckbs_blkdiag_mul(pinvr, x);
temp = ckbs_blkdiag_mul_t(dh, temp);
temp = ckbs_diag_solve_b(gghh, ggo, temp);
temp = ckbs_blkdiag_mul(dh, temp);
v = ckbs_blkdiag_mul(pinvr, x - temp);

end
