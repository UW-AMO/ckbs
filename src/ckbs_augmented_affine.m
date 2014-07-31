% augmented affine solver: takes in matrices B and D.
% xOut: estimates of time-dependent variables
% yOut: estimates of time-independent variables

% only new input is P, the measurement matrix for non-time varying
% variables. P = 1 \kron Pone.
% written by Aleksandr Aravkin and Karthikeyan Natesan Ramamurthy
% ----------------------------------------------------------------------------
function [xOut, yOut, yVar] = ckbs_augmented_affine(max_itr, epsilon, ...
    z, b, g, h, db, dg, dh, qinv, rinv, Pone)
 %   if nargin ~= 12
 %       error('ckbs_affine: improper number of input arguments');
 %   end
%    if nargout ~= 2
%        error('ckbs_affine: improper number of return values');
%    end
    % size of problem
    n     = size(g, 1);
    N     = size(z, 2);
    m     = size(z,   1);
    ell   = size(b,   1);
    %
    % check sizes
    if N ~= size(b,2) | N ~= size(h,2) | N ~= size(db,3) | ...
            N ~= size(dg,3) | N ~= size(dh,3) | N ~= size(qinv,3) | N ~= size(rinv,3)
        N
        size(z,2)
        size(b,2)
        size(h,2)
        size(db,3)
        size(dg,3)
        size(dh,3)
        size(qinv,3)
        size(rinv,3)
        error('ckbs_affine: argument sizes do not agree');
    end
    if n ~= size(db,2) | n ~= size(dg,2) | n ~= size(dh,2) | ...
            n ~= size(qinv,1) | n ~= size(qinv,2)
        n
        size(g,1)
        size(db,2)
        size(dg,2)
        size(dh,2)
        size(qinv,1)
        size(qinv,2)
        error('ckbs_affine: argument sizes do not agree');
    end
    if m ~= size(h,1) | m ~= size(dh,1) | m ~= size(rinv,1) | m ~= size(rinv,2)
        m
        size(h,1)
        size(dh,1)
        size(rinv,1)
        size(rinv,2)
        error('ckbs_affine: argument sizes do not agree');
    end
    if ell ~= size(db,1)
        ell
        size(db,1)
        error('ckbs_affine: argument sizes do not agree');
    end
    %
    % Other usefull sizes
    r        = ell * N;
    p        = n   * N;
    %
    % vector of zero state values
    xZero = zeros(n, N);
    %
    % diagonal and off diagonal blocks for Hessian of the objective
    [D, A] = ckbs_sumsq_hes(dg, dh, qinv, rinv);
    %
    % gradient of the objective that corresponds to xIn
    d      = ckbs_sumsq_grad(xZero, z, g, h, dg, dh, qinv, rinv);
  
    
    % gradient of the objective that corresponds to y_in
    invRz = ckbs_blkdiag_mul(rinv, z(:));
    grady = sum(Pone'*invRz, 1); 
    
%    h1tr = dh;
    B = zeros(n*N, size(Pone,1));
    h2rh2t = zeros(size(Pone,2));
    for ind = 1:N
%        h1tr(:,:,ind) = h1tr(:,:,ind)'*rinv(:,:,ind);
        inds = n*(ind-1)+1:n*ind;
        B(inds, :) = dh(:,:,ind)'*rinv(:,:,ind)*Pone;
        h2rh2t = h2rh2t + Pone'*rinv(:,:,ind)*Pone;
    end
%    B = ckbs_blkdiag_mul_t(h1tr, repmat(Pone, N, 1));
    %B = h1trh2';
    DD = h2rh2t;
    
    [e1 e2 var] = ckbs_augmented_solve(D, A, B, DD, -d(:), grady(:));
    
    xOut        = reshape(e1, n, N);
    yOut        = e2; 
    yVar        = var;

    return
end
