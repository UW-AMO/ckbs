% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2011
% Authors: Bradlely Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
%          Aleksandr Aravkin:    saravkin at eos dot ubc dot ca
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin L2L1_obj_ok.m$$ $newlinech %$$
% $spell
%       ckbs
%       dg
%       dh
%       obj obj
%       qinv
%       rinv
%       sumsq
%       tmp
%       xk
%       xkm
%       xres
%       zres
%       L2L1
%       chol
%       sqrt
% $$
%
% $section ckbs_L2L1_obj Example and Test$$
%
% $index ckbs_L2L1_obj, example and test$$
% $index L2L1_obj, example and test$$
% $index example, L2L1_obj$$
% $index test, L2L1_obj$$
%
% $head Source Code$$
% $newlinech $$ $codep
function [ok] = L2L1_obj_ok()
    ok = true;
    % --------------------------------------------------------
    % You can change these parameters
    m    = 1;   % number of measurements per time point
    n    = 2;   % number of state vector components per time point
    N    = 3;   % number of time points
    % ---------------------------------------------------------
    %  Define the problem
    rand('seed', 123);
    x    = rand(n, N);
    z    = rand(m, N);
    h    = rand(m, N);
    g    = rand(n, N);
    dh   = zeros(m, n, N);
    dg   = zeros(n, n, N);
    qinv = zeros(n, n, N);
    rinv = zeros(m, m, N);
    for k = 1 : N
        dh(:, :, k)   = rand(m, n);
        dg(:, :, k)   = rand(n, n);
        tmp           = rand(m, m);
        rinv(:, :, k) = (tmp + tmp') / 2 + 2 * eye(m);
        tmp           = rand(n, n);
        qinv(:, :, k) = (tmp + tmp') / 2 + 2 * eye(n);
    end
    % ---------------------------------------------------------
    % Compute the Objective using ckbs_L2L1_obj
    obj  = ckbs_L2L1_obj(x, z, g, h, dg, dh, qinv, rinv);
    % ---------------------------------------------------------
    L2L1 = 0;
    xk    = zeros(n, 1);
    for k = 1 : N
        xkm   = xk;
        xk    = x(:, k);
        xres  = xk      - g(:, k) - dg(:,:, k) * xkm;
        zres  = z(:, k) - h(:, k) - dh(:,:, k) * xk;
        L2L1 = L2L1 + 0.5*norm(chol(qinv(:,:,k))*xres)^2;
        L2L1 = L2L1 + sqrt(2)*norm(chol(rinv(:,:,k))*zres, 1);
    end
    ok = ok & ( abs(obj - L2L1) < 1e-10 );
    return
end
% $$ $newlinech %$$
% $end
