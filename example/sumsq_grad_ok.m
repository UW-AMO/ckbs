% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006
% Authors: Bradlely Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin sumsq_grad_ok.m$$ $newlinech %$$
% $spell
%       ckbs
%       dg
%       dh
%       diff
%       end end
%       grad grad
%       obj
%       qinv
%       rinv
%       sm
%       sumsq
%       tmp
%       xm
%       xp
% $$
%
% $section ckbs_sumsq_grad Example and Test$$
%
% $index ckbs_sumsq_grad, example and test$$
% $index sumsq_grad, example and test$$
% $index example, sumsq_grad$$
% $index test, sumsq_grad$$
%
% $head Source Code$$
% $newlinech $$ $codep
function [ok] = sumsq_grad_ok()
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
    dg   = zeros(n, n, N);
    dh   = zeros(m, n, N);
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
    % Compute the gradient using ckbs_sumsq_grad
    grad = ckbs_sumsq_grad(x, z, g, h, dg, dh, qinv, rinv);
    % ---------------------------------------------------------
    % Use finite differences to check gradient
    step   = 1;
    for k  = 1 : N
        for i  = 1 : n
            % Check second partial w.r.t x(i,k)
            xm        = x;
            xm(i, k)  = xm(i, k) - step;
            Sm        = ckbs_sumsq_obj(xm, z, g, h, dg, dh, qinv, rinv);
            %
            xp        = x;
            xp(i, k)  = xp(i, k) + step;
            Sp        = ckbs_sumsq_obj(xp, z, g, h, dg, dh, qinv, rinv);
            %
            check  = (Sp - Sm) / ( 2 * step);
            diff   = grad(i, k) - check;
            ok     = ok & ( abs(diff) < 1e-10 );
        end
    end
    return
end
% $$ $newlinech %$$
% $end
