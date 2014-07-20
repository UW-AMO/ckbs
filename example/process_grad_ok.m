% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2011
% Authors: Bradlely Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
%          Aleksandr Aravkin:    saravkin at eos dot ubc dot ca
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin process_grad_ok.m$$ $newlinech %$$
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
%       process
%       tmp
%       xm
%       xp
%       sumsq
% $$
%
% $section ckbs_process_grad Example and Test$$
%
% $index ckbs_process_grad, example and test$$
% $index process_grad, example and test$$
% $index example, process_grad$$
% $index test, process_grad$$
%
% $head Source Code$$
% $newlinech $$ $codep
function [ok] = process_grad_ok()
    ok = true;
    % --------------------------------------------------------
    % You can change these parameters
    m    = 1;   % number of measurements per time point
    n    = 2;   % number of state vector components per time point
    N    = 3;   % number of time points
    % ---------------------------------------------------------
    %  Define the problem
    rand('seed', 123);

    % Measurement terms are fixed at zero
    z    = zeros(m, N);
    h    = zeros(m, N);
    dh   = zeros(m, n, N);
    rinv = zeros(m, m, N);


    % Initialization of process terms
    x    = rand(n, N);
    g    = rand(n, N);
    dg   = zeros(n, n, N);
    qinv = zeros(n, n, N);
    for k = 1 : N
        dg(:, :, k)   = rand(n, n);
        tmp           = rand(n, n);
        qinv(:, :, k) = (tmp + tmp') / 2 + 2 * eye(n);
    end
    % ---------------------------------------------------------
    % Compute the gradient using ckbs_process_grad
    grad = ckbs_process_grad(x, g,dg, qinv);
    % ---------------------------------------------------------
    % Use finite differences to check gradient
    step   = 1;
    for k  = 1 : N
        for i  = 1 : n
            % Check second partial w.r.t x(i,k), setting measurement
            % piece to 0
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
