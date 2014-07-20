% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2011
% Authors: Bradlely Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
%          Aleksandr Aravkin:    saravkin at eos dot ubc dot ca
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin process_hes_ok.m$$ $newlinech %$$
% $spell
%       ckbs
%       dg
%       dh
%       diff
%       end end
%       gradm
%       gradp
%       hes
%       qinv
%       rinv
%       process
%       sumsq
%       tmp
%       xm
%       xp
% $$
%
% $section ckbs_process_hes Example and Test$$
%
% $index ckbs_process_hes, example and test$$
% $index process_hes, example and test$$
% $index example, process_hes$$
% $index test, process_hes$$
%
% $head Source Code$$
% $newlinech $$ $codep
function [ok] = process_hes_ok()
    ok = true;
    % --------------------------------------------------------
    % You can change these parameters
    m    = 1;   % number of measurements per time point
    n    = 2;   % number of state vector components per time point
    N    = 3;   % number of time points
    % ---------------------------------------------------------
    %  Define the problem
    rand('seed', 123);
    dg   = zeros(n, n, N);
    qinv = zeros(n, n, N);
    for k = 1 : N
        dg(:, :, k)   = rand(n, n);
        tmp           = rand(n, n);
        qinv(:, :, k) = (tmp + tmp') / 2 + 2 * eye(n);
    end
    % ---------------------------------------------------------
    % Compute the Hessian using ckbs_process_hes
    [D, A] = ckbs_process_hes(dg, qinv);
    % ---------------------------------------------------------
    H    = zeros(n * N , n * N );
    for k = 1 : N
        index           = (k - 1) * n + (1 : n);
        H(index, index) = D(:, :, k);
        if k > 1
            H(index - n, index) = A(:, :, k)';
            H(index, index - n) = A(:, :, k);
        end
    end
    %
    % Use finite differences to check Hessian
    x    = rand(n, N);
    g    = rand(n, N);
    %
    step   = 1;
    for k = 1 : N
        for i = 1 : n
            % Check second partial w.r.t x(i, k)
            xm       = x;
            xm(i, k) = xm(i, k) - step;
            gradm = ckbs_process_grad(xm, g, dg, qinv);
            %
            xp       = x;
            xp(i, k) = xp(i, k) + step;
            gradp = ckbs_process_grad(xp, g, dg, qinv);
            %
            check  = (gradp - gradm) / ( 2 * step );
            for k1 = 1 : N
                for i1 = 1 : n
                    value = H(i + (k-1)*n, i1 + (k1-1)*n);
                    diff  = check(i1, k1) - value;
                    ok     = ok & ( abs(diff) < 1e-10 );
                end
            end
        end
    end
    return
end
% $$ $newlinech %$$
% $end
