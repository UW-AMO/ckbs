% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2013
% Authors: Aleksandr Y. Aravkin: saravkin at us dot ibm dot com
%          Bradley Bell:         bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin t_hess_ok.m$$ $newlinech %$$
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
%       t
%       tmp
%       xm
%       xp
%       hess
%       params
%       pos
%       vel
%       df
%       inds
%       blktridiag
%       mul
%       dt
%       
% $$
%
% $section ckbs_t_hess Example and Test$$
%
% $index ckbs_t_hess, example and test$$
% $index t_hess, example and test$$
% $index example, t_hess$$
% $index test, t_hess$$
%
% $head Source Code$$
% $newlinech $$ $codep
function [ok] = t_hess_ok()
    ok = true;
% --------------------------------------------------------
% You can change these parameters
m    = 1;   % number of measurements per time point
n    = 2;   % number of state vector components per time point
N    = 1;   % number of time points
% ---------------------------------------------------------
%  Define the problem
rand('seed', 123);
%  Define the problem
rand('seed', 123);
x    = rand(n, N);
z    = rand(m, N);
qinv = zeros(n, n, N);
rinv = zeros(m, m, N);
for k = 1 : N
    tmp           = rand(m, m);
    rinv(:, :, k) = (tmp + tmp') / 2 + 2 * eye(m);
    tmp           = rand(n, n);
    qinv(:, :, k) = (tmp + tmp') / 2 + 2 * eye(n);
end
params.direct_h_index = 1;
h_fun = @(k,x) direct_h(k,x,params);

params.pos_vel_g_dt = .05;
params.pos_vel_g_initial = zeros(n,1);
g_fun = @(k,x) pos_vel_g(k,x,params);

df_meas = 4;
df_proc = 4;
params.df_proc = df_proc;
params.df_meas = df_meas;
params.inds_proc_st = 2;
params.inds_meas_st = 1;
    % ---------------------------------------------------------
    % Compute the Hessian using ckbs_t_hess
    [D, A] = ckbs_t_hess(x, z, g_fun, h_fun, qinv, rinv, params);
    % ---------------------------------------------------------
    H = ckbs_blktridiag_mul(D, A, eye(n*N));
    %
    % Use finite differences to check Hessian
    x    = rand(n, N);
    z    = rand(m, N);
    %
    step   = 1e-5;
    for k = 1 : N
        for i = 1 : n
            % Check second partial w.r.t x(i, k)
            xm       = x;
            xm(i, k) = xm(i, k) - step;
            gradm = ckbs_t_grad(xm, z, g_fun, h_fun, qinv, rinv, params);
            %
            xp       = x;
            xp(i, k) = xp(i, k) + step;
            gradp = ckbs_t_grad(xp, z, g_fun, h_fun, qinv, rinv, params);
            %
            check  = (gradp - gradm) / ( 2 * step );
            for k1 = 1 : N
                for i1 = 1 : n
                    value = H(i + (k-1)*n, i1 + (k1-1)*n);
                    diff  = check(i1, k1) - value
                    ok     = ok & ( abs(diff) < 1e-10 );
                end
            end
        end
    end
    return
end
% $$ $newlinech %$$
% $end
