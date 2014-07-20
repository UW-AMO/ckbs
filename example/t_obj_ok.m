% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2013
% Authors: Aleksandr Y. Aravkin: saravkin at us dot ibm dot com
%          Bradley Bell:         bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin t_obj_ok.m$$ $newlinech %$$
% $spell
%       ckbs
%       dg
%       dh
%       obj obj
%       qinv
%       rinv
%       t
%       tmp
%       xk
%       xkm
%       xres
%       zres
%       params
%       pos
%       vel
%       df
%       inds
%       res
%       gk
%       feval
%       dt
% $$
%
% $section ckbs_t_obj Example and Test$$
%
% $index ckbs_t_obj, example and test$$
% $index t_obj, example and test$$
% $index example, t_obj$$
% $index test, t_obj$$
%
% $head Source Code$$
% $newlinech $$ $codep
function [ok] = t_obj_ok()
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
params.inds_proc_st = [];
params.inds_meas_st = 1;

% ---------------------------------------------------------
% Compute the Objective using ckbs_t_obj
obj  = ckbs_t_obj(x, z, g_fun, h_fun, qinv, rinv, params);
% ---------------------------------------------------------
res_z = zeros(m, N);
for k = 1 : N
    xk    = x(:, k);
    res_z(:,k)    = z(:,k) - feval(h_fun, k, xk);
end

res_x = zeros(n, N);
xk = zeros(n, 1);
for k = 1 : N
	xkm   = xk;
	xk    = x(:, k);
	gk    = feval(g_fun, k, xkm);
	res_x(:,k)  = xk - gk;
end

obj2 = 0;
for k = 1 : N
    qinv_k = qinv(:, :, k);
    rinv_k = rinv(:, :, k);

    res_x_k = res_x(:,k);
    res_z_k = res_z(:,k);
    
    obj2   = obj2 + 0.5*(df_meas)*log(1 + (1/df_meas)*res_z_k' * rinv_k * res_z_k) + ...
        0.5*res_x_k' * qinv_k * res_x_k;

end
ok = ok & ( abs(obj - obj2) < 1e-10 );
return
end
% $$ $newlinech %$$
% $end
