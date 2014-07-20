% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2013
% Authors: Aleksandr Y. Aravkin: saravkin at us dot ibm dot com
%          Bradley Bell:         bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
% License: GNU General Public License Version 2
% -------------------------------------------------------------------

% $begin t_grad_ok.m$$ $newlinech %$$
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
%       t
%       tmp
%       xm
%       xp
%       params
%       pos
%       vel
%       dt
%       df
%       inds
%       eps
%       randn
% $$
%
% $section ckbs_t_grad Example and Test$$
%
% $index ckbs_t_grad, example and test$$
% $index t_grad, example and test$$
% $index example, t_grad$$
% $index test, t_grad$$
%
% $head Source Code$$
% $newlinech $$ $codep
function [ok] = t_grad_ok()
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
params.direct_h_index = 2;
h_fun = @(k,x) direct_h(k,x,params);

params.pos_vel_g_dt = .01;
params.pos_vel_g_initial = zeros(n,1);
g_fun = @(k,x) pos_vel_g(k,x,params);

df_meas = 4;
df_proc = 4;
params.df_proc = df_proc;
params.df_meas = df_meas;
eps = 1e-5;

% ---------------------------------------------------------
% Compute the gradient using ckbs_t_grad


params.inds_proc_st = [];
params.inds_meas_st = [];
%%%%%%%
g1 = ckbs_t_grad(x, z, g_fun, h_fun, qinv, rinv, params);
f1 = ckbs_t_obj(x, z, g_fun, h_fun, qinv, rinv, params);

diff = eps*randn(size(x));
x2 = x + diff;
g2 = ckbs_t_grad(x2, z, g_fun, h_fun, qinv, rinv, params);
f2 = ckbs_t_obj(x2, z, g_fun, h_fun, qinv, rinv, params);
temp = 0.5 *(g1(:) + g2(:))'*(diff(:))/(f2 - f1);
y = abs(1-temp);

ok = ok & (abs(y) < 1e-8);
% ---------------------------------------------------------

params.inds_proc_st = 1;
params.inds_meas_st = [];
%%%%%%%
g1 = ckbs_t_grad(x, z, g_fun, h_fun, qinv, rinv, params);
f1 = ckbs_t_obj(x, z, g_fun, h_fun, qinv, rinv, params);

diff = eps*randn(size(x));
x2 = x + diff;
g2 = ckbs_t_grad(x2, z, g_fun, h_fun, qinv, rinv, params);
f2 = ckbs_t_obj(x2, z, g_fun, h_fun, qinv, rinv, params);
temp = 0.5 *(g1(:) + g2(:))'*(diff(:))/(f2 - f1);
y = abs(1-temp);

ok = ok & (abs(y) < 1e-8);

% ---------------------------------------------------------

params.inds_proc_st = 1:2;
params.inds_meas_st = [];
%%%%%%%
g1 = ckbs_t_grad(x, z, g_fun, h_fun, qinv, rinv, params);
f1 = ckbs_t_obj(x, z, g_fun, h_fun, qinv, rinv, params);

diff = eps*randn(size(x));
x2 = x + diff;
g2 = ckbs_t_grad(x2, z, g_fun, h_fun, qinv, rinv, params);
f2 = ckbs_t_obj(x2, z, g_fun, h_fun, qinv, rinv, params);
temp = 0.5 *(g1(:) + g2(:))'*(diff(:))/(f2 - f1);
y = abs(1-temp);

ok = ok & (abs(y) < 1e-8);

% ---------------------------------------------------------

params.inds_proc_st = 2;
params.inds_meas_st = [];
%%%%%%%
g1 = ckbs_t_grad(x, z, g_fun, h_fun, qinv, rinv, params);
f1 = ckbs_t_obj(x, z, g_fun, h_fun, qinv, rinv, params);

diff = eps*randn(size(x));
x2 = x + diff;
g2 = ckbs_t_grad(x2, z, g_fun, h_fun, qinv, rinv, params);
f2 = ckbs_t_obj(x2, z, g_fun, h_fun, qinv, rinv, params);
temp = 0.5 *(g1(:) + g2(:))'*(diff(:))/(f2 - f1);
y = abs(1-temp);

ok = ok & (abs(y) < 1e-8);
% ---------------------------------------------------------

params.inds_proc_st = [];
params.inds_meas_st = 1;
%%%%%%%
g1 = ckbs_t_grad(x, z, g_fun, h_fun, qinv, rinv, params);
f1 = ckbs_t_obj(x, z, g_fun, h_fun, qinv, rinv, params);

diff = eps*randn(size(x));
x2 = x + diff;
g2 = ckbs_t_grad(x2, z, g_fun, h_fun, qinv, rinv, params);
f2 = ckbs_t_obj(x2, z, g_fun, h_fun, qinv, rinv, params);
temp = 0.5 *(g1(:) + g2(:))'*(diff(:))/(f2 - f1);
y = abs(1-temp);

ok = ok & (abs(y) < 1e-8);
% ---------------------------------------------------------

params.inds_proc_st = 1;
params.inds_meas_st = 1;
%%%%%%%
g1 = ckbs_t_grad(x, z, g_fun, h_fun, qinv, rinv, params);
f1 = ckbs_t_obj(x, z, g_fun, h_fun, qinv, rinv, params);

eps = 1e-5;
diff = eps*randn(size(x));
x2 = x + diff;
g2 = ckbs_t_grad(x2, z, g_fun, h_fun, qinv, rinv, params);
f2 = ckbs_t_obj(x2, z, g_fun, h_fun, qinv, rinv, params);
temp = 0.5 *(g1(:) + g2(:))'*(diff(:))/(f2 - f1);
y = abs(1-temp);

ok = ok & (abs(y) < 1e-8);

% ---------------------------------------------------------

params.inds_proc_st = 2;
params.inds_meas_st = 1;
%%%%%%%
g1 = ckbs_t_grad(x, z, g_fun, h_fun, qinv, rinv, params);
f1 = ckbs_t_obj(x, z, g_fun, h_fun, qinv, rinv, params);

diff = eps*randn(size(x));
x2 = x + diff;
g2 = ckbs_t_grad(x2, z, g_fun, h_fun, qinv, rinv, params);
f2 = ckbs_t_obj(x2, z, g_fun, h_fun, qinv, rinv, params);
temp = 0.5 *(g1(:) + g2(:))'*(diff(:))/(f2 - f1);
y = abs(1-temp);

ok = ok & (abs(y) < 1e-8);

% ---------------------------------------------------------

params.inds_proc_st = 1:2;
params.inds_meas_st = 1;
%%%%%%%
g1 = ckbs_t_grad(x, z, g_fun, h_fun, qinv, rinv, params);
f1 = ckbs_t_obj(x, z, g_fun, h_fun, qinv, rinv, params);

diff = eps*randn(size(x));
x2 = x + diff;
g2 = ckbs_t_grad(x2, z, g_fun, h_fun, qinv, rinv, params);
f2 = ckbs_t_obj(x2, z, g_fun, h_fun, qinv, rinv, params);
temp = 0.5 *(g1(:) + g2(:))'*(diff(:))/(f2 - f1);
y = abs(1-temp);

ok = ok & (abs(y) < 1e-8);

return
end
% $$ $newlinech %$$
% $end
