% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2013
% Authors: Aleksandr Aravkin:    saravkin at us dot ibm dot com
%          Bradley Bell:         bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin vanderpol_experiment_simple.m$$ $newlinech %$$
% $spell
%       ko
%       hk
%       Hk
%       sumsq
%       clf
%       src
%       dg
%       dh
%       xkm
%       uk
%       gk
%       Gk
%       linspace
%       randn
%       rinv
%       qinv
%       nonlinearity
%       Runge-Kutta
%       sim
%       ckbs
%       mu
%       est
%       itr
%       xi
%       optimizer
%       dt
%       xk
%       tk
%       dxk
%       nargin
%       vanderpol
%       Gillijns
%       Kandepu Foss Imsland
%       params
%       end end
%       Var
%       Freq
%       mse
%       qa
%       sqrt
%       strcmp
%       vdp
%       nnz
%       exprnd
%       binornd
%       nums
%       bi
%       trnd
%       componentwise
%       addpath
%       Linewidth
%       bo
%       hleg
%       xlabel
%       ylabel
%       unidrnd
%       gca
%       savefig
%       pwd
% $$
%
% $section ckbs_nonlinear:
% Unconstrained Nonlinear Transition Model Example and Test$$
%
% $index ckbs_nonlinear, unconstrained$$
% $index example, unconstrained$$
% $index test, unconstrained$$
%
% $index vanderpol, unconstrained$$
% $index unconstrained, ckbs_nonlinear$$
%
% $index example, unconstrained$$
% $index test, unconstrained$$
%
% $index example, transition nonlinear$$
% $index test, transition nonlinear$$
%
% $head Syntax$$
% $codei%[%ok%] = vanderpol_ok(%draw_plot%)%$$
%
% $head draw_plot$$
% If this optional argument is true, a plot is drawn showing the results.
% If $icode draw_plot$$ is not present or false, no plot is drawn.
%
% $head ok$$
% If the return value $icode ok$$ is true, the test passed,
% otherwise the test failed.
%
% $head Reference$$
% $index reference$$
% This example appears in many places; for example
%
% $list number$$
% Section 4.1.1 of
% $italic
% Applying the unscented Kalman filter for nonlinear state estimation
% $$, R. Kandepu, B. Foss , L. Imsland,
% Journal of Process Control,
% vol. 18, pp. 753, 2008.
%
% $lnext
% Section V of
% $italic
% What Is the Ensemble Kalman Filter and How Well Does it Work?
% $$,
% S. Gillijns et-al.
% Proceedings of the 2006 American Control Conference,
% Minneapolis, June 14-16, 2006
%
% $lend
%
% $head Simulated State Vector$$
% We use $latex x1 (t)$$ and $latex x2 (t)$$
% to denote the oscillator position and velocity as a function of time.
% The deterministic ordinary differential equation for the
% Van der Pol oscillator is
% $latex \[
% \begin{array}{rcl}
%       x1 '(t) & = & x2 (t)
%       \\
%       x2 '(t) & = & \mu [ 1 - x1(t)^2 ] x2 (t) - x1(t)
% \end{array}
% \] $$
% The simulated state vector values are the solution to this ODE
% with nonlinear parameter $latex \mu = 2$$ and initial conditions
% $latex x1 (0) = 0$$, $latex x2 (0) = -5.$$.
% The routine $cref vanderpol_sim$$ calculates these values.
%
% $head Transition Model$$
% We simulate imperfect knowledge in the dynamical system
% by using Euler's approximation for the ODE with a step size of
% $latex \Delta t = .1$$. In other words,
% given $latex x_{k-1}$$, the value of the state at time $latex t_{k-1}$$,
% we model $latex x_k$$, the state value at time
% $latex t_k = t_{k-1} + \Delta t$$,
% by
% $latex \[
% \begin{array}{rcl}
%       x_{1,k} & = & x_{1,k-1} + x_{2,k-1} \Delta t + w_{1,k}
%       \\
%       x_{2,k} & = & x_{2,k-1} +
%       [ \mu ( 1 - x_{1,k-1}^2 ) x_{2,k-1} - x_{1,k-1} ] \Delta t
%       + w_{2,k}
% \end{array}
% \] $$
% where $latex w_k \sim \B{N}( 0 , Q_k )$$
% and the variance $latex Q_k = \R{diag} ( 0.01, 0.01 )$$.
% The routine $cref vanderpol_g.m$$ calculates the mean for $latex x_k$$
% given $latex x_{k-1}$$.
% Note that the initial state estimate is calculated by this routine as
% $latex \hat{x} = g_1 ( x_0 )$$ and has $latex \hat{x}_1 = \hat{x}_2 = 0$$
% and $latex Q_1 = \R{diag} ( 100. , 100. )$$.
%
% $head Measurements$$
% For this example, the measurements are noisy direct observations of position
% $latex \[
% \begin{array}{rcl}
%       z_{1,k} & = & x_{1,k} + v_{1,k}
% \end{array}
% \] $$
% where $latex v_k \sim \B{N}( 0 , R_k )$$ and
% the measurement variance $latex R_k = \R{diag} ( 0.01, 0.01 )$$.
% The routine $cref direct_h.m$$ calculates the mean for $latex z_k$$,
% given the value of $latex x_k$$; i.e. $latex x_{1,k}$$.
%
% $head Constraints$$
% There are no constraints for this example so we define
% $latex \[
%       f_k ( x_k ) = - 1
% \] $$
% The routine $cref no_f.m$$ is used for these calculations.
%
% $head Call Back Functions$$
% $table
% $rref vanderpol_sim$$
% $rref vanderpol_g.m$$
% $rref direct_h.m$$
% $rref no_f.m$$
% $tend
%
% $head Source Code$$
% $newlinech $$ $codep

function [err]= vanderpol_experiment_simple(draw_plot)

close all;

conVar = 0;
conFreq = 0;
noise = 'normal';

mseKS = zeros(1, 1);
mseRKS = zeros(1, 1);

% if nargin < 1
%     draw_plot = false;
% end
% ---------------------------------------------------------------
% Simulation problem parameters
max_itr = 100;           % Maximum number of optimizer iterations
epsilon = 1e-3;         % Convergence criteria
N       = 164;           % Number of time points (N > 1)
T       = 16.;           % Total time
ell     = 1;            % Number of constraints (never active)
n       = 2;            % Number of components in the state vector (n = 2)
m       = 1;            % Number of measurements at each time point (m <= n)
xi      = [ 0. , -5.]'; % True initial value for the state at time zero
x_in    = zeros(n, N);  % Initial sequence where optimization begins
% SASHA: changed 0, .5 to 0 0
x_est   = [ .1 , -.4]'; % Estimate of initial state (at time index one)
mu      = 2.0;          % ODE nonlinearity parameter
sigma_r = 1.;           % Standard deviation of measurement noise
sigma_q = .1;  % Standard deviation of the transition noise
sigma_qa = .1; % 'Actual' noise


params.vanderpol_g_initial  = xi;
params.vanderpol_g_dt       = T/N;
params.vanderpol_g_mu       = mu;

%direct_h_info
params.direct_h_index = 1;  %which index are we directly measuring. 
   
 g_fun = @ (k,x) vanderpol_g(k,x,params);
 h_fun = @(k,x) direct_h(k,x,params);

% SASHA: changed .1 to 1
sigma_e = sqrt(.1);           % Standard deviation of initial state estimate
rand('seed', 4321);     % Random number generator seed
% Level of tracing during optimization
if draw_plot
    level = 1;
else
    level = 0;
end



if (strcmp(noise,'student'))
    sigma_r = sqrt(2);
end

% Plotting parameters
h_min   = 0;       % minimum horizontal value in plots
h_max   = T;       % maximum horizontal value in plots
v_min   = -5.0;    % minimum vertical value in plots
v_max   = 5.0;    % maximum vertical value in plots



% ---------------------------------------------------------------
% Rest of the information required by ckbs_nonlinear
%
% Step size for fourth order Runge-Kutta method in vanderpol_sim
% is the same as time between measurement points
step    = T / (N-1);
% Simulate the true values for the state vector
x_vdp  = vanderpol_sim(mu, xi, N, step);
time    = linspace(0., T, N);

% x_trueLong = zeros(n, 10*N);
% x_true = zeros(n, N);
% vanderpol_info.N = 10*N;
% x_trueLong(:,1) = xi; % + sigma_e*randn(n, 1);
% for k = 2:10*N
%     x_trueLong(:,k) = vanderpol_g(k, x_trueLong(:,k-1)) + sqrt(.1*sigma_q)*randn(n,1);
% end
% x_true = x_trueLong(:, [1 10:10:10*N-10]);
% vanderpol_info.N = N;

 x_true = zeros(n, N);
 x_true(:,1) = xi;% + .1*sigma_e*randn(n, 1);
for k = 2:N
     x_true(:,k) = g_fun(k, x_true(:,k-1)) + sigma_qa*randn(n,1);
 end

% Inverse covariance of the measurement and transition noise
rinv    = zeros(m, m, N);
qinv    = zeros(n, n, N);
for k = 1 : N
    if(mod(k,4)==0)
        rinv(:, :, k) = eye(m, m) / sigma_r^2;
    end
    qinv(:, :, k) = eye(n, n) / sigma_q^2;
end
qinv(:, :, 1) = eye(n, n) / sigma_e^2;
% call back functions

f_fun = @ no_f;

v_true = zeros(m,N);
temp = zeros(m,N);
bin = zeros(n,N);

nnzMeas = (mod(1:N, 4) == 0); 
    
switch(noise)
    case{'normal'}
        v_true  = sigma_r * randn(m, N);

    case{'laplace'}
        temp = sigma_r*exprnd(1/sigma_r, m, N);
        bin  = 2*(binornd(1, 0.5, m, N)-0.5);
        v_true = temp.*bin/sqrt(2);

    case{'student'}
        v_true = trnd(4,m,N);
end

% Create contaminating distribution

nums = rand(1, N);
biNums = nums < conFreq;
%    nI= ones(m,1)*binornd(1, conFreq, 1, N); % contaminate entire measurement, not just componentwise
nI = ones(m,1)*biNums;
newNoise = sqrt(conVar)*randn(m,N); % Contaminating noise

%M = 50;

%newNoise = unidrnd(M, 1, N) - M/2;

% Simulate the measurement values
%measurement = x_true(1:m,:).*x_true(1:m, :); % measurement is x^2
measurement =x_true(1,:); % measurement is x^2

if(draw_plot)
    z_Nom =  measurement + v_true;
end
z       =  measurement + v_true + nI.*newNoise;
%addNoise = v_true + nI.*newNoise;
%z       = addNoise + x_true(1, :) + 0.5*(x_true(2,:) - 0).^2  ;
%z1 = addNoise(1, :) + x_true(1, :);s
%z2 = addNoise(2, :) + 0.5*x_true(1, :).*x_true(1, :) + 0.5*x_true(2,:).*x_true(2,:);
%z = [z1; z2];
%z = addNoise + (1/3)^3*(x_true(2,:)).^3 + x_true(2);


% ---------------------------------------------------------------
% call the optimizer
%addpath('../src');

if(draw_plot)
    ckbs_level = 1;
    [ x_out_Gauss , u_out, info ] = ckbs_nonlinear( ...
        f_fun ,    ...
        g_fun ,    ...
        h_fun ,    ...
        max_itr ,  ...
        epsilon ,  ...
        x_in ,     ...
        z_Nom ,        ...
        qinv ,     ...
        rinv ,     ...
        ckbs_level      ...
        );


end
    % ----------------------------------------------------------------------
    % Check that x_out is optimal
    % (this is an unconstrained case, so there is no need to check f).
    % ok     = size(info, 1) <= max_itr;
    % g_out  = zeros(n,   N);
    % h_out  = zeros(m,   N);
    % dg_out = zeros(n, n,   N);
    % dh_out = zeros(m, n,   N);
    % xk     = zeros(n, 1);
    % for k = 1 : N
    %     xkm      = xk;
    %     %
    %     xk       = x_out_Gauss(:, k);
    %     uk       = u_out(:, k);
    %     [gk, Gk] = g_fun(k, xkm);
    %     [hk, Hk] = h_fun(k, xk);
    %     %
    %     g_out(:, k) = gk - xk;
    %     h_out(:, k) = hk;
    %     dg_out(:,:, k) = Gk;
    %     dh_out(:,:, k) = Hk;
    % end
    % dx    = zeros(n, N);
    % grad  = ckbs_sumsq_grad(dx, z, g_out, h_out, dg_out, dh_out, qinv, rinv);
    % ok    = max( max( abs(grad) ) ) < epsilon;
    
    % ----------------------------------------------------------------------
    
    
    
 norm(x_out_Gauss - x_true, 2)




% % ----------------------------------------------------------------------
if draw_plot
   
    
    figure(1);
    
    clf    
    hold on
  %  plot(time, x_vdp(1,:), 'r-.', 'Linewidth', 2, 'MarkerEdgeColor', 'k'); 
    plot(time, x_true(1,:), 'k-', 'Linewidth', 3.5, 'MarkerEdgeColor', 'k');
    plot(time, x_out_Gauss(1,:), 'b--', 'Linewidth', 3.5, 'MarkerEdgeColor', 'k');
   
    
    
    meas = z(1,:).*nnzMeas;
    meas(meas > v_max) = v_max;
    meas(meas < v_min) = v_min;
    meas(meas == 0) = v_max + 2; % don't plot these

    %plot(time, z(1,:), 'o');
    plot(time, meas, 'bo',  'MarkerFaceColor', [.49 1 .63], 'MarkerSize', 6.2);
    
    axis([h_min, h_max, v_min, v_max]);

    %hleg = legend('vdp solution', 'stochastic realization', 'Gaussian Smoother', 'L1 Smoother', 'Data', 'Location', 'NorthWest');
    xlabel('Time (s)', 'FontSize', 16);
    ylabel('X_1-component', 'FontSize', 16);
%    set(hleg, 'Box', 'off');
    set(gca,'FontSize',16);
   % set(hleg, 'FontSize', 14);
   %     title('20% of measurement errors are N(0, 100)');

    hold off
  %  savefig(1, [pwd '/VP1']);

    %plot(time, z(1,:), 'o');
  
    figure(2)
    
    hold on
    
  %  plot(time, x_vdp(2,:), 'r-.', 'Linewidth', 2, 'MarkerEdgeColor', 'k'); 
   plot(time, x_true(2,:), 'k-', 'Linewidth', 3.5, 'MarkerEdgeColor', 'k');
    plot(time, x_out_Gauss(2,:), 'b--', 'Linewidth', 3.5, 'MarkerEdgeColor', 'k');
  %  plot(time, x_out_L1(2,:), 'm-.', 'Linewidth', 3.5, 'MarkerEdgeColor', 'k');
    axis([h_min, h_max, v_min, v_max]);
    
 %hleg = legend('stochastic realization', 'Gaussian Smoother',  'Location', 'NorthWest');
 xlabel('Time(s)', 'FontSize', 16);
 ylabel('X_2-component', 'FontSize', 16);
 %set(hleg, 'Box', 'off');
 set(gca,'FontSize',16);
 %set(hleg, 'FontSize', 14);
 
   
    hold off
    
%     savefig(2, [pwd '/VP2']);
   
    
 
    
    
    return
end

end

% $$ $newlinech %$$
% $end
