% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2013
% Authors: Aleksandr Y. Aravkin: saravkin at us dot ibm dot com
%          Bradley Bell:         bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin t_general_noisy_jump.m$$ $newlinech %$$
% $spell
%	itr
%	ko
%	nargin
%	clf
%	fid
%	fopen
%	wt
%	fprintf
%	fclose
%	disp
%	diff
%	ckbs
%	bk
%	cos
%	dg
%	dh
%	dk
%	dt
%	inv
%	qinv
%	qinvk
%	qk
%	randn
%	rinv
%	rinvk
%	rk
%	sk
%	sumsq
%	uk
%	xk
%   Var
%   Freq
%   strcmp
%   sqrt
%   params
%   derivs
%   cumsum
%   params
%   pos
%   vel
%   df
%   inds
%   exp
%   trnd
%   binornd
%   exprnd
%   end end
%   Linewidth
%   bd
%   rs
%   xfilt
%   xsmooth
%   xlabel
%   ylabel
%   Fontsize
%   gca
%   ind
%   rvec
%   res
%   xkm
%   feval
%   gk 
%   res 
% $$
%
% $section ckbs_t_general Jump Tracking Example and Test$$
%
% $index smoothing, jump tracking example and test$$
% $index example, smoothing and jump tracking$$
% $index test, smoothing spline with jump$$
%
% $index ckbs_t_general, example and test$$
% $index t_general, example and test$$
% $index example, t_general$$
% $index test, Student's t$$
%
% $head Syntax$$
% $codei%[%ok%] = affine_ok(%draw_plot%, %conVar%, %conFreq%, %proc_ind%,
% %meas_ind%)%$$
%
% $head draw_plot$$
% If this argument is true, a plot is drawn showing the results
%
% $head conVar$$
% This argument gives the variance of the contaminating signal. 
% The larger the variance, the bigger the outliers. 
% 
% $head conFreq$$
% This argument gives outlier frequency. For example, .1 means 1/10
% measurements will be an outlier. 
% 
% $head proc_ind$$
% This selects which process indices the smoother will model using Student's t. 
% The choices are [], 1, 2, and 1:2. 
%
% $head meas_ind$$
% This selects which measurement indices the smoother will model using Student's t. 
% The choices are [] and 1. 
%
% $head ok$$
% If the return value $icode ok$$ is true, the test passed,
% otherwise the test failed.
%
% $head State Vector$$
% $table
% $latex x1 (t)$$ $cnext derivative of function we are estimating $rnext
% $latex x2 (t)$$ $cnext value of function we are estimating
% $tend
%
% $head Measurement$$
% $latex z1 (t)$$ value of $latex x2 (t)$$ plus noise
%
%
% $head Source Code$$
%
% $newlinech $$ $codep


function [ok]= t_general_noisy_jump(draw_plot, conVar, conFreq)

noise = 'normal';
rand('seed', 1234);
randn('seed', 1234);

% --------------------------------------------------------

N     = 100;        % number of measurement time points

dt    = 4*pi / N;  % time between measurement points

%dt = 1;

gamma =  1;        % transition covariance multiplier

sigma =  0.5;       % standard deviation of measurement noise

if (strcmp(noise,'student'))
    sigma = sqrt(2);
end


max_itr = 100;      % maximum number of iterations

epsilon = 1e-5;    % convergence criteria

h_min   = 0;       % minimum horizontal value in plots

h_max   = 14;       % maximum horizontal value in plots

v_min   = -7.0;    % minimum vertical value in plots

v_max   = 7.0;    % maximum vertical value in plots

% ---------------------------------------------------------


if nargin < 1

     draw_plot = false;

end


% number of measurements per time point
m     = 1;
% number of state vector components per time point
n     = 2;

params.level = 1;



% simulate the true trajectory and measurement noise

t1       =  (1 : N/2 - 2) * dt;
t2       =  (N/2 - 1: N/2 + 1)*dt;
t3       =  (N/2 +2 : N)*dt;

t = [t1, t2, t3];

% alpha = 8;
% beta = 4;
% gammaOne = 20;

jump = 10;
derivs = jump/(4*dt) * ones(1,3);

    temp1 = [-cos(t1), derivs, -cos(t3)];
    temp2 = [-jump/2-sin(t1), -jump/2-sin(t2) + cumsum(derivs)*dt, jump/2-sin(t3)];

   x1_true = temp1(1:N);
   x2_true = temp2(1:N);



% velocity is in x_1
% position is in x_2
x_true  = [ x1_true ; x2_true ];

% direct measurement on position
params.direct_h_index = 2;
h_fun = @(k,x) direct_h(k,x,params);

params.pos_vel_g_dt = dt;  %why was this set to .01??
params.pos_vel_g_initial = x_true(:,1); % initialize to true values for now
g_fun = @(k,x) pos_vel_g(k,x,params);

df_meas = 4;
df_proc = 4;
params.df_proc = df_proc;
params.df_meas = df_meas;

% measurement values and model


%exp_true = random('exponential', 1/sigma, 1, N);
%bin = 2*random('binomial', 1, 0.5, 1, N) - 1;
%exp_true = exp_true .* bin;
%z        = x2_true + exp_true;


% Normal Noise

v_true = zeros(1,N);
%temp = zeros(1,N);
%bin = zeros(1,N);

randn('state', sum(100*clock))
rand('twister', sum(100*clock));


if (strcmp(noise,  'normal'))
    v_true  = sigma * randn(1, N);
end


%Laplace Noise
if (strcmp(noise,'laplace'))
    temp = sigma*exprnd(1/sigma, 1,100);
    bin  = 2*(binornd(1, 0.5, 1, 100)-0.5);
    v_true = temp.*bin/sqrt(2);
end

if (strcmp(noise,'student'))
  v_true = trnd(4,1,100);
end

nI= binornd(1, conFreq, 1, N);
newNoise = sqrt(conVar)*randn(1,N);


% measurements obtained from second coordinate
z = x2_true + v_true + nI.*newNoise;

rk      = sigma * sigma;

rinvk   = 1 / rk;

rinv    = zeros(m, m, N);


for k = 1 : N
    if (mod(k, 3) ~= 0)
        rinv(:, :, k) = rinvk;
   end
end

% transition model

qinv    = zeros(n, n, N);
qk      = gamma * [ dt , dt^2/2 ; dt^2/2 , dt^3/3 ];
%qk       = gamma * [ dt^3/3 , dt^2/2 ; dt^2/2 , dt ];
qinvk   = inv(qk);

for k = 2 : N
     qinv(:,:, k)  = qinvk;
end

qinv(:,:,1) = 100*eye(n);
%


xZero = zeros(n, N);
% L1 Kalman-Bucy Smoother

w = 1;


params.inds_proc_st = [];
params.inds_meas_st = [];
[xOut, info] = ckbs_t_general(g_fun, h_fun, max_itr, epsilon, xZero, z, qinv, w*rinv, params);

params.inds_proc_st = 1:2;
params.inds_meas_st = [];
[xOutStProc, infoStProc] = ckbs_t_general(g_fun, h_fun, max_itr, epsilon, xZero, z, qinv, w*rinv, params);




params.inds_proc_st = [];
params.inds_meas_st = 1;
[xOutStMeas, infoStMeas] = ckbs_t_general(g_fun, h_fun, max_itr, epsilon, xZero, z, qinv, w*rinv, params);

params.inds_proc_st = 1:2;
params.inds_meas_st = 1;
[xOutTT, infoTT] = ckbs_t_general(g_fun, h_fun, max_itr, epsilon, xZero, z, qinv, w*rinv, params);

% params.inds_proc_st = 2;
% params.inds_meas_st = 1;
% [xOutStPosTT, infoStPosTT] = ckbs_t_general(g_fun, h_fun, max_itr, epsilon, xZero, z, qinv, w*rinv, params);



% % --------------------------------------------------------------------------
% 
if draw_plot
    msize = 4;
    lwdThin = 2;
    lwdThick = 3;
    meas = z(1,:);
    meas(meas > v_max) = v_max;
    meas(meas < v_min) = v_min;
    rvec = rinv(1,:);
    meas(rvec==0) = v_max + 5; % take non-measurements out

    % All Gaussian
  % % --------------------------------------------------------------------------  
    figure(1);

    clf

    hold on

    plot(t', x_true(2,:)', '-k', 'Linewidth', lwdThin, 'MarkerEdgeColor', 'k' );
    plot(t', meas, 'bd',  'MarkerFaceColor', [.49 1 .63], 'MarkerSize', msize);

    plot(t', xOut(2,:)', '--r', 'Linewidth', lwdThick, 'MarkerEdgeColor', 'k' );
    hold off

    axis([h_min, h_max, v_min-.5, v_max+.5]);
    
   % %
   % Student's t Process, gaussian measurements
   % --------------------------------------------------------------------------
    figure(2);
    
    clf

    hold on

    plot(t', x_true(2,:)', '-k', 'Linewidth', lwdThin, 'MarkerEdgeColor', 'k' );


    plot(t', meas, 'bd',  'MarkerFaceColor', [.49 1 .63], 'MarkerSize', msize);
    plot(t', xOutStProc(2,:)', '--r', 'Linewidth', lwdThick, 'MarkerEdgeColor', 'k');
    hold off
    axis([h_min, h_max, v_min-.5, v_max+.5]);
   
    % %
    % Student's t measurements, gaussian process
    % --------------------------------------------------------------------------
    figure(3);

    clf

    hold on

    plot(t', x_true(2,:)', '-k', 'Linewidth', lwdThin, 'MarkerEdgeColor', 'k' );


    plot(t', meas, 'bd',  'MarkerFaceColor', [.49 1 .63], 'MarkerSize', msize);
    plot(t', xOutStMeas(2,:)', '--r', 'Linewidth', lwdThick, 'MarkerEdgeColor', 'k');
    hold off
    axis([h_min, h_max, v_min-.5, v_max+.5]);
    
    % %
    % Student's t all
    % --------------------------------------------------------------------------
    figure(4);

    clf

    hold on

    plot(t', x_true(2,:)', '-k', 'Linewidth', lwdThin, 'MarkerEdgeColor', 'k' );


    plot(t', meas, 'bd',  'MarkerFaceColor', [.49 1 .63], 'MarkerSize', msize);
    plot(t', xOutTT(2,:)', '--r', 'Linewidth', lwdThick, 'MarkerEdgeColor', 'k');
  %  plot(t', xOutStPosTT(2,:)', '-.m', 'Linewidth', lwdThick, 'MarkerEdgeColor', 'k');

    hold off
    axis([h_min, h_max, v_min-.5, v_max+.5]);
    
    % %
    % Student's t position, gaussian velocity and measurements
    % ---------------------------------------------------------------------
%     % -----
%     figure(5);
% 
%     clf
% 
%     hold on
% 
%     plot(t', x_true(2,:)', '-k', 'Linewidth', lwdThin, 'MarkerEdgeColor', 'k' );
% 
% 
%     plot(t', meas, 'bd',  'MarkerFaceColor', [.49 1 .63], 'MarkerSize', msize);
%  %   plot(t', xOutTT(2,:)', '--r', 'Linewidth', lwdThick, 'MarkerEdgeColor', 'k');
%     plot(t', xOutStPosTT(2,:)', '--r', 'Linewidth', lwdThick, 'MarkerEdgeColor', 'k');
% 
%     hold off
%     axis([h_min, h_max, v_min-.5, v_max+.5]);
    

    %   plot(t', xOut5(2,:)', '-b', 'Linewidth', 2, 'MarkerEdgeColor', 'k');

    dir = '/Users/aleksand/Dropbox/ckbs/trunk/results/';
    nameBase = 'NoisyJumpCompare';
    
  %  savefig(1:4, '/Users/aleksand/Dropbox/AnnalsELQP/TStudentIEEEtran/TFigures/NoisyJumpCompare');
       savefig(1:4, [dir nameBase 'var' num2str(conVar) 'freq' num2str(conFreq*100)]);

   %plot(t', xsmooth(2,:)', '-m', 'Linewidth', 2, 'MarkerEdgeColor', 'k');
  
     %plot(t', xfilt(2,:)', '-b', 'Linewidth', 2, 'MarkerEdgeColor', 'k');
 
%     plot(t', - ones(N,1), 'b-');

%     plot(t', ones(N,1), 'b-');

%    xlabel('Time', 14);
 %   ylabel('Position', 'Fontsize', 14);
  %  set(gca,'FontSize',14);
     

  axis([h_min, h_max, v_min, v_max]);

  res_x = zeros(n, N);
  xk = zeros(n, 1);
  for k = 1 : N
      xkm   = xk;
      xk    = xOut(:, k);
      gk    = feval(g_fun, k, xkm);
      res_x(:,k)  = xk - gk;
  end
  
%   hold off
%   figure(2);
%   clf
%   hold on
%   plot(t', x_true(2,:)', '-k', 'Linewidth', 2, 'MarkerEdgeColor', 'k' );
%   plot(t', xOut(2,:), '-m', 'Linewidth', 2, 'MarkerEdgeColor', 'k' );
%   hold off
%   
%   figure(3);
%   clf 
%   hold on
%   
%   plot(t', res_x(2,:)', '-m', 'Linewidth', 2, 'MarkerEdgeColor', 'k' );
% 
%   hold off
% 
%   figure(4)
%   clf
%   hold on
%   plot(t', res_x(2,:)', '-b', 'Linewidth', 2, 'MarkerEdgeColor', 'k' );
%   hold off
  
  
end

% $$ $newlinech %$$
% $end
