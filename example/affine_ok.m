% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006
% Authors: Bradlely Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
%          Aleksandr Aravkin:    saravkin@eos.ubc.ca
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin affine_ok.m$$ $newlinech %$$
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
% $$
%
% $section ckbs_affine Smoothing Spline Example and Test$$
%
% $index smoothing, spline example and test$$
% $index example, smoother spline$$
% $index test, smoothing spline$$
%
% $index ckbs_affine, example and test$$
% $index affine, example and test$$
% $index example, affine$$
% $index test, affine$$
%
% $head Syntax$$
% $codei%[%ok%] = affine_ok(%draw_plot%)%$$
%
% $head draw_plot$$
% If this argument is true, a plot is drawn showing the results
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
function [ok] = affine_ok(draw_plot)
    % --------------------------------------------------------
    % You can change these parameters
    N     = 100;        % number of measurement time points
    dt    = 4*pi / N;  % time between measurement points
    gamma =  1;        % transition covariance multiplier
    sigma =  .35;       % standard deviation of measurement noise
    max_itr = 30;      % maximum number of iterations
    epsilon = 1e-5;    % convergence criteria
    h_min   = 0;       % minimum horizontal value in plots
    h_max   = 7;       % maximum horizontal value in plots
    v_min   = -2.0;    % minimum vertical value in plots
    v_max   = +2.0;    % maximum vertical value in plots
    % ---------------------------------------------------------
    ok = true;
    if nargin < 1
        draw_plot = false;
    end
    %  Define the problem
    rand('seed', 1234);
    %
    % number of constraints per time point
    ell   = 4;
    %
    % number of measurements per time point
    m     = 1;
    %
    % number of state vector components per time point
    n     = 2;
    %
    % simulate the true trajectory and measurement noise
    t       =  (1 : N) * dt;
    x1_true = - cos(t);
    x2_true = - sin(t);
    x_true  = [ x1_true ; x2_true ];
    v_true  = sigma * rand(1, N);
    %
    % measurement values and model
    v_true  = sigma * randn(1, N);
    z       = x2_true + v_true;
    rk      = sigma * sigma;
    rinvk   = 1 / rk;
    rinv    = zeros(m, m, N);
    h       = zeros(m, N);
    dh      = zeros(m, n, N);
    for k = 1 : N
        rinv(:, :, k) = rinvk;
        h(:, k)       = 0;
        dh(:,:, k)    = [ 0 , 1 ];
    end
    %
    % transition model
    g       = zeros(n, N);
    dg      = zeros(n, n, N);
    qinv    = zeros(n, n, N);
    qk      = gamma * [ dt , dt^2/2 ; dt^2/2 , dt^3/3 ];
    qinvk   = inv(qk);
    for k = 2 : N
        g(:, k)       = 0;
        dg(:,:, k)    = [ 1 , 0 ; dt , 1 ];
        qinv(:,:, k)  = qinvk;
    end
    %
    % initial state estimate
    g(:, 1)      = x_true(:, 1);
    qinv(:,:, 1) = 100 * eye(2);
    %
  % constraints

    b       = zeros(0, N);

    db      = zeros(0, n, N);


    % --------------------------------------------------------------------
    %
    % -------------------------------------------------------------------------
    [xOut, uOut, info] = ...
        ckbs_affine(max_itr, epsilon, z, b, g, h, db, dg, dh, qinv, rinv);
    % --------------------------------------------------------------------------
    ok   = ok & all( info(end, 1:3) <= epsilon);
    d    = ckbs_sumsq_grad(xOut, z, g, h, dg, dh, qinv, rinv);
    for k = 1 : N
        xk = xOut(:, k);
        uk = uOut(:, k);
        bk = b(:, k);
        Bk = db(:,:, k);
        dk = d(:, k);
        sk = - bk - Bk * xk;
        %
        ok = ok & (min(uk) >= 0.);
        ok = ok & (max (bk + Bk * xk) <= epsilon);
        ok = ok & (max ( abs( Bk' * uk + dk ) ) <= epsilon);
        ok = ok & (max ( uk .* sk ) <= epsilon );
    end
    if draw_plot
        figure(1);
        clf
        hold on
        plot(t', x_true(2,:)', 'r-' );
        plot(t', z(1,:)', 'ko' );
        plot(t', xOut(2,:)', 'b--' );
   %     plot(t', - ones(N,1), 'b-');
  %      plot(t', ones(N,1), 'b-');
        axis([h_min, h_max, v_min, v_max]);
        title('Affine smoother');
        %legend('Truth', 'Measurements', 'Estimate', 'Location', 'SouthEast' );
        hold off
        %
        % constrained estimate
        x_con = xOut;
    end
    savefig(1, [pwd '/Experiments/affine']);
    %
    
end
% $$ $newlinech %$$
% $end
