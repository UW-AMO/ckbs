% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006
% Authors: Bradlely Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin affine_ok_box.m$$ $newlinech %$$
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
% $section ckbs_affine Box Constrained Smoothing Spline Example and Test$$
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
% $codei%[%ok%] = affine_ok_box(%draw_plot%)%$$
%
% $head draw_plot$$
% If this argument is true, a plot is drawn showing the results
% and the $code nonlinear_ok_box.out$$ file is written for use with 
% the program $code nonlinear_ok_box.r$$.
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
% $head Constraint$$
% $latex \[
%	\begin{array}{c}
%	-1 \leq x1 (t) \leq +1
%	\\
%	-1 \leq x2 (t) \leq +1
%	\end{array}
% \]$$.
%
% $head Source Code$$
%
% $newlinech $$ $codep
function [ok] = tridiag_solver_break(length)
    % --------------------------------------------------------
    % You can change these parameters
    N     = length;        % number of measurement time points
    dt    = 2*pi /N;  % time between measurement points
    gamma =  1;        % transition covariance multiplier
    sigma =  .5;       % standard deviation of measurement noise
    max_itr = 30;      % maximum number of iterations
    epsilon = 1e-5;    % convergence criteria
    h_min   = 0;       % minimum horizontal value in plots
    h_max   = 7;       % maximum horizontal value in plots
    v_min   = -10.0;    % minimum vertical value in plots
    v_max   = +10.0;    % maximum vertical value in plots
    % ---------------------------------------------------------
    ok = true;
%    if nargin < 3
%        draw_plot = false;
%    end
draw_plot = true;
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
   
    %rinvk   = 1 / rk;
    
    rinv    = zeros(m, m, N);
    h       = zeros(m, N);
    dh      = zeros(m, n, N);
    for k = 1 : N
    %    rinv(:, :, k) = rinvk;
        h(:, k)       = 0;
        dh(:,:, k)    = [ 0 , 1 ];
    end
    %
    % transition model
    g       = zeros(n, N);
    dg      = zeros(n, n, N);
    qinv    = zeros(n, n, N);
    %qk      = gamma * [ dt , dt^2/1.75 ; dt^2/1.75 , dt^3/3 ];
    qk      = gamma * [ dt , dt^2/2 ; dt^2/2 , dt^3/3 ];

    qinvk   = inv(qk);
    for k = 2 : N
        g(:, k)       = 0;
        dg(:,:, k)    = [ 1.1 , 0 ; dt , 1.1 ];
        qinv(:,:, k)  = qinvk;
    end
    %
    % initial state estimate
    g(:, 1)      = x_true(:, 1);
    qinv(:,:, 1) =  .01* eye(2);
    %
    
    % Unconstrained Case
    b           = zeros(0, N);
    db          = zeros(0, n, N);
    % -------------------------------------------------------------------------
disp(sprintf('\nForward solver min eigenvalues:\n'));
    [xOut, uOut, info] = ...
        ckbs_affine(max_itr, epsilon, z, b, g, h, db, dg, dh, qinv, rinv);
disp(sprintf('\nBackward solver min eigenvalues:\n'));
    [xOutb, uOut, info] = ...
        ckbs_affine_b(max_itr, epsilon, z, b, g, h, db, dg, dh, qinv, rinv);


    % --------------------------------------------------------------------------
    ok   = ok & all( info(end, 1:3) <= epsilon);
    d    = ckbs_sumsq_grad(xOut, z, g, h, dg, dh, qinv, rinv);
    for k = 1 : N
        xk = xOut(:, k);
        dk = d(:, k);
        %
        ok = ok & (min(dk) <= epsilon);
    end
    if draw_plot
        figure(1);
        clf
        hold on
        plot(t', x_true(2,:)', 'r-' );
        plot(t', z(1,:)', 'ko' );
        plot(t', xOut(2,:)', 'b-' );
        plot(t', xOutb(2,:)', 'm-' );
        axis([h_min, h_max, v_min, v_max]);
        title('Forward');
        hold off
                figure(2);
        clf
        hold on
        plot(t', x_true(2,:)', 'r-' );
        plot(t', z(1,:)', 'ko' );
        plot(t', xOutb(2,:)', 'b-' );
        axis([h_min, h_max, v_min, v_max]);
        title('Backward');
        hold off

        %
        % unconstrained estimate
        x_free = xOut;
        %
        % write out constrained and unconstrained results
        [fid, msg] = fopen('affine_ok_box.out', 'wt');
        if size(msg, 2) > 0 
            disp(['affine_ok: ', msg]);
        end
        %                      123456789012345678901234'
        heading =             '       t';
        heading = [ heading , ' x2_true  x2_con x2_free' ];
        heading = [ heading , '      z1\n'               ];
        fprintf(fid, heading);
        for k = 1 : N
            fprintf(fid,'%8.3f%8.3f%8.3f%8.3f\n', ...
                t(k), ...
                x_true(2,k),x_free(2,k), ...
                z(1,k) ...
            );
        end
        fclose(fid);
    end
    return
end
% $$ $newlinech %$$
% $end
