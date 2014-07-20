% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2011
% Authors: Bradlely Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
%          Aleksandr Aravkin:    saravkin at eos dot ubc dot ca
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin L1_affine_ok.m$$ $newlinech %$$
% $spell
%       itr
%       ko
%       nargin
%       clf
%       fid
%       fopen
%       wt
%       fprintf
%       fclose
%       disp
%       diff
%       ckbs
%       bk
%       cos
%       dg
%       dh
%       dk
%       dt
%       inv
%       qinv
%       qinvk
%       qk
%       randn
%       rinv
%       rinvk
%       rk
%       sk
%       sumsq
%       uk
%       xk
%       Var
%       Freq
%       strcmp
%       sqrt
%       exp
%       binornd
%       exprnd
%       trnd
%       Linewidth
%       bd
%       end end       
% $$
%
% $section ckbs_L1_affine Robust Smoothing Spline Example and Test$$
%
% $index robust smoothing, spline example and test$$
% $index example, robust smoother spline$$
% $index test, robust smoothing spline$$
%
% $index ckbs_L1_affine, example and test$$
% $index robust affine, example and test$$
% $index example, robust affine$$
% $index test, robust affine$$
%
% $head Syntax$$
% $codei%[%ok%] = L1_affine_ok(%draw_plot%)%$$
%
% $head draw_plot$$
% If this argument is true, a plot is drawn showing the results
% comparing $cref ckbs_L1_affine$$ with $cref ckbs_affine$$.
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
% $head Source Code$$
%
% $newlinech $$ $codep
function [ok] = L1_affine_ok(draw_plot)

    % --------------------------------------------------------

    ok = true;
    % You can change these parameters

    conVar = 100;         % variance of contaminating distribution

    conFreq = .3;         % frequency of outliers

    noise = 'normal';     % type of noise

    N     = 30;        % number of measurement time points

    dt    = 4*pi / N;  % time between measurement points

    %dt = 1;

    gamma =  1;        % transition covariance multiplier

    sigma =  0.5;       % standard deviation of measurement noise

    if (strcmp(noise,'student'))
        sigma = sqrt(2);
    end


    max_itr = 30;      % maximum number of iterations

    epsilon = 1e-5;    % convergence criteria

    h_min   = 0;       % minimum horizontal value in plots

    h_max   = 14;       % maximum horizontal value in plots

    v_min   = -5.0;    % minimum vertical value in plots

    v_max   = 5.0;    % maximum vertical value in plots

    % ---------------------------------------------------------


    if nargin < 1

        draw_plot = false;

    end

    %  Define the problem

    %rand('seed', 12);

    %

    % number of constraints per time point

    % ell   = 0;

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

    % measurement values and model


    %exp_true = random('exponential', 1/sigma, 1, N);
    %bin = 2*random('binomial', 1, 0.5, 1, N) - 1;
    %exp_true = exp_true .* bin;
    %z        = x2_true + exp_true;


    % Normal Noise

    v_true = zeros(1,N);
    temp = zeros(1,N);
    bin = zeros(1,N);

    randn('state', sum(100*clock));
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

    
    temp = rand(1, N);
    nI = temp < conFreq;
    newNoise = sqrt(conVar)*randn(1,N);


    z = x2_true + v_true + nI.*newNoise;

    rk      = sigma * sigma;

    rinvk   = 1 / rk;

    rinv    = zeros(m, m, N);

    h       = zeros(m, N);

    % Covariance between process and measurement?
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


    % -------------------------------------------------------------------------

    % Linear unconstrained Kalman-Bucy Smoother
    [xOut, uOut, info] = ckbs_affine(max_itr, epsilon, z, b, g, h, db, dg, dh, qinv, rinv);

    % -------------------------------------------------------------------------

    % L1 Kalman-Bucy Smoother
    [xOut2, rOut, sOut, pPlusOut, pMinusOut, info] = ckbs_L1_affine(max_itr, epsilon, z, g, h, dg, dh, qinv, rinv);

    ok = (ok) && (info(end, 2) < epsilon);


    % % --------------------------------------------------------------------------
    %
    if draw_plot

        figure(1);

        clf

        hold on

        plot(t', x_true(2,:)', '-k', 'Linewidth', 2, 'MarkerEdgeColor', 'k' );

        meas = z(1,:);
        meas(meas > v_max) = v_max;
        meas(meas < v_min) = v_min;
        plot(t', meas, 'bd',  'MarkerFaceColor', [.49 1 .63], 'MarkerSize', 6);


        plot(t', xOut(2,:)', '-m', 'Linewidth', 2, 'MarkerEdgeColor', 'k' );

        plot(t', xOut2(2,:)', '-r', 'Linewidth', 2, 'MarkerEdgeColor', 'k');



        axis([h_min, h_max, v_min, v_max]);

        title('Gaussian and L1 smoothers. Black=truth, Blue = meas, Magenta = Gaussian, Red = L1.');

        hold off

        %

        % constrained estimate

    end
end
% $$ $newlinech %$$
% $end
