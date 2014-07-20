% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2011
% Authors: Bradlely Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
%          Aleksandr Aravkin:    saravkin at eos dot ubc dot ca
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin affine_singular_ok.m$$ $newlinech %$$
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
%       singular Singular
%       iter
%       alg
%       end end       
% $$
%
% $section ckbs_affine_singular Singular Smoothing Spline Example and Test$$
%
% $index singular smoothing, spline example and test$$
% $index example, singular smoother spline$$
% $index test, singular smoothing spline$$
%
% $index ckbs_affine_singular, example and test$$
% $index singular affine, example and test$$
% $index example, singular affine$$
% $index test, singular affine$$
%
% $head Syntax$$
% $codei%[%ok%] = affine_singular_ok(%draw_plot%)%$$
%
% $head draw_plot$$
% If this argument is true, a plot is drawn showing the results
% comparing $cref ckbs_affine_singular$$ with $cref ckbs_affine$$.
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
function [ok] = affine_singular_ok(draw_plot)

    % --------------------------------------------------------

    ok = true;
    % You can change these parameters

    conVar = 100;         % variance of contaminating distribution

    conFreq = .3;         % frequency of outliers

    noise = 'normal';     % type of noise

    N     = 1000;        % number of measurement time points

    dt    = 10*pi / N;  % time between measurement points

    %dt = 1;

    gamma =  1;        % transition covariance multiplier

    sigma =  1;       % standard deviation of measurement noise
    sigma_alg = 1;

    if (strcmp(noise,'student'))
        sigma = sqrt(2);
    end


    max_itr = 30;      % maximum number of iterations

    epsilon = 1e-4;    % convergence criteria

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

    v_true = zeros(m,N);
    temp = zeros(m,N);
    bin = zeros(m,N);

    randn('state', sum(100*clock));
    rand('twister', sum(100*clock));


    if (strcmp(noise,  'normal'))
        v_true  = sigma * randn(m, N);
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



    z = [x2_true] + v_true;

    rk      = sigma_alg * sigma_alg * eye(m);

    rinvk   = inv(rk);

    rinv    = zeros(m, m, N);
    r       = zeros(m, m, N);

    h       = zeros(m, N);

    % Covariance between process and measurement?
    dh      = zeros(m, n, N);

    for k = 1 : N

        rinv(:, :, k) = rinvk;
        r(:,:,k) = rk;

        h(:, k)       = 0;

        dh(:,:, k)    = [ 0 , 1 ];

    end

    %

    % transition model

    g       = zeros(n, N);

    dg      = zeros(n, n, N);

    qinv    = zeros(n, n, N);
    q       = zeros(n, n, N);

    qk      = gamma * [ dt , 0 ; 0 , 0 ];

    
    qinvk   = inv(qk);

    for k = 2 : N

        g(:, k)       = 0;

        dg(:,:, k)    = [ 1 , 0 ; dt , 1 ];

        qinv(:,:, k)  = qinvk;
        
        q(:,:,k) = qk;

    end

    %

    % initial state estimate

    g(:, 1)      = x_true(:, 1);

    qinv(:,:, 1) = 100 * eye(2);
    
    q(:,:,1) = inv(qinv(:,:,1));

    %

    % constraints

    b       = zeros(0, N);

    db      = zeros(0, n, N);


    % -------------------------------------------------------------------------


    % -------------------------------------------------------------------------

    % L1 Kalman-Bucy Smoother
    % Linear unconstrained Kalman-Bucy Smoother
    [xOut2, info] = ckbs_affine_singular(z, g, h, dg, dh, q, r);
    
    %qk      = gamma * [ dt , dt^2/2 ; dt^2/2 , dt^3/3 ];

    %for k = 2:N
    %    q(:,:,k) = qk;
    %end
    
    %[xOut3, info] = ckbs_affine_singular(z, g, h, dg, dh, q, r);
    
    iter = info(5);
    info(5) = 0;

    iter
    ok = (ok) && (max(info) < epsilon);


    
    
    

    % % --------------------------------------------------------------------------
    %
    if draw_plot

        qk      = gamma * [ dt , dt^2/2 ; dt^2/2 , dt^3/3 ];
        
        qinvk   = inv(qk);

        for k = 2 : N
            qinv(:,:, k)  = qinvk;       
        end
        
        
        [xOut3, uOut, info] = ckbs_affine(max_itr, epsilon, z, b, g, h, db, dg, dh, qinv, rinv);
        
        figure(1);

        clf

        hold on

        plot(t', x_true(2,:)', '-k', 'Linewidth', 2, 'MarkerEdgeColor', 'k' );

        meas = z(1,:);
        meas(meas > v_max) = v_max;
        meas(meas < v_min) = v_min;
        plot(t', meas, 'bd',  'MarkerFaceColor', [.49 1 .63], 'MarkerSize', 6);


        plot(t', xOut2(2,:)', '-r', 'Linewidth', 2, 'MarkerEdgeColor', 'k');
%        plot(t', xOut3(2,:)', '-b', 'Linewidth', 2, 'MarkerEdgeColor', 'k');


        axis([h_min, h_max, v_min, v_max]);

        title('Position Plot. Black=truth, Blue = meas, Magenta = Singular Smoother.');

        hold off

        figure(2);
        clf

        hold on

        plot(t', x_true(1,:)', '-k', 'Linewidth', 2, 'MarkerEdgeColor', 'k' );



        plot(t', xOut2(1,:)', '-r', 'Linewidth', 2, 'MarkerEdgeColor', 'k');

        plot(t', xOut3(1,:)', '-b', 'Linewidth', 2, 'MarkerEdgeColor', 'k');


        axis([h_min, h_max, v_min, v_max]);

        title('Derivative plot. Black=truth, Blue = SDE smoother, Magenta = Singular Smoother.');

        
        %

        % constrained estimate

    end
   % [norm(xOut2(2,:) - x_true(2,:)) norm(xOut3(2,:) - x_true(2,:))]
end
% $$ $newlinech %$$
% $end
