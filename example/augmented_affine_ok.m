% augmented affine test
% creates plots (when draw_plot set to 1) for the MLSP 2014 submission. 
% Contributed by Karthikeyan Natesan Ramamurthy

function [ok] = augmented_affine_ok(draw_plot)
    % --------------------------------------------------------
    % You can change these parameters
    N     = 100;        % number of measurement time points
    dt    = 4*pi / N;  % time between measurement points
    gamma =  1;        % transition covariance multiplier
    sigma =  .55;       % standard deviation of measurement noise
    max_itr = 30;      % maximum number of iterations
    epsilon = 1e-5;    % convergence criteria
    h_min   = 0;       % minimum horizontal value in plots
    h_max   = 7;       % maximum horizontal value in plots
    v_min   = -2.0;    % minimum vertical value in plots
    v_max   = +6.0;    % maximum vertical value in plots
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
    
    bias = 1; 
    Pone    =  1; % measure bias directly

    z       = x2_true + v_true + Pone*bias;
    
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
    [xOutOld, uOut, info] = ...
        ckbs_affine(max_itr, epsilon, z, b, g, h, db, dg, dh, qinv, rinv);
    [xOutNew, yOut] = ...
        ckbs_augmented_affine(max_itr, epsilon, z, b, g, h, db, dg, dh, qinv, rinv, Pone);

    fprintf('yOut is %5.1f\n', yOut);

    if draw_plot
        figure(1);
        clf
        hold on
        plot(t', x_true(2,:)', 'k-', 'linewidth', 1.5 );
        plot(t', z(1,:)', 'ko' );
        plot(t', xOutOld(2,:)', 'b--', 'linewidth', 1.5 );
        plot(t', xOutNew(2,:)', 'r-.', 'linewidth', 1.5 );
   %     plot(t', - ones(N,1), 'b-');
  %      plot(t', ones(N,1), 'b-');
        axis([h_min, h_max, v_min, v_max]);
        title('Affine smoother');
        %legend('Truth', 'Measurements', 'Naive Estimate', 'Debiased Estimate', 'Location', 'SouthEast' );
        hold off
        %
        % constrained estimate
%        x_con = xOut;
    end
    savefig(1, ['/Users/saravkin/Dropbox/KalmanCovariates/tex/figures/sine']);
    %
    
end
% $$ $newlinech %$$
% $end
