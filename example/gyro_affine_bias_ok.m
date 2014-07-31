% gyro affine test with bias estimation
% creates figures for MLSP submission 2014 when draw_plot set to 1
% Contributed by Karthikeyan Natesan Ramamurthy

function [ok] = gyro_affine_bias_ok(draw_plot)
close all; clc;

%% Simulated paramers
Fs = 100; % Sampling frequency of gyoscope (in Hz)
dt = 1/Fs; % Sampling interval
T = 1; % Stop time
t = (0:dt:T-dt)';

% Ground truth ang. vel.
F = 3; % Frequency of sine wave (Ang. Vel.)
av = sin(2*pi*F*t); % Ang. Vel.
figure; plot(t,av); title(['True angular velocity']);
avb = 0.3; % Angular velocity bias

ang = zeros(length(t),1);
% Generate noiseless (true) angles
for i1 = 2:length(t)
    ang(i1) = ang(i1-1)+dt*av(i1-1);
end
figure; plot(t,ang); title(['True orientation']);


%%

N = length(t);
gamma =  1;        % transition covariance multiplier
sigma =  0.1;       % standard deviation of measurement noise

max_itr = 30;      % maximum number of iterations
epsilon = 1e-5;    % convergence criteria
h_min   = 0;       % minimum horizontal value in plots
h_max   = 1;       % maximum horizontal value in plots
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
%
% number of measurements per time point
m     = 1;
%
% number of state vector components per time point
n     = 2;
%
% simulate the true trajectory and measurement noise
x1_true = av'; % velocity on top!
x2_true = ang'; % orientatino on bottom!
x_true  = [ x1_true ; x2_true ];
%

% measurement values and model
v_true  = sigma * randn(1, N);

bias = avb;
Pone    =  1; % measure bias directly

z       = x1_true + v_true + Pone*bias;

rk      = sigma * sigma;
rinvk   = 1 / rk;
rinv    = zeros(m, m, N);
h       = zeros(m, N);
dh      = zeros(m, n, N);
for k = 1 : N
    rinv(:, :, k) = rinvk;
    h(:, k)       = 0;
    dh(:,:, k)    = [ 1 , 0 ]; % measuring av
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

fprintf('Old model!!\n');
[xOutOld, uOut, info] = ...
    ckbs_affine(max_itr, epsilon, z, b, g, h, db, dg, dh, qinv, rinv);

fprintf('New model!!\n');
[xOutNew, yOut, yVar] = ...
    ckbs_augmented_affine(max_itr, epsilon, z, b, g, h, db, dg, dh, qinv, rinv, Pone);


%% now set up the 3-state smoother. 
clear qinvk;
n     = 3;
invDelta = 1e5;
g       = zeros(n, N);
dg      = zeros(n, n, N);
qinv    = zeros(n, n, N);
qk      = gamma * [ dt , dt^2/2 ; dt^2/2 , dt^3/3 ];
qinvk   = [inv(qk) zeros(2,1); 0 0 invDelta];
for k = 2 : N
    g(:, k)       = 0;
    dg(:,:, k)    = [ 1 , 0 0; dt , 1, 0; 0 0 1 ];
    qinv(:,:, k)  = qinvk;
end
%
% initial state estimate
g(:, 1)      = [x_true(:, 1); 10];
qinv(:,:, 1) = [100 0 0; 0 100 0; 0 0 1e-5];
%
% constraints


h       = zeros(m, N);
dh      = zeros(m, n, N);
for k = 1 : N
    rinv(:, :, k) = rinvk;
    h(:, k)       = 0;
    dh(:,:, k)    = [ 1 , 0, 1 ]; % measuring av
end


b       = zeros(0, N);

db      = zeros(0, n, N);

fprintf('Old Augmented model!!\n');
[xOutBiasState, ~, ~] = ...
    ckbs_affine(max_itr, epsilon, z, b, g, h, db, dg, dh, qinv, rinv);



fprintf('yOut is %5.3f\n', yOut);

if draw_plot
    figure(1);
    clf
    hold on
    plot(t', x_true(2,:)', 'k-', 'linewidth', 1.5 );
    plot(t', xOutOld(2,:)', 'b--', 'linewidth', 1.5 );
    plot(t', xOutNew(2,:)', 'r-.', 'linewidth', 1.5 );
    plot(t', xOutBiasState(2,:)', 'g--', 'linewidth', 1.5 );

        
    %     plot(t', - ones(N,1), 'b-');
    %      plot(t', ones(N,1), 'b-');
    axis([h_min, h_max, v_min, v_max]);
    title('Affine smoother');
 %   legend('Truth', 'Naive Estimate', 'Debiased Estimate', 'Location', 'SouthEast' );
    hold off
    axis tight;
    
    figure(2);
    clf
    hold on
    plot(t', x_true(1,:)', 'k-.', 'linewidth', 1.5 );
    plot(t', z(1,:)', 'ko' );
    plot(t', xOutOld(1,:)', 'b--', 'linewidth', 1.5 );
    plot(t', xOutNew(1,:)', 'r-.', 'linewidth', 1.5 );
    plot(t', xOutBiasState(1,:)', 'g--', 'linewidth', 1.5 );

        axis tight;
    %     plot(t', - ones(N,1), 'b-');
    %      plot(t', ones(N,1), 'b-');
%    axis([h_min, h_max, v_min, v_max]);
    title('Affine smoother');
    %legend('Truth', 'Measurements', 'Naive Estimate', 'Debiased Estimate', 'Location', 'SouthEast' );
    hold off
    
    figure(3)
    plot(t', xOutBiasState(3,:)', 'g--', 'linewidth', 1.5);
    hold on
    plot(t', yOut*ones(N,1), 'r-.', 'linewidth', 1.5);
        plot(t', yOut*ones(N,1) + yVar, 'b.', 'linewidth', 2.5);
    plot(t', yOut*ones(N,1) - yVar, 'b.', 'linewidth', 2.5);

    hold off;
    %
    % constrained estimate
    %        x_con = xOut;
end
 savefig(1:3, ['/Users/saravkin/Dropbox/KalmanCovariates/tex/figures/gyro']);
%

end
% $$ $newlinech %$$
% $end
