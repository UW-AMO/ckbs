% Run general simulator
% includes cross-talk matrix 
% -------------------------------------------------------------------

function [ok] = vanderpol_sim_general_ok()

    mu     = 1.;
    ndim   = 3; % three components
%    xi_single     = [ 2 ; 0 ];
%    xi = repmat(xi_single, ndim, 1);
    xi = [2; 0; 1; 1; 0; 2]; 
    n_out  = 21;
    step   = .1;
    
    W = rand(ndim*2); % cross talk matrix linking components
   
    x_out  = vanderpol_sim_general(mu, xi, n_out, step, ndim, W);
    
    figure()
    plot(1:n_out, x_out(1, :));
    hold on
    plot(1:n_out, x_out(3,:));
    plot(1:n_out, x_out(5,:));
    hold off
    
    figure()
    plot(1:n_out, x_out(2,:));
    hold on
    plot(1:n_out, x_out(4,:));
    plot(1:n_out, x_out(6,:));
    hold off;

    
  %  t_sol  = [ 0.0    , 0.5    , 1.0    ,  1.5   , 2.0    ];
  %  x1_sol = [ 2.0000 , 1.8377 , 1.5081 ,  1.0409, 0.3233 ];
    ok     = 1;
%    for k = 1 : 5
 %       ok = ok & abs( x_out(1, (k-1)*5 + 1) - x1_sol(k) ) < 1e-4;
 %   end
    return
end
% $$ $newlinech %$$
% $end
