% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2011
% Authors: Bradlely Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
%          James Burke:          burke at math dot washington dot edu
%          Aleksandr Aravkin:    saravkin at eos dot ubc dot ca
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin vanderpol_sim$$ $newlinech %$$
% $spell
%       sim
%       Runge
%       mu
%       xi
%       xk
%       tk
%       dxk
%       nargin
%       vanderpol
%       Gillijns
%       Kandepu Foss Imsland
% $$
%
% $section Van der Pol Oscillator Simulation (No Noise)$$
% $index noiseless, vanderpol oscillator $$
% $index vanderpol, oscillator no noise$$
% $index oscillator, vanderpol no noise$$
%
% $head Syntax$$
% $codei%[%x_out%] = vanderpol_sim(%$$
% $icode%mu%, %xi%, %n_out%, %step%)%$$
%
% $head Differential Equation$$
% We use $latex x_1 (t)$$ and $latex x_2 (t)$$
% to denote the oscillator position and velocity as a function of time.
% The ordinary differential equation for the
% Van der Pol oscillator with no noise satisfies the differential equation
% $latex \[
% \begin{array}{rcl}
%       x_1 '(t) & = & x_2 (t)
%       \\
%       x_2 '(t) & = & \mu [ 1 - x_1(t)^2 ] x_2 (t) - x_1(t)
% \end{array}
% \] $$
%
% $head mu$$
% Is a scalar specifying the value of $latex \mu$$ is the
% differential equation above.
%
% $head xi$$
% is a column vector with two elements specifying the initial value for
% $latex x(t) \in \B{R}^2 $$.
% To be specific, $latex x_1 (0)$$ is equal to $icode%xi%(1)%$$ and
% $latex x_2(0)$$ is equal to $icode%xi%(2)%$$.
%
% $head n_out$$
% is an integer scalar that specifies the number of time points at which
% the approximate solution to the ODE is returned.
%
% $head step$$
% is a scalar that specifies the difference in time between points at which
% the solution of the ODE is approximated.
%
% $head x_out$$
% is a matrix with row size two and column size $icode%n_out%$$ that
% contains the approximation solution to the ODE.
% To be specific, for $icode k = 1 , ... , n_out$$,
% $icode%x_out%(i,k)%$$ is an approximation for
% $latex x_i [ (k-1) \Delta t ]$$.
%
% $head Method$$
% A fourth-order Runge method with step size $icode step$$ is used
% to approximate the solution of the ODE.
%
% $children%
%       example/nonlinear/vanderpol_sim_ok.m
% %$$
% $head Example$$
% The file $cref vanderpol_sim_ok.m$$ is an example and test
% of $code vanderpol_sim$$.
% It returns true, if the test passes, and false otherwise.
%
% $end
% ------------------------------------------------------------------------
function [x_out] = vanderpol_sim(mu, xi, n_out, step)
    % anonymous function for the ODE
    ode        = @(x) [ x(2) ; mu * ( 1 - x(1)*x(1) ) * x(2) - x(1) ];
    x_out      = zeros(2, n_out);
    x_out(:,1) = xi;
    step_2     = step / 2;
    step_6     = step / 6;
    for k = 2 : n_out
        x1 = x_out(:,k-1);
        f1 = ode(x1);
        x2 = x1 + step_2 * f1;
        f2 = ode(x2);
        x3 = x1 + step_2 * f2;
        f3 = ode(x3);
        x4 = x1 + step   * f3;
        f4 = ode(x4);
        x_out(:,k) = x1 + step_6 * (f1 + 2 * f2 + 2 * f3 + f4);
    end
    return
end
% $$ $newlinech %$$
% $end
