% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006-2013
% Authors: Aleksandr Y. Aravkin: saravkin at us dot ibm dot com
%          Bradley Bell:         bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin ckbs_newton_step_L1$$ $newlinech %$$
% $spell
%       ckbs
%       blkdiag
%       Hdia
%       Hlow
%       Bdia
%       tri
%       ds
%       dy
%       du
%       mu
%       Kuhn
%   tridiagonal
%     cccccc
%       sqrt
% $$
%
% $section Affine Robust L1 Kalman Bucy Smoother Newton Step$$
%
% $index ckbs_newton_step_L1$$
% $index newton_step_L1$$
%
% $index step, Newton$$
%
% $head Syntax$$
% $codei/[/dPPlus/, /dPMinus/, /dR/, /dS/, /dY/] = ckbs_newton_step_L1(/
%       mu/, /s/, /y/, /r/, /b/, /d/, /Bdia/, /Hdia/, /Hlow/,
%       /pPlus/, /pMinus/)/$$
%
% $head Purpose$$
% This routine computes one step of Newton's method for solving
% the non-linear Kuhn-Tucker conditions for
% the $latex \mu$$-relaxed affine robust L1 Kalman-Bucy smoother problem.
%
%
% $head Problem$$
% Given
% $latex \mu \in \B{R}_+$$,
% $latex s \in \B{R}^{m \times N}$$,
% $latex y \in \B{R}^{n \times N}$$,
% $latex r \in \B{R}^{m \times N}$$,
% $latex b \in \B{R}^{m \times N}$$,
% $latex d \in \B{R}^{n \times N}$$,
% $latex B \in \B{R}^{m \times n \times N}$$,
% $latex H \in \B{R}^{n \times n \times N}$$,
% $latex p^+ \in \B{R}^{m \times N}$$,
% $latex p^- \in \B{R}^{m \times N}$$,
% the $latex \mu$$-relaxed affine L1 robust Kalman-Bucy smoother problem is:
%$latex \[
%\begin{array}{ll}
%{\rm minimize}
%& \frac{1}{2} y^\R{T} H(x) y + d(x)^\R{T} y
%       + \sqrt{\B{2}}^\R{T} (p^+ + p^-) -\mu \sum_{i =
%       1}^{mN} \log(p_i^+) - \sum_{i=1}^{mN} \mu \log(p_i^-)
%    \;{\rm w.r.t} \;  y \in \B{R}^{nN}\; , \; p^+ \in \B{R}_+^{M} \; , \; p^- \in \B{R}_+^{M}
%\\
%{\rm subject \; to} &  b(x) + B(x) y - p^+ + p^- = 0
%\end{array}
%\] $$
% In addition, $latex H$$ is symmetric block tri-diagonal with each block of
% size $latex n \times n$$ and
% $latex B$$ is block diagonal with each block of size $latex m \times n$$
%
%
% $head Lagrangian$$
% We use $latex r, \; s \in \B{R}^{m \times N}$$
% to denote $latex \mu /p^+\;,\; \mu/p^-$$, respectively, and
% we denote by $latex q$$ the lagrange multiplier associated to the
% equality constraint. We also use  $latex \B{1}$$
% to denote the vector of length $latex mN$$ with all its components
% equal to one, and $latex \B{\sqrt{2}}$$ to denote the vector of
% length $latex mN$$ with all its components equal to $latex
% \sqrt{2}$$.
% The corresponding Lagrangian is
% $latex \[
% L(p^+, p^-, y, q)  =
% \frac{1}{2} y^\R{T} H y + d^\R{T} y + \B{\sqrt{2}}^T(p^+ + p^-)
% - \mu \sum_{i=1}^{mN} \log(p_i^+) - \mu\sum_{i=1}^{mN}\log(p_i^-)
% + q^\R{T} (b + B y - p^+ + p^-)
% \] $$
% The partial gradients of the Lagrangian are given by
% $latex \[
% \begin{array}{rcl}
% \nabla_p^+ L(p^+, p^-, y, q ) & = & \B{\sqrt{2}} - q - r \\
% \nabla_p^- L(p^+, p^-, y, q) & = & \B{\sqrt{2}} + q - s \\
% \nabla_y L(p^+, p^-, y, q ) & = & H y + c + B^\R{T} q \\
% \nabla_q L(p^+, p^-, y, q ) & = & b + B y - p^+ + p^- \\
% \end{array}
% \] $$
% From the first two of the above equations,
% we have $latex q = (r - s)/2$$.
%
% $head Kuhn-Tucker Conditions$$
% We use $latex D(s)$$
% to denote the diagonal matrix with $latex s$$
% along its diagonal.
% The Kuhn-Tucker Residual function
% $latex F : \B{R}^{4mN + nN} \rightarrow \B{R}^{4mN + nN}$$
% is defined by
% $latex \[
% F(p^+, p^-, r, s, y)
% =
% \left(
% \begin{array}{c}
% p^+ - p^- - b - B y       \\
% D(p^-) D(s) \B{1} - \tau \B{1} \\
% r + s - 2 \B{\sqrt{2}} \\
% D(p^+) D(r ) \B{1} - \tau \B{1}   \\
% H y + d + B^\R{T} (r - s)/2
% \end{array}
% \right)
% \] $$
% The Kuhn-Tucker conditions for a solution of the
% $latex \mu$$-relaxed constrained affine Kalman-Bucy smoother problem is
% $latex F(p^+, p^-, r, s, y) = 0 $$.
%
% $subhead Newton Step$$
% Given a value for $latex (p^+, p^-, r, s, y)$$, the Newton step
% $latex ( (\Delta p^+)^\R{T} , (\Delta p^-)^\R{T} , \Delta
% r^\R{T}, \Delta s^\R{T}, \Delta y^\R{T} )^\R{T}$$ solves the problem:
% $latex \[
% F_\mu^{(1)} (p^+, p^-, r, s, y)
% \left( \begin{array}{c} \Delta p^+ \\ \Delta p^- \\ \Delta r \\
% \Delta s \\ \Delta y \end{array} \right)
% =
% - F(p^+, p^-, r, s, y)
% \] $$
%
% $head mu$$
% The argument $icode mu$$ is a positive scalar specifying the
% relaxation parameter $latex \mu$$.
%
% $head s$$
% The argument $icode s$$ is an array of size $latex m \times
% N$$.
% All the elements of $icode s$$ are greater than zero.
%
% $head y$$
% The argument $icode y$$ is an array of size $latex n \times
% N$$
%
% $head r$$
% The argument $icode r$$ is an array of size $latex m \times
% N$$. All the elements of $icode r$$ are greater than zero.
%
% $head b$$
% The argument $icode b$$ is an array of size $latex m \times N$$.
%
% $head d$$
% The argument $icode d$$ is an array of size $latex n \times
% N$$.
%
% $head Bdia$$
% The argument $icode Bdia$$ is an $latex m \times n \times N$$ array.
% For $latex k = 1 , \ldots , N$$ we define
% $latex B_k \in \B{R}^{m \times n}$$ by
% $latex \[
%       B_k = Bdia(:, :, k)
% \] $$
%
% $head B$$
% The matrix $latex B$$ is defined by
% $latex \[
% B
% =
% \left( \begin{array}{cccc}
% B_1 & 0      & 0      &           \\
% 0   & B_2    & 0      & 0         \\
% 0   & 0      & \ddots & 0         \\
%     & 0      & 0      & B_N
% \end{array} \right)
% \] $$
%
% $head Hdia$$
% The argument $icode Hdia$$ is an $latex n \times n \times N$$ array.
% For $latex k = 1 , \ldots , N$$ we define
% $latex H_k \in \B{R}^{n \times n}$$ by
% $latex \[
%       H_k = Hdia(:, :, k)
% \] $$
%
%
% $head Hlow$$
% The argument $icode Hlow$$ is an $latex n \times n \times N$$ array.
% For $latex k = 1 , \ldots , N$$ we define
% $latex L_k \in \B{R}^{n \times n}$$ by
% $latex \[
%       L_k = Hlow(:, :, k)
% \] $$
%
% $head H$$
% The matrix $latex H$$ is defined by
% $latex \[
% H
% =
% \left( \begin{array}{cccc}
% H_1 & L_2^\R{T} & 0         &           \\
% L_2 & H_2       & L_3^\R{T} & 0         \\
% 0   & \ddots    & \ddots    & \ddots    \\
%     & 0         & L_N       & H_N
% \end{array} \right)
% \] $$
%
% $head dPPlus$$
% The result $icode dPPlus$$ is an array of size $latex m \times N$$
% equal to the $latex \Delta p^+$$ components of the Newton step.
%
% $head dPMinus$$
% The result $icode dPMinus$$ is an array of size $latex m \times N$$
% equal to the $latex \Delta p^-$$ components of the Newton step.
%
% $head dR$$
% The result $icode dR$$ is an array of size $latex m \times N$$
% equal to the $latex \Delta r$$ components of the Newton step.
%
% $head dS$$
% The result $icode dS$$ is an array of size $latex m \times N$$
% equal to the $latex \Delta s$$ components of the Newton step.
%
% $head dY$$
% The result $icode dY$$ is an array of size $latex n \times N$$
% equal to the $latex \Delta y$$ components of the Newton step.
%
%
% $children#
%       example/newton_step_L1_ok.m
% #$$
%
% $head Example$$
% The file $cref newton_step_L1_ok.m$$ contains an example and test of
% $code ckbs_newton_step_L1$$.
% It returns true if $code ckbs_newton_step_L1$$ passes the test
% and false otherwise.
%
% $head Method$$
% The derivative of $latex F$$ is given by
% $latex \[
% F_\mu^{(1)} (p^+, p^-, r, s, y) =
% \left(
% \begin{array}{ccccccc}
% I &  -I  & 0 & 0 & -B  \\
% 0 & D( s^- )   & 0   & D(p^-) & 0 \\
% 0 & 0 & I & I & 0 \\
% D( s^+) & 0 & D( p^+ ) & 0 & 0 \\
% 0 & 0 & - 0.5 B^\R{T} & 0.5 B^\R{T} & C
% \end{array}
% \right)
% \] $$

% Given the inputs
% $latex p^+ , p^-, s^+ ,  s^- , y , B , b , C $$ and $latex c $$,
% the following algorithm solves the Newton System for
% $latex \Delta p^+ , \Delta p^- , \Delta s^+ , \Delta s^- $$,
% and $latex \Delta y$$:
% \[
% \begin{array}{cccccc}
%\bar{d}  &= &  \tau \B{1} / s^+  - \tau \B{1} / s^- - b - B y + p^+
% \\
% \bar{e}  &= & B^\R{T} ( \sqrt{\B{2}} - s^- ) - C y - c
% \\
% \bar{f}  &= & \bar{d} - D( s^+ )^{-1} D( p^+ ) ( 2 \sqrt{\B{2}} - s^- )
% \\
% T        &= & D( s^+ )^{-1} D( p^+ ) + D( s^- )^{-1} D( p^- )
% \\
% \Delta y &= &
%       [ C + B^\R{T}  T^{-1} B ]^{-1} ( \bar{e} + B^\R{T} T^{-1} \bar{f} )
% \\
% \Delta s^- &= &
%       T^{-1} B \Delta y - T^{-1} \bar{f}
% \\
% \Delta s^+ &= &
%       - \Delta s^- +  2 \sqrt{\B{2}} - s^+ - s^-
% \\
% \Delta p^- &= &
%       D( s^- )^{-1} [  \tau \B{1} - D( p^- ) \Delta s^- ] - p^-
% \\
% \Delta p^+ &= &
%       \Delta p^- + B \Delta y + b + B y - p^+ + p^-
% \end{array}
% \]
% where the matrix $latex T$$ is diagonal with positive elements,
% and the matrix
% $latex C + B^\R{T} T^{-1} B \in \B{R}^{N n \times N n }$$ is
% block tridiagonal positive definite.
% $end
% ---------------------------------------------------------------------------



function [dPPlus, dPMinus, dR, dS, dY] = ...
    ckbs_newton_step_L1(mu, s, y, r, b, d, Bdia, Hdia, Hlow, pPlus, pMinus)


    % pPlus     m x N
    % pMinus    m x N
    % y         n x N
    % r         m x N
    % s         m x N
    % d         n x N
    % Hdia     n x n x N
    % Hlow  n x n x N
    % b         m x N
    % Bdia     m x n x N
    % mu        scalar


    %%
    % Determine the size of the problem
    n = size(y, 1);
    m = size(b, 1);
    N = size(y, 2);
    %

    %
    % check sizes
    if N ~= size(pPlus,2) || N ~= size(pMinus,2) || N ~= size(y,2) || ...
            N ~= size(r,2) || N ~= size(s,2) || N ~= size(d, 2) || ...
            N ~= size(Hdia, 3)|| N ~= size(Hlow, 3) || N ~= size(b, 2) ||...
            N ~= size(Bdia, 3)

        N

        size(pPlus,2)
        size(pMinus,2)
        size(y,2)
        size(r,2)
        size(s,2)
        size(d,2)
        size(Hdia,3)
        size(Hlow, 3)
        size(b, 2)
        size(Bdia, 3)


        error('L1_kalman_solver: argument sizes do not agree');
    end
    if n ~= size(y,1) || n~= size(d, 1) || n ~= size(Hdia,1) || ...
            n ~= size(Hdia,2) || n ~= size(Hlow,1) || n ~= size(Hlow,2) || ...
            n ~= size(Bdia,2)
        n

        size(y,1)
        size(d,1)
        size(Hdia,1)
        size(Hdia,2)
        size(Hlow,1)
        size(Hlow,2)
        size(Bdia,2)

        error('L1_kalman_solver: argument sizes do not agree');
    end
    if m ~= size(pPlus,1) || m ~= size(pMinus,1) || m ~= size(r,1) || m ~= size(s,1)||...
            m ~= size(b, 1) || m ~= size(Bdia, 1)
        m

        size(pPlus,1)
        size(pMinus,1)
        size(r,1)
        size(s,1)
        size(b, 1)
        size(Bdia, 1)

        error('L1_kalman_solver: argument sizes do not agree');
    end

    %%
    % Create required datastructures


    if min(r) <= 0
        r
        error('L1_kalman_solver: min(r) <= 0 ');
    end
    if min(s) <= 0
        s
        error('L1_kalman_solver: min(s) <= 0 ');
    end
    if min(pPlus) <= 0
        pPlus
        error('L1_kalman_solver: min(pPlus) <= 0 ');
    end
    if min(pMinus) <= 0
        pMinus
        error('L1_kalman_solver: min(pMinus) <= 0 ');
    end


    rs = r.*s;
    rpm = r.*pMinus;
    spp = s.*pPlus;
    tinv = rs./(rpm + spp);
    tinvir = s./(rpm + spp);
    tinvis = r./(rpm + spp);

    % Create a new array for modified diagonal blocks
    modCDiag = Hdia;

    %Form B'T^{-1}B and add it to the diagonal C terms
    for iter = 1:N
        modCDiag(:,:,iter) = modCDiag(:,:,iter) +...
            Bdia(:,:,iter)'*diag(tinv(:,iter))*Bdia(:,:,iter);
    end


    % Compute D(p^+) and D(p^-)
    Dpp = pPlus(:);
    Dpm = pMinus(:);

    Cy = ckbs_blktridiag_mul(Hdia, Hlow, y(:));
    By = ckbs_blkdiag_mul(Bdia, y(:));

    temp = s(:) - sqrt(2) + tinv(:).*(b(:) + By - Dpp) + mu*tinvis(:) + ...
        tinvir(:).*(2*sqrt(2)*Dpp - mu - Dpp.*s(:));
    BTtemp = ckbs_blkdiag_mul_t(Bdia, temp);

    R5 = Cy + d(:) + BTtemp;

    dY = -ckbs_tridiag_solve_b(modCDiag, Hlow, R5);


    %Bymdy = ckbs_blkdiag_mul(Bdia, y(:) - dY);
    Bypdy = ckbs_blkdiag_mul(Bdia, y(:) + dY);

    dS = tinv(:).*(b(:) + Bypdy - Dpp) + tinvis(:)*mu + tinvir(:).*(2*sqrt(2)*Dpp - mu - Dpp.*s(:));

    dR = 2*sqrt(2) - r(:) - s(:) - dS;

    dPMinus = (mu - Dpm.*dS)./s(:) - Dpm;

    dPPlus = dPMinus + Bypdy + b(:) + Dpm - Dpp;




    %%

    %%
    % Modify output: back into multi-dim arrays.
    %%
    dPPlus = reshape(dPPlus, m, N);
    dPMinus = reshape(dPMinus, m, N);
    dY = reshape(dY, n, N);
    dR = reshape(dR, m, N);
    dS = reshape(dS, m, N);
    %dS = zeros(m, N);
    %%
    % Output error check:

    if N ~= size(dPPlus, 2) || N ~= size(dPMinus, 2) || N ~= size(dY, 2)|| ...
            N ~= size(dR, 2) || N ~= size(dS, 2)

        N

        size(dPPlus, 2)
        size(dPMinus, 2)
        size(dY, 2)
        size(dR, 2)
        size(dS, 2)

        error('L1_kalman_solver: argument sizes do not agree');
    end

    if n ~= size(dY, 1)
        n
        size(dY, 1)

        error('L1_kalman_solver: argument sizes do not agree');
    end

    if m ~= size(dPPlus, 1) || m ~= size(dPMinus, 1) ||...
            m ~= size(dR, 1) || m ~= size(dS, 1)
        m
        size(dPPlus, 1)
        size(dPMinus, 1)
        size(dR, 1)
        size(dS, 1)


        error('L1_kalman_solver: argument sizes do not agree');
    end

    %%


    return
end

