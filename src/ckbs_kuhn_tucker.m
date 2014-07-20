% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006
% Authors: Bradlely Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin ckbs_kuhn_tucker$$ $newlinech %$$
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
% $$
%
% $section Compute Residual in Kuhn-Tucker Conditions$$
%
% $index ckbs_kuhn_tucker$$
% $index kuhn_tucker$$
%
% $index Kuhn-Tucker, residual$$
% $index Tucker, Kuhn residual$$
% $index residual, Kuhn-Tucker$$
%
% $head Syntax$$
% $codei/[/F/] = ckbs_kuhn_tucker(/
%       mu/, /s/, /y/, /u/, /b/, /d/, /Bdia/, /Hdia/, /Hlow/)/$$
%
% $head Purpose$$
% This routine computes the residual in the Kuhn-Tucker conditions for
% the $latex \mu$$-relaxed affine constrained Kalman-Bucy smoother problem.
%
% $head Problem$$
% Given
% $latex \mu \in \B{R}_+$$,
% $latex H \in \B{R}^{p \times p}$$,
% $latex d \in \B{R}^p$$,
% $latex b \in \B{R}^r$$, and
% $latex B \in \B{R}^{r \times p}$$,
% the $latex \mu$$-relaxed affine constrained Kalman-Bucy smoother problem is:
% $latex \[
% \begin{array}{rl}
% {\rm minimize} & \frac{1}{2} y^\R{T} H y + d^\R{T} y
% - \mu \sum_{i=1}^r \log(s_i)
% \; {\rm w.r.t} \; y \in \B{R}^p \; , \; s \in \B{R}_+^r
% \\
% {\rm subject \; to} & s + b + B y  = 0
% \end{array}
% \] $$
% In addition, $latex H$$ is symmetric block tri-diagonal with each block of
% size $latex n \times n$$ and
% $latex B$$ is block diagonal with each block of size $latex m \times n$$
% (there is an integer $latex N$$ such that
% $latex p = n * N$$ and $latex r = m * N$$).
%
%
% $head Lagrangian$$
% We use $latex u \in \B{R}^r$$
% to denote the Lagrange multipliers corresponding to the constraint equation.
% The corresponding Lagrangian is
% $latex \[
% L(y, s, u)  =
% \frac{1}{2} y^\R{T} H y + d^\R{T} y
% - \mu \sum_{i=1}^r \log(s_i)
% + u^\R{T} (s + b + B y)
% \] $$
% The partial gradients of the Lagrangian are given by
% $latex \[
% \begin{array}{rcl}
% \nabla_y L(y, s, u ) & = & H y + B^\R{T} u + d  \\
% \nabla_s L(y, s, u ) & = & u - \mu / s \\
% \nabla_u L(y, s, u ) & = & s + b + B y \\
% \end{array}
% \] $$
% where $latex \mu / s $$ is the component by component division of
% $latex \mu $$ by the components of the $latex s$$.
% Note, from the second equation, that we only need consider
% $latex u \geq 0$$ because $latex s \geq 0$$.
%
% $head Kuhn-Tucker Residual$$
% We use $latex D(s)$$
% to denote the diagonal matrix with $latex s$$
% along its diagonal and $latex 1_r$$
% to denote the vector, of length $latex r$$ with all its components
% equal to one.
% The Kuhn-Tucker Residual function
% $latex F : \B{R}^{r + p + r} \rightarrow \B{R}^{r + p + r}$$
% is defined by
% $latex \[
% F(s, y, u)
% =
% \left(
% \begin{array}{c}
% s + b + B y       \\
% H y + B^\R{T} u + d   \\
% D(s) D(u) 1_r - \mu 1_r
% \end{array}
% \right)
% \] $$
% The Kuhn-Tucker conditions for a solution of the
% $latex \mu$$-relaxed constrained affine Kalman-Bucy smoother problem is
% $latex F(s, y, u) = 0 $$.
%
% $head mu$$
% The argument $icode mu$$ is a positive scalar specifying the
% relaxation parameter $latex \mu$$.
%
% $head s$$
% The argument $icode s$$ is a column vector of length $latex r$$.
% All the elements of $icode s$$ are greater than zero.
%
% $head y$$
% The argument $icode y$$ is a column vector of length $latex p$$
%
% $head u$$
% The argument $icode u$$ is a column vector of length $latex r$$.
% All the elements of $icode s$$ are greater than zero.
%
% $head b$$
% The argument $icode b$$ is a column vector of length $latex r$$.
%
% $head d$$
% The argument $icode d$$ is a column vector of length $latex p$$
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
% $head F$$
% The result $icode F$$ is a column vector of length $latex r + p + r$$
% containing the value of the
% $cref/Kuhn-Tucker residual/ckbs_kuhn_tucker/Kuhn-Tucker Residual/$$; i.e.,
% $latex F(s, y, u)$$.
%
% $children#
%       example/kuhn_tucker_ok.m
% #$$
%
% $head Example$$
% The file $cref kuhn_tucker_ok.m$$ contains an example and test of
% $code ckbs_kuhn_tucker$$.
% It returns true if $code ckbs_kuhn_tucker$$ passes the test
% and false otherwise.
%
% $end
% ---------------------------------------------------------------------------
function [F] = ckbs_kuhn_tucker(mu, s, y, u, b, d, Bdia, Hdia, Hlow)
    m = size(Bdia, 1);
    n = size(Bdia, 2);
    N = size(Bdia, 3);
    p = n * N;
    r = m * N;
    %
    % check sizes
    if size(s,1) ~= r | size(u,1) ~= r | size(b,1) ~= r
        r
        size(s,1)
        size(u,1)
        size(b,1)
        error('ckbs_kuhn_tucker: argument sizes do not agree');
    end
    if size(y,1) ~= p | size(d,1) ~= p
        p
        size(y,1)
        size(d,1)
        error('ckbs_kuhn_tucker: argument sizes do not agree');
    end
    if size(Hdia,1) ~= n | size(Hdia,2) ~= n | ...
            size(Hlow,1) ~= n | size(Hlow,2) ~= n
        n
        size(Hdia,1)
        size(Hdia,2)
        size(Hlow,1)
        size(Hlow,2)
        error('ckbs_kuhn_tucker: argument sizes do not agree');
    end
    if size(Hdia,3) ~= N | size(Hlow,3) ~= N
        N
        size(Hdia,3)
        size(Hlow,3)
        error('ckbs_kuhn_tucker: argument sizes do not agree');
    end
    %
    % compute the necessary matrix products; i.e.,
    % B * y, H * y, and u^T * B
    blk_m = 1 : m;
    blk_n = 1 : n;
    B_y   = zeros(r, 1);
    H_y   = zeros(p, 1);
    Bt_u  = zeros(p, 1);

    B_y = ckbs_blkdiag_mul(Bdia, y);
    H_y = ckbs_blktridiag_mul(Hdia, Hlow, y);
    Bt_u = ckbs_blkdiag_mul_t(Bdia, u);


    %
    F = [ ...
        s + b + B_y ; ...
        H_y + Bt_u + d; ...
        s .* u - mu ...
        ];
    return
end
