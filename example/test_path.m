% -------------------------------------------------------------------
% ckbs: Constrained Kalman-Bucy Smoother Program: Copyright (C) 2006
% Authors: Bradlely Bell:        bradbell at washington dot edu
%          Gianluigi Pillonetto: giapi at dei dot unipd dot it
% License: GNU General Public License Version 2
% -------------------------------------------------------------------
% $begin test_path.m$$ $newlinech %$$
% $spell
%       ind
%       findstr
%       addpath
%       src
% $$
%
% $section Set Up Path for Testing$$
%
% $index test_path$$
% $index path, for testing$$
% $index test, path$$
%
% $head Syntax$$
% $code test_path$$
%
% $head Purpose$$
% The directory $code ../src$$ is added to the current path setting.
% This is needed when running the tests in the $code test$$ subdirectory.
%
% $head Source Code$$
% $newlinech $$ $codep
function test_path()
    addpath('../src')
    addpath('nonlinear')
end
% $$ $newlinech %$$
% $end
