%Objective function for the fixed point problem
function [f,g,h] = phiP(x,d)
%   f = x'*d.*(log(x)-1);
%   g = d.*log(x);
    g = d.*log(x);
    f = x'*(g-d);
    h = d./x;
end
