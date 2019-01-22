
function [sys, h, hh] = uy2ss_pk(u, y, n, delta, e, f, ep, fp)
if ~exist('e'),  e  = []; f  = []; end
if ~exist('ep'), ep = []; fp = []; end

p = size(y, 2); ell = ceil(n / p);
if ~exist('delta'), delta = 2 * ell + 1; end

hh  = uy2h_pk(u, y, ell, delta, e, f, ep, fp); 
sys = h2ss(hh, ell); 
h   = impulse(sys, delta - 1);
