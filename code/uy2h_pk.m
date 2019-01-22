
function h = uy2h_pk(u, y, ell, delta, e, f, ep, fp)
[T, m] = size(u); [T, p] = size(y);
if (~exist('e') || isempty(e)) && (~exist('ep') || isempty(ep)) % no prior
  h = uy2hblk(u, y, ell, delta); return % Algorithm 1
end

Y  = blkhank(y, ell + delta); 
Yp = Y(1:ell * p, :); 
Yf = Y(ell * p + 1:end, :);

a = [blkhank(u, ell + delta); Yp];
b = zeros(size(a, 1), m); b(ell * m + 1:(ell + 1) * m, :) = eye(m);

% analytical solution
eYf = e * Yf; Gp = pinv(eYf) * f; N = null(eYf);
g = Gp + N * pinv(a * N) * (b - a * Gp);
% quadratic programming 
if (exist('ep') && ~isempty(ep))    
  g = lsqlin(a, b, ep * Yf, fp, eYf, f, [], [], g);
end
h = Yf * g;
