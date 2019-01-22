clear all, warning off

T = 50; delta = 10; N = 100; s = 0.01; eiv = 1;
nu = 1; ny = 1; n = 2;

m = 1; k = 10; d = 1;
A = [0 1; -k/m -d/m]; B = [0; 1/m]; C = [1 0]; D = 0;
sys = c2d(ss(A, B, C, D), 0.2); sys.ts = -1;
H = impulse(sys, delta - 1);

u0 = rand(T, nu); y0 = lsim(sys, u0); w0 = [u0 y0]; 
e = ones(1, delta); f = e * H; ep = []; fp = [];
disp('Running Monte-Carlo simulation ...')
for j = 1:N
  if eiv
    wt = randn(T, nu + ny); w = w0 + s * norm(w0) * wt / norm(wt); 
  else
    wt = [zeros(T, nu) randn(T, ny)]; w = w0 + s * norm(w0(:, nu + 1:end)) * wt / norm(wt); 
  end
  u = w(:, 1:nu); y = w(:, nu+1:end);
  [sysh1, Hh1, Hh1_] = uy2ss_pk(u, y, n, delta, e, f, ep, fp);
  [sysh2, Hh2, Hh2_] = uy2ss_pk(u, y, n, delta);
  sysh3 = n4sid(iddata(y, u), n); Hh3 = impulse(sysh3, delta - 1);
  esys(j, :) = [norm(sys - sysh1) norm(sys - sysh2) norm(sys - sysh3)] / norm(sys);
  eh(j, :)   = [norm(vec(H - Hh1)) norm(vec(H - Hh2)) norm(vec(H - Hh3))] / norm(vec(H));
  ehh(j, :)  = [norm(vec(H - Hh1_)) norm(vec(H - Hh2_)) NaN] / norm(vec(H));
  ehp(j, :)   = [norm(vec(e * Hh1 - f)) norm(vec(e * Hh2 - f)) norm(vec(e * Hh3 - f))];
  ehhp(j, :)  = [norm(vec(e * Hh1_ - f)) norm(vec(e * Hh2_ - f)) NaN];
end

[{'' 'uy2ss_pk' 'uy2ss' 'n4sid'}
['esys' num2cell(mean(esys))]
['eh'   num2cell(mean(eh))]
['ehh'  num2cell(mean(ehh))]
['ehp'   num2cell(mean(ehp))]
['ehhp'  num2cell(mean(ehhp))]]