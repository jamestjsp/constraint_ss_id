clear all, warning off

T = 100; delta = 15; s = 0.04; eiv = 1;
nu = 1; ny = 1; n = 2;

m = 1; k = 10; d = 1;
A = [0 1; -k/m -d/m]; B = [0; 1/m]; C = [1 0]; D = 0;
sys = c2d(ss(A, B, C, D), 0.2); sys.ts = -1;
H = impulse(sys, delta - 1);
u0 = rand(T, nu); y0 = lsim(sys, u0); w0 = [u0 y0]; 

ep = [eye(delta); - eye(delta)]; fp = ep * H + 0.005;
e  = zeros(0, delta); f = zeros(0, 1);

  if eiv
    wt = randn(T, nu + ny); w = w0 + s * norm(w0) * wt / norm(wt); 
  else
    wt = [zeros(T, nu) randn(T, ny)]; w = w0 + s * norm(w0(:, nu + 1:end)) * wt / norm(wt); 
  end
  u = w(:, 1:nu); y = w(:, nu+1:end);
  [sysh1, Hh1, Hh1_] = uy2ss_pk(u, y, n, delta, e, f, ep, fp);
  [sysh2, Hh2, Hh2_] = uy2ss_pk(u, y, n, delta);
  sysh3 = n4sid(iddata(y, u), n); Hh3 = impulse(sysh3, delta - 1);
  esys = [norm(sys - sysh1) norm(sys - sysh2) norm(sys - sysh3)] / norm(sys);
  eh   = [norm(vec(H - Hh1)) norm(vec(H - Hh2)) norm(vec(H - Hh3))] / norm(vec(H));
  ehh  = [norm(vec(H - Hh1_)) norm(vec(H - Hh2_)) NaN] / norm(vec(H));

figure, hold on
plot(H, 'k', 'linewidth', 2)
plot(Hh1_, 'b-.', 'linewidth', 2)
plot(Hh2_, 'r-.', 'linewidth', 2)
plot(fp(1:delta), ':k', 'linewidth', 2)
plot(- fp(delta + 1:end), ':k')
ax = axis; axis([1 delta ax(3:4)])
title('estimated impulse response'), xlabel('time'), ylabel('h')
legend('true impulse response', 'uy2h\_pk', 'uy2h', 'location', 'NorthEast'), 
set(gca, 'fontsize', 15)
print -dpng '../results/f3a'

figure, hold on
plot(H, 'k', 'linewidth', 2)
plot(Hh1, 'b--', 'linewidth', 2), 
plot(Hh2, 'r-', 'linewidth', 2),  
plot(Hh3, 'r-.', 'linewidth', 2)
ax = axis; axis([1 delta ax(3:4)])
title('impulse response of sysh'), xlabel('time'), ylabel('h')
legend('true impulse response', 'uy2ss\_pk', 'uy2ss', 'n4sid', 'location', 'NorthEast'), 
set(gca, 'fontsize', 15)
print -dpng '../results/f3b'

[{'' 'uh2ss_pk' 'uh2ss' 'n4sid'}
['esys' num2cell(esys)]
['ehh'  num2cell(ehh) ]]
