clear all, warning off

T = 50; delta = 10; N = 100; s = 0.01; eiv = 1;
nu = 1; ny = 1; n = 2;

m = 1; k = 10; d = 1;
A = [0 1; -k/m -d/m]; B = [0; 1/m]; C = [1 0]; D = 0;
sys = c2d(ss(A, B, C, D), 0.2); sys.ts = -1;
H = impulse(sys, delta - 1);

u0 = rand(T, nu); y0 = lsim(sys, u0); w0 = [u0 y0]; 
E = rand(10, delta); F = E * H; ep = []; fp = [];
disp('Running Monte-Carlo simulation ...')
for i = 1:10
  e = E(1:i, :); f = F(1:i, :);
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
  Esys(i, :) = mean(esys); Eh(i, :) = mean(eh); Ehh(i, :) = mean(ehh);
  Ehp(i, :) = mean(ehp); Ehhp(i, :) = mean(ehhp); 
end

figure, hold on, 
plot(Esys(:, 1), 'b.', 'MarkerSize', 30), h(1) = plot(Esys(:, 1), 'b--', 'linewidth', 2); 
plot(Esys(:, 2), 'r.', 'MarkerSize', 30), h(2) = plot(Esys(:, 2), 'r-', 'linewidth', 2); 
plot(Esys(:, 3), 'r.', 'MarkerSize', 30), h(3) = plot(Esys(:, 3), 'r-.', 'linewidth', 2);
title(''), xlabel('number of equality constraints'), ylabel('esys')
legend(h, {'uy2ss\_pk', 'uy2ss', 'n4sid'}, 'location', 'SouthWest'), 
ax = axis; set(gca, 'xTick', 1:10), set(gca, 'yTick', [0 ax(4)])
ax = axis; axis([1 i ax(3:4)])
set(gca, 'fontsize', 15)
print -dpng '../results/f1a'

figure, hold on, 
plot(Ehh(:, 1), 'b.', 'MarkerSize', 30), h(1) = plot(Ehh(:, 1), 'b--', 'linewidth', 2); 
plot(Ehh(:, 2), 'r.', 'MarkerSize', 30), h(2) = plot(Ehh(:, 2), 'r-', 'linewidth', 2);  
title(''), xlabel('number of equality constraints'), ylabel('ehh')
legend(h(1:2), {'uy2ss\_pk', 'uy2ss'}, 'location', 'SouthWest'), 
ax = axis; set(gca, 'xTick', 1:10), set(gca, 'yTick', [0 ax(4)])
ax = axis; axis([1 i ax(3:4)])
set(gca, 'fontsize', 15)
print -dpng '../results/f1b'

figure, hold on, 
plot(Ehhp(:, 1), 'b.', 'MarkerSize', 30), h(1) = plot(Ehhp(:, 1), 'b--', 'linewidth', 2); 
plot(Ehhp(:, 2), 'r.', 'MarkerSize', 30), h(2) = plot(Ehhp(:, 2), 'r-', 'linewidth', 2);  
title(''), xlabel('number of equality constraints'), ylabel('ehhp')
legend(h(1:2), {'uy2ss\_pk', 'uy2ss'}, 'location', 'NorthWest'), 
ax = axis; set(gca, 'xTick', 1:10), set(gca, 'yTick', [0 ax(4)])
ax = axis; axis([1 i ax(3:4)])
set(gca, 'fontsize', 15)
print -dpng '../results/f2a'

figure, hold on, 
plot(Ehp(:, 1), 'b.', 'MarkerSize', 30), h(1) = plot(Ehp(:, 1), 'b--', 'linewidth', 2);
plot(Ehp(:, 2), 'r.', 'MarkerSize', 30), h(2) = plot(Ehp(:, 2), 'r-', 'linewidth', 2); 
plot(Ehp(:, 3), 'r.', 'MarkerSize', 30), h(3) = plot(Ehp(:, 3), 'r-.', 'linewidth', 2);   
title(''), xlabel('number of equality constraints'), ylabel('ehp')
legend(h, {'uy2ss\_pk', 'uy2ss', 'n4sid'}, 'location', 'NorthWest'), 
ax = axis; set(gca, 'xTick', 1:10), set(gca, 'yTick', [0 ax(4)])
ax = axis; axis([1 i ax(3:4)])
set(gca, 'fontsize', 15)
print -dpng '../results/f2b'