% UY2HBLK - From data to impulse response. 
%           (Block algorithm)
%
% [h, delta] = uy2hblk(u, y, lmax, delta)
%
% U,Y - input (TxM matrix) and output (TxP matrix)
% LMAX  - upper bound for the system lag
% DELTA - desired length of the impulse response
%  default: the maximal possible number is used
% H - first DELTA samples of the impulse response
%     H = [H0; H1; ...; Hdelta]

function [h, delta] = uy2hblk(u,y,lmax,delta)

% Constants
[T,m] = size(u);
p = size(y,2);

if ( nargin < 4 | isempty(delta) )
  delta = ceil((T+1)/(m+1) - lmax * (1 + p));
  if delta < 1
    error(sprintf('Not enough data for DELTA=%d.', ...
                  delta));
  end
end

L = lmax+delta; % # of block rows of the Hankel matrix
nrows = L*m+lmax*p; % # of rows in the system (7.4)

% QR of the Hankel matrix of the data (LHS of (7.5))
r = triu(qr([blkhank(u,L); blkhank(y,L)]'))';

% Select the R11 and R21 blocks of R
r11 = r(1:nrows,1:nrows);
r21 = r(nrows+1:end,1:nrows);

% Right-hand side of the system (7.4)
B = zeros(nrows,m); 
B(lmax*m+1:(lmax+1)*m,:) = eye(m);

% Solve the system (7.4) for G and compute H = Yf * G
h = r21 * pinv(r11) * B;

%  for MIMO (M-inputs, P-outputs) system, h is PxMxT
% if ( m > 1 | p > 1 ) % MIMO
%   % Construct a 3d output array
%   h_ = h;
%   h = zeros(p,m,delta);
%   for i = 1:delta
%     h(:,:,i) = h_((i-1)*p+1:i*p,:);
%   end
% else % SISO
%   h = h(1:delta,:);
% end