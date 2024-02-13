function Y = addnoise(X , snr)
sigma = sqrt(norm(X , 'fro').^2 / snr);
sz = size(X);
w = wgn(sz(1),sz(2) , 0);
w = w / norm(w , 'fro');
Y = X + sigma * w;
