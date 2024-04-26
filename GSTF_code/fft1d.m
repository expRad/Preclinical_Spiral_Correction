function sig=fft1d(img,dim);

% 1d fourier transformation along dimension dim
% by Felix Breuer

sig=fftshift(fft(fftshift(img,dim),[],dim),dim);