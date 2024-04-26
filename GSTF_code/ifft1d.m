function out=ifft1d(in,dim);

out=ifftshift(ifft(ifftshift(in,dim),[],dim),dim);