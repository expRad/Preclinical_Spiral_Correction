function sig=fft2d(img);

% 2d fourier transformation along last 2 dimensions 
% by Felix Breuer

dim=size(img);

shift=zeros(1,size(dim,2));
shift(1,end-1)=dim(end-1)/2;
shift(1,end)=dim(end)/2;



sig=cmshiftnd(fft(fft(cmshiftnd(img,shift),[],size(dim,2)-1),[],size(dim,2)),shift);