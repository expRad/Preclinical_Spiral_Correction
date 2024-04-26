function img=ifft2d(sig)

dim=size(sig);

shift=zeros(1,size(dim,2));
shift(1,end-1)=dim(end-1)/2;
shift(1,end)=dim(end)/2;
%shift


img=icmshiftnd(ifft(ifft(icmshiftnd(sig,shift),[],size(dim,2)-1),[],size(dim,2)),shift); %% changed from cmshiftnd to icmshiftnd %% DLSTAEB
