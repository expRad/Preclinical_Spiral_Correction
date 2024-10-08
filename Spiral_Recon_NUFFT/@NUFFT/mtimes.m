% This script is from the SparseMRI V2.0 collection by Michael Lustig, which can be
% downloaded here: https://people.eecs.berkeley.edu/~mlustig/Software.html

function ress = mtimes(a,bb)
% performs the normal nufft

for n=1:size(bb,3)
b = bb(:,:,n);

if a.adjoint
	b = b(:).*a.w(:);
	res = nufft_adj(b, a.st)/sqrt(prod(a.imSize));
	res = reshape(res, a.imSize(1), a.imSize(2));
	res = res.*conj(a.phase);
	if a.mode==1
		res = real(res);
	end

else
	b = reshape(b,a.imSize(1),a.imSize(2));
	if a.mode==1
        b = real(b);
	end
	b = b.*a.phase;
	res = nufft(b, a.st)/sqrt(prod(a.imSize));
	res = reshape(res,size(a.w,1),size(a.w,2));

end
ress(:,:,n) = res;
end

