% This script is from the SparseMRI V2.0 collection by Michael Lustig, which can be
% downloaded here: https://people.eecs.berkeley.edu/~mlustig/Software.html

function res = ctranspose(a)
a.adjoint = xor(a.adjoint,1);
res = a;

