function [weight,cdfcheck] = betaweight(zmax,dz,a,b)
%
% computes weights according to the beta distribution
% for uniform distribution of dz
%
z = ((0:dz:zmax)/zmax)';
z(1) = .5/zmax;
coef = gamma(a+b)/(gamma(a)*gamma(b));
weight = coef*(z.^(a-1)).*((1-z).^(b-1));
cdfcheck = sum(weight.*dz/zmax);

