function re = randexp(alpha,row,col);
%
re = -alpha*log(1-rand(row,col));
