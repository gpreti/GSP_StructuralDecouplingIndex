%PERCENTILE: The kth percentile Pk is that value of X, say Xk, which
%  corresponds to a cumulative frequency of Nk/100.   ( cited from 
%  Eric W. Weisstein. "Percentile." From MathWorld--A Wolfram Web
%  Resource. http://mathworld.wolfram.com/Percentile.html )
%
%  Usage: Pk = percentile(X, Nk);
%
%  In which, X can be an 1D vector or a 2D matrix. and Nk can be a
%  scalar or a vector.
%
%  If Nk is a scalar and X is an 1D vector, Pk will be the percentile
%  of this 1D vector. If Nk is a scalar and X is a 2D matrix, Pk will
%  be a row vector. Each will represent the percentile of each col in
%  this 2D matrix. If Nk is a vector, the result Pk will be stacked
%  together.
%
%  Nathasha fixed a bug when 2D matrix is not a square one,
%  see notes below
%
function Pk = percentile(X, Nk)

   X = squeeze(X);

   if length(Nk)~=length(Nk(:))				% wrong Nk
      error('Nk must either be a scalar or a vector');
   else
      if (length(X)==length(X(:))) & (size(X,1)==1)	% X is 1D row
         X = X';
      end
   end

   [r c] = size(X);					% NK fix

   x = [0, ([0.5 : (r-0.5)] ./ r) * 100, 100];
   y = [min(X); sort(X); max(X)];

   Pk = interp1(x, y, Nk);

   return;						% percentile

