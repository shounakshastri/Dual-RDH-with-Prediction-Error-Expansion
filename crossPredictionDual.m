function [Ic, er, p] = crossPredictionDual(I, o)

% o indicates whether the embedding is to be done in even rows or odd rows.
% if o == 0, embedding will be done in even rows. If o == 1, embedding will
% be done in odd rows.

[m, n] = size(I);
I1 = padarray(I, [1,1]);
p=zeros(size(I1));
Ic = I1;
for ii = 2:m+1
    for jj = 2:n+1
        if mod(ii+jj, 2) == 0
            if mod(ii, 2) == o
                Ic(ii, jj) = floor((I1(ii+1, jj) + I1(ii-1, jj) + I1(ii, jj+1) + I1(ii, jj-1))/4);
                p(ii,jj)=1;
            end
        end
    end
end
Ic = Ic(2:end-1, 2:end-1);
p = p(2:end-1, 2:end-1);
er = I - Ic;
end