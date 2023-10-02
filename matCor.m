function [r] = matCor(X, Y)


r = trace(X'*Y) / sqrt(trace(X'*X) * trace(Y'*Y));



end