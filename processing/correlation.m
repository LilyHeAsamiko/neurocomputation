function [H] = MVAR(X)
% S = 5;
% T = 3;
% sig1 = X(1,:);
% sig2 = X(2,:);
% sig3 = X(3,:);
for i = 1:size(X,2)
A = corr(X(:,i),X(:,i));
A_H = fft(A);
E = A*X;
H{i} = inv(A);
end
end




