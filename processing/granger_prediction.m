function GC = granger_prediction(X,i,t)
Xn = X(i,1:t);
Xn_0 = X(i-1,1:t);
Vn = var(Xn,0,2)
Vn_0=var(Xn_0,0,2)
GC = log(Vn/Vn_0);
end