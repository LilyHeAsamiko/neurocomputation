function [E] = effective_connectivity(sig1,sig2,s,t)
for s =1:S
    for t = 1:N-s-T+1
        R(s) = mean(sig1(t:t+T-1).*sig2(t+s:t+s+T-1));
    end
    for t = N-s-T+2:N-T+1
        R(s) = mean(sig1(t:t+T-1).*sig2(1:T));
    end
end
end