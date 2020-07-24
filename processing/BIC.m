function [bicR,s] = BIC(R)
%R = H_mul;
k = length(R);
bics = zeros(1,k);
    for i =1:k
        n = size(R{i},1);
        m = mean(R{i},2);
        RSS = mean(sum(((R{i}-m).^2),2));
        bics(i) = n*log(RSS/n)*k-k*log(n);
    end
    [bicR, s] = min(bics(:));  
    figure,
    plot(1:k,bics,'--');
    hold on     
    plot(s,bicR,'o')
    xlabel('order p')
    ylabel('BIC criteria')
    title(['best order ',int2str(s),'found by BIC'])
end