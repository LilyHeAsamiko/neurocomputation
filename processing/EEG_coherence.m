%%EEG
% direct
figure,
plot(zef_EEG.sensors(:,1),zef_EEG.sensors(:,2),'o');

% index1 = find(zef_EEG.sensors(:,1)>25);
% index2 = find(zef_EEG.sensors(:,1)<50);
% index3 = find(zef_EEG.sensors(:,2)>-20);
% index4 = find(zef_EEG.sensors(:,2)>70);
% for i = 1:length(index2)
%      index5(i,:) = find(index1==index2(i));
% end
% 
% +for i = 1:length(index4) 
%      index6(i,:) = find(index3==index4(i));
% end
% 
% for i = 1:length(index6) 
%      index7(i,:) = find(index5==index6(i));
% end
% 
% for i = 1:length(index7) 
%     X(i,:) = EEG_T(index7(i),:);
% end
X(1,:) =  EEG_T(5,:);
X(2,:) =  EEG_T(33,:);
X(3,:) =  EEG_T(18,:);

fs = 1200;

%R(1,:) = coherence(S(1,:),S(2,:),3,5);
C(1,:) = mscohere(X(1,:),X(2,:),hamming(100),[],[],fs,'mimo');
%Ra1b1 = wcoherence(Xa1,Xb1);
%R(2,:) = correlation(S(2,:),S(3,:),3,5);
C(2,:) = mscohere(X(2,:),X(3,:),hamming(100),[],[],fs,'mimo');
%Ra1I1 = wcoherence(Xa1,XI1);
%R(3,:) = correlation(S(3,:),S(3,:),3,5);
C(3,:) = mscohere(X(3,:),X(1,:),hamming(100),[],[],fs,'mimo');
[h,p,ci,stats] = ttest(C(1,:),C(2,:))
%partial
%RI1b1 = wcoherence(XI1,Xb1);
% normalized_partial(remove1Ra1b1 = correlation(Fa1,Fb1,3,5);
Cpa1 = (C(1,:)-C(3,:).*C(2,:))./sqrt((1-C(3,:).^2).*(1-C(2,:).^2));
%Rpa1 = (R(1,:)-R(3,:).*R(2,:))./sqrt((1-R(3,:).^2).*(1-R(2,:).^2));

figure,
area = ['a','b','1','a'];
for i =1:3
    subplot(4,1,i)
    plot(C(i,:));
    title(['coherence between area ',area(i),'& ',area(i+1)])
end
subplot(4,1,4)
plot(Cpa1);
title(['partial correlation remove ',area(1)])


[H_mul,~] = MVAR(EEG_T,1,360);
figure,
for i = 1:4
    subplot(4,1,i)
    im1agesc(abs(H_mul{90*i}));
    title(['transform matrix for coherence analysis on time ',int2str(90*i)])
end
%pspectrogram((H_mul{360})^2,'interp');
R = H_mul;
k = length(R);
bics = zeros(1,k);
    for i =1:k
        n = size(R{i},1);
        m = mean(R{i},2);
        RSS = mean(sum(((R{i}-m).^2),2));
        bics(i) = n*log(RSS/n)*k-k*log(n);
    end
    plot(1:k,bics,'--');
    hold on    
    [bicR, s] = min(bics(1:k));    
    plot(s,bicR,'o')
    xlabel('order p')
    ylabel('BIC criteria')
%granger_prediction
figure,
for i = 1:16
subplot(4,4,i)
DTFX = DTF(EEG_T,i,j);
area(DTFX)
title(['D1TF',int2str(1),'to',int2str(i)])
end

for i = 1:4
    for j = 1:i
        PDCX{j,i} = PDC(EEG_T,j,i)
        figure,
        area(PDCX{j,i})
        title(['PDC of ', int2str(j),' to ',int2str(i)])
    end
end
GP = granger_prediction(EEG_T,4,150);
