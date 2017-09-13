function [S1]=normalization(S)

T=sqrt(sum(S,2)*sum(S,2)');
T(T==0)=eps;
S1 = S./T;

%idx2 = find(sum(S1) > 0);
%S2 = S1;
%S2(:) = 0;
%for ii = 1 : length(idx2)
%     if sum(S1(:,idx2(ii)))~=0
%    S2(:,idx2(ii)) = S1(:,idx2(ii))/sum(S1(:,idx2(ii)));
%     end
%end
