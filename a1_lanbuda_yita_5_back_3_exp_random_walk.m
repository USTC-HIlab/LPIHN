tic
clc
clear

%��Ҫ����������������
load('data/pro_delet_lnc_matrix.mat');%lncRNA-protein�໥���þ��󣬸�����Ϊɾ��lncRNA�����ڵ��ľ���
B=pro_delet_lnc_matrix;
load('data/lnc_delet_corr_coeff.mat');%����PCC�������lncRNA����������Ծ���
load('data/PPI_delet_original_matrix.mat');%protein�໥���ã�PPI������

AP=normalization(PPI_delet_original_matrix);%PPI�Ĺ�һ����AP��PPI����Ȩ�ؾ���
AL=lnc_delet_corr_coeff;% AL:lncRNA����Ȩ�أ������ԣ�����
[n,m]=size(B);
lanbuda=0.5;
yita=0.5;
backprobability=0.3;%lanbuda,yita,backprobability��Ϊ���������Լ����ã�0.1-0.9֮�䡣

%���¿�ʼ���ܵ�ת�ƾ���M(��MP2,MPL,MLP,ML2���).
MP2=(1-lanbuda)*bsxfun(@rdivide,AP,sum(AP,2));
index11=find(0==sum(B,2));
MP2(index11,:)=bsxfun(@rdivide,AP(index11,:),sum(AP(index11,:),2));

ML2=(1-lanbuda)*bsxfun(@rdivide,AL,sum(AL,2));
index22=find(0==sum(B));
ML2(index22,:)=bsxfun(@rdivide,AL(index22,:),sum(AL(index22,:),2));

MPL=lanbuda*bsxfun(@rdivide,B,sum(B,2));
MPL(find(isnan(MPL)==1))=0;%p����ת�Ƶ�l�ģ�������ĸ���ΪNAN����Ϊsum(B,2)=0��,�ⲽ��NAN��0�滻,�ò������NAN������û���⣬MLPͬ��
MLP=lanbuda*bsxfun(@rdivide,B,sum(B));
MLP(find(isnan(MLP)==1))=0;
MLP=MLP';
M1=[MP2,MPL];
M2=[MLP,ML2];
M=[M1;M2];%�ܵ�ת�ƾ���M

%�������ܵĳ�ʼ����P0
u0=bsxfun(@rdivide,B,sum(B));
v0=eye(m,m);
p0=[(1-yita)*u0;yita*v0];%�˴��õ���p0Ϊ�ܵĳ�ʼ���ʾ���ÿһ�д���һ��lncRNA��Ӧ�ĳ�ʼ����

%���¿�ʼ�������
for k=1:m
    k
p1=p0(:,k);
p=(1-backprobability)*M'*p1+backprobability*p0(:,k);
while max(p-p1)>10^(-10)
    p1=p;
    p=(1-backprobability)*M'*p1+backprobability*p0(:,k);
end

    original_random_walk_score((n*(k-1)+1:n*k),1)=p(1:n,1);
end
correct_exp_random_walk_corr_score_5_5_3=original_random_walk_score;
save result/correct_exp_random_walk_corr_score_5_5_3 correct_exp_random_walk_corr_score_5_5_3
toc


