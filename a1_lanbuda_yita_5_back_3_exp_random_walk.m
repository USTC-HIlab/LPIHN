tic
clc
clear

%需要输入以下三种数据
load('data/pro_delet_lnc_matrix.mat');%lncRNA-protein相互作用矩阵，该数据为删除lncRNA孤立节点后的矩阵
B=pro_delet_lnc_matrix;
load('data/lnc_delet_corr_coeff.mat');%利用PCC计算出的lncRNA表达谱相似性矩阵
load('data/PPI_delet_original_matrix.mat');%protein相互作用（PPI）矩阵。

AP=normalization(PPI_delet_original_matrix);%PPI的归一化。AP：PPI网络权重矩阵。
AL=lnc_delet_corr_coeff;% AL:lncRNA网络权重（相似性）矩阵
[n,m]=size(B);
lanbuda=0.5;
yita=0.5;
backprobability=0.3;%lanbuda,yita,backprobability均为参数，可自己设置，0.1-0.9之间。

%以下开始求总的转移矩阵M(由MP2,MPL,MLP,ML2组成).
MP2=(1-lanbuda)*bsxfun(@rdivide,AP,sum(AP,2));
index11=find(0==sum(B,2));
MP2(index11,:)=bsxfun(@rdivide,AP(index11,:),sum(AP(index11,:),2));

ML2=(1-lanbuda)*bsxfun(@rdivide,AL,sum(AL,2));
index22=find(0==sum(B));
ML2(index22,:)=bsxfun(@rdivide,AL(index22,:),sum(AL(index22,:),2));

MPL=lanbuda*bsxfun(@rdivide,B,sum(B,2));
MPL(find(isnan(MPL)==1))=0;%p不能转移到l的，计算出的概率为NAN（因为sum(B,2)=0）,这步将NAN用0替换,该步骤出现NAN正常，没问题，MLP同理
MLP=lanbuda*bsxfun(@rdivide,B,sum(B));
MLP(find(isnan(MLP)==1))=0;
MLP=MLP';
M1=[MP2,MPL];
M2=[MLP,ML2];
M=[M1;M2];%总的转移矩阵M

%以下求总的初始概率P0
u0=bsxfun(@rdivide,B,sum(B));
v0=eye(m,m);
p0=[(1-yita)*u0;yita*v0];%此处得到的p0为总的初始概率矩阵，每一列代表一个lncRNA对应的初始概率

%以下开始随机游走
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


