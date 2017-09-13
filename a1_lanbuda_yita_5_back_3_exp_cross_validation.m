clc
clear

%��һ������֤������Ҫ����������������
load('result/correct_exp_random_walk_corr_score_5_5_3.mat'); % ����a1_lanbuda_yita_5_back_3_exp_random_walk.m����ó�������
load('data/pro_delet_lnc_matrix.mat'); %lncRNA��protein�໥���þ����д���protein���д���lncRNA��������Ϊɾ��lncRNA�����ڵ��ľ���
load('data/lnc_delet_corr_coeff.mat'); %��PCC�������lncRNA����������Ծ���
load('data/PPI_delet_original_matrix.mat'); %PPI����

AP=normalization(PPI_delet_original_matrix);%AP:��һ����ĵ������໥���þ���
AL=lnc_delet_corr_coeff;  %AL:lncRNA����������Ծ���
coeff=correct_exp_random_walk_corr_score_5_5_3;
B=pro_delet_lnc_matrix;
[n,m]=size(B);
gold_data=reshape(B,n*m,1);
lanbuda=0.5;
yita=0.5; %lanbuda��yita�����������������ʼ���ʣ�һ��������ѡȡ0.5����
backprobability=0.3; %���ظ���

%��һ������֤����
for h=1:m
    h
    coeff_final=coeff(((h-1)*n+1):h*n,1);
    index=find(1==B(:,h));
    if ~isempty(index)
        for z=1:length(index)
            B(index(z),h)=0;
			
            MP2=(1-lanbuda)*bsxfun(@rdivide,AP,sum(AP,2));
            index11=find(0==sum(B,2));
            MP2(index11,:)=bsxfun(@rdivide,AP(index11,:),sum(AP(index11,:),2));
            MP2(find(0==sum(AP,2)),:)=0;

            ML2=(1-lanbuda)*bsxfun(@rdivide,AL,sum(AL,2));
            index22=find(0==sum(B));
            ML2(index22,:)=bsxfun(@rdivide,AL(index22,:),sum(AL(index22,:),2));

            MPL=lanbuda*bsxfun(@rdivide,B,sum(B,2));
            MPL(find(isnan(MPL)==1))=0;%protein����ת�Ƶ�lncRNA�ģ�������ĸ���ΪNAN����Ϊsum(B,2)=0��,�ⲽ��NAN��0�滻���ò������NAN������û���⣬MLPͬ��
            MLP=lanbuda*bsxfun(@rdivide,B,sum(B));
            MLP(find(isnan(MLP)==1))=0;
            MLP=MLP';
            M1=[MP2,MPL];
            M2=[MLP,ML2];
            M=[M1;M2];%ת�ƾ���
			
            u0=B(:,h);
            v0=zeros(m,1);
            v0(h,1)=1;
            if sum(u0)~=0
               u0=u0/sum(u0);
               p0=[(1-yita).*u0;yita.*v0];
            else
                p0=[u0;v0];  %��ʼ���ʾ���
            end 
			
             P=randomwalk(p0,backprobability,M);  %�˴��õ���randomwalk�����ڸ��ļ�����
             coeff_final(index(z),1)=P(index(z),1); 
             B(index(z),h)=1;
        end
    end
    coeff_all((n*(h-1)+1):n*h,1)=coeff_final;
end

save result/coeff_all coeff_all
auc_2=roc(coeff_all,gold_data,'r');