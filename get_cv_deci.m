function [deci,label_y] = get_cv_deci(prob_y,prob_x,param,nr_fold)
l=length(prob_y);
deci = ones(l,1);
label_y = ones(l,1);
% rand_ind = randperm(l);
K_fold_indx = CV_split_data(prob_y,nr_fold);
for i=1:nr_fold % Cross training : folding
    train_ind=K_fold_indx{i,2};
    test_ind=K_fold_indx{i,1};
%     train_data=prob_x(train_index,:);
%     test_data=prob_x(test_index,:);
%     train_label=label(train_index);
    
%     test_ind=rand_ind((floor((i-1)*l/nr_fold)+1:floor(i*l/nr_fold))');
%     train_ind =(1:l)';
%     train_ind(test_ind) = [];
    model = svmtrain(prob_y(train_ind),prob_x(train_ind,:),param);
    [predict_label,mse,subdeci] = svmpredict(prob_y(test_ind),prob_x(test_ind,:),model);
    if model.Label(1)==1
        mod=1;
    else 
        mod=-1;
    end
    deci(test_ind) = subdeci*mod;
    label_y(test_ind) = prob_y(test_ind);
end
end