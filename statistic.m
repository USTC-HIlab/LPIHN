function [tp,tn,fp,fn]=statistic(label,pre_label)
tp=0;
tn=0;
fp=0;
fn=0;
for i=1:length(label)
    if((label(i)==1)&&(pre_label(i)==1))
        tp=tp+1;
    elseif((label(i)==1)&&(pre_label(i)==0))
        fn=fn+1;
    elseif((label(i)==0)&&(pre_label(i)==1))
        fp=fp+1;
    elseif((label(i)==0)&&(pre_label(i)==0))
        tn=tn+1;
    else
        disp('error');
    end
end