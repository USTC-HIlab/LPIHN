function p=randomwalk(p0,backprobability,M)
p1=p0;
p=(1-backprobability)*M'*p1+backprobability*p0;
while max(p-p1)>10^(-10)
    p1=p;
    p=(1-backprobability)*M'*p1+backprobability*p0;
end