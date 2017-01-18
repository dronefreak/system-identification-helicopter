for i=2:20
    for j=1:13
        if rand>0.5
    population(i,j)=(population(i,j)+randsample(50,1))*rand;
        else
    population(i,j)=(population(i,j)-randsample(50,1))*rand;
        end
    end
end