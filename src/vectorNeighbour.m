function NN = vectorNeighbour(v,times)
    v = v(:)';
    n = length(v);
    N = repmat(v,times,1);
    R = round(10*rand(times,n));
    N = N + R;
    NN = N(1:times,:);
    NN = unique(NN,'rows');
    NN(1,:) = v;
end