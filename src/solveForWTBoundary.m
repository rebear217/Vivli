function y = solveForWTBoundary(P,I,J,guess)
    %find where the posterior of cluster I
    % is bigger than for cluster J
    foptions = optimoptions('fsolve','Display','none');
    y=fsolve(@(x)(makeFunction(P,I,x) - makeFunction(P,J,x)),guess,foptions);
end

function y=makeFunction(P,j,x)
    y = P(x(:));
    y = y(j);
end