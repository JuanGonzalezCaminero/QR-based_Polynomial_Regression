
function [sol]=qr_polynomial_regression_matlab(x, y, n)
    M=[];
    for i=0:n
        M(:,n+1-i)=x.^i;
    end

    [Q,R]=qr(M);

    sol=R\(Q'*y);
    
    disp("R =")
    disp(R)
    disp("Q =")
    disp(Q)
    disp("Rx = Q'b")
    disp(sol)
end














