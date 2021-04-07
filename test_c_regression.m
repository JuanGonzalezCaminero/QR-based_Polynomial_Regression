
%Calls the C regression function with different degree polynomials 
%for the same set of points
function []=test_c_regression()
    axis([0 1 0 1]);
    [x,y]=ginput();
    
    for degree=1:7
        
        sol=qr_polynomial_regression_c(x,y,degree);
        
        %Plot the points and the polynomial
        xp=0:0.01:1;
        yp=polyval(sol, xp);
        
        figure();
        plot(xp, yp);
        axis([0 1 0 1]);
        hold on;
        plot(x,y,'*');
        title("Polynomial of degree " + degree);
    end