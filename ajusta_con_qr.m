
function [sol]=ajusta_con_qr()
    axis([0 1 0 1]);
    [x,y]=ginput();
    %plot(x,y,'*')
    
    for grado=1:7
        M=[];
        for i=0:grado
            M(:,grado+1-i)=x.^i;
        end
        
        [Q,R]=qr(M);
    
        sol=R\(Q'*y);
        
        disp("Coeficientes del polinomio de grado " + grado + " generado:")
        sol
        
        %Para plottear la curva generada y nuestros puntos
        xp=0:0.01:1;
        yp=polyval(sol, xp);
        
        figure();
        plot(xp, yp);
        axis([0 1 0 1]);
        hold on;
        plot(x,y,'*');
    end
end














