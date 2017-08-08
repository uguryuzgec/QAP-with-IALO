%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPAP104
% Project Title: Solving QAP using PSO and Firefly Algorithm (FA)
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

function PlotSolution(s, model)

    [~, p] = sort(s);

    n=model.n;
    x=model.x;
    y=model.y;
    
    plot(x(p(1:n)),y(p(1:n)),'kd',...
        'MarkerSize',10,...
        'MarkerFaceColor','m');                 % Assigned Locations
    
    hold on;

    plot(x(p(n+1:end)),y(p(n+1:end)),'bo');     % Not-Assigned Locations
    
    for i=1:n
        text(x(p(i))+1,y(p(i))-2,num2str(i),'FontSize',9);
    end
    xlabel('x');
    ylabel('y');
    hold off;
    grid on;
    axis equal;
    
end