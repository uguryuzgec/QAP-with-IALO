clc
figure
[ugur, index_ugur]=min(GlobalBest_GA(end,:));
x=model.x;
y=model.y;
        plot(x((1:end)),y((1:end)),'ro');     % Not-Assigned Locations
        hold on
        plot(x(GlobalSolution_GA(index_ugur,:)),y(GlobalSolution_GA(index_ugur,:)),'ks',...
        'MarkerSize',10,...
        'MarkerFaceColor','c'); 
        for i=1:model.n
        text(x(GlobalSolution_GA(index_ugur,i))+1,y(GlobalSolution_GA(index_ugur,i))-2,num2str(i),'FontSize',9);
        end
    xlabel('x');
    ylabel('y');
    hold off;
    grid on;
    axis equal;
