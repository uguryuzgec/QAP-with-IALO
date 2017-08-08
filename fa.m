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

%% Initialization

% Empty Firefly Structure
firefly.Position=[];
firefly.Cost=[];
firefly.Sol=[];

% Initialize Population Array
pop=repmat(firefly,nPop,1);

% Initialize Best Solution Ever Found
BestSol.Cost=inf;

% Create Initial Fireflies
for i=1:nPop
%    pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
    pop(i).Position=general_pop(i).Position;
   [pop(i).Cost, pop(i).Sol]=CostFunction(pop(i).Position);
   
   if pop(i).Cost<=BestSol.Cost
       BestSol=pop(i);
   end
end

% Array to Hold Best Cost Values
BestCost=zeros(MaxIt,1);

%% Firefly Algorithm Main Loop

for it=1:MaxIt
    
    newpop=repmat(firefly,nPop,1);
    for i=1:nPop
        newpop(i).Cost = inf;
        for j=1:nPop
            if pop(j).Cost < pop(i).Cost
                rij=norm(pop(i).Position-pop(j).Position)/dmax;
                beta=beta0*exp(-gamma*rij^m);
                e=delta*unifrnd(-1,+1,VarSize);
                %e=delta*randn(VarSize);
                
                newsol.Position = pop(i).Position ...
                                + beta*rand(VarSize).*(pop(j).Position-pop(i).Position) ...
                                + alpha*e;
                
                newsol.Position=max(newsol.Position,VarMin);
                newsol.Position=min(newsol.Position,VarMax);
                
                [newsol.Cost, newsol.Sol]=CostFunction(newsol.Position);
                
                if newsol.Cost <= newpop(i).Cost
                    newpop(i) = newsol;
                    if newpop(i).Cost<=BestSol.Cost
                        BestSol=newpop(i);
                    end
                end
                
            end
        end
        
        % Perform Mutation
        for k=1:nMutation
            newsol.Position = Mutate(pop(i).Position);
            [newsol.Cost, newsol.Sol]=CostFunction(newsol.Position);
            if newsol.Cost <= newpop(i).Cost
                newpop(i) = newsol;
                if newpop(i).Cost<=BestSol.Cost
                    BestSol=newpop(i);
                end
            end
        end
                
    end %  for i=1:nPop...
    
    % Merge
    pop=[pop
         newpop];  %#ok
    
    % Sort
    [~, SortOrder]=sort([pop.Cost]);
    pop=pop(SortOrder);
    
    % Truncate
    pop=pop(1:nPop);
    
    % Store Best Cost Ever Found
    BestCost(it)=BestSol.Cost;
    
    % Show Iteration Information
     if (rem(it,refresh) == 0)
        disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
        figure(1);
        PlotSolution(BestSol.Position, model);
        str = sprintf('QAP with FA : %.3f',BestSol.Cost');
%         title(str);
        text(-5,15,str);
        pause(0.01);
     end
    % Damp Mutation Coefficient
    alpha = alpha*alpha_damp;
end %for it=1:MaxIt...
