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
% ALO, IALO,TALO, GA, PSO and FA based QAP works, 23.05.2017
% by U.Yuzgec

clc;
clear;
close all;
refresh = 50;

%% Problem Definition

model=CreateModel();

CostFunction=@(s) MyCost(s, model);        % Cost Function

nVar=model.m;       % Number of Decision Variables

VarSize=[1 nVar];   % Size of Decision Variables Matrix

VarMin=0;         % Lower Bound of Variables
VarMax=1;         % Upper Bound of Variables

MaxIt=1000;         % Maximum Number of Iterations
nPop=20;            % Population Size

empty_general_pop.Position=[];
general_pop=repmat(empty_general_pop,nPop,1);
for i=1:nPop
   general_pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
end

%% Firefly Algorithm Parameters

gamma=1;                               % Light Absorption Coefficient
beta0=2;                                  % Attraction Coefficient Base Value
alpha=0.2;                                % Mutation Coefficient
alpha_damp=0.98;                    % Mutation Coefficient Damping Ratio
delta=0.05*(VarMax-VarMin);     % Uniform Mutation Range
m=2;                                        % using for updating beta

if isscalar(VarMin) && isscalar(VarMax)
    dmax = (VarMax-VarMin)*sqrt(nVar);
else
    dmax = norm(VarMax-VarMin);
end

nMutation = 2;      % Number of Additional Mutation Operations

%% GA Parameters
pc=0.4;                         % Crossover Percentage
nc=2*round(pc*nPop/2);  % Number of Offsprings (Parents)

pm=0.8;                         % Mutation Percentage
nm=round(pm*nPop);      % Number of Mutants

beta=5;                          % Selection Pressure

%% PSO Parameters

w=1;                  % Inertia Weight
wdamp=0.99;     % Inertia Weight Damping Ratio
c1=1.5;              % Personal Learning Coefficient
c2=2.0;              % Global Learning Coefficient

% If you would like to use Constriction Coefficients for PSO,
% uncomment the following block and comment the above set of parameters.

% % Constriction Coefficients
% phi1=2.05;
% phi2=2.05;
% phi=phi1+phi2;
% chi=2/(phi-2+sqrt(phi^2-4*phi));
% w=chi;          % Inertia Weight
% wdamp=1;        % Inertia Weight Damping Ratio
% c1=chi*phi1;    % Personal Learning Coefficient
% c2=chi*phi2;    % Global Learning Coefficient

% Velocity Limits
VelMax=0.1*(VarMax-VarMin);
VelMin=-VelMax;

nParticleMutation = 1;      % Number of Mutations Performed on Each Particle
nGlobalBestMutation = 3;    % Number of Mutations Performed on Global Best

Options = {'Genetic Algorithm', 'Particle Swarm Optimization', 'Firefly Optimization',...
                'Antlion Optimization', 'Improved Antlion Optimization', 'Tournament Antlion Optimization','All Heuristic Algorithms'};

[Selection, Ok] = listdlg('PromptString', 'Select training method for ANFIS:', ...
                          'SelectionMode', 'single', ...
                          'ListString', Options);

pause(0.01);
          
if Ok==0
    return;
end

switch Selection
    case 1,
        %% GA Results
        fprintf('\nGenetic Algorithm Results\n');
        ga
        figure(11);
        plot(BestCost,'LineWidth',2, 'Color', 'm');
        xlabel('Iteration');
        ylabel('Best Cost');
        save sonuc_ga.mat      
    case 2, 
        %% PSO Results
        fprintf('\nPSO Results\n');
        pso
        figure(11);
        plot(BestCost,'LineWidth',2, 'Color', 'r');
        xlabel('Iteration');
        ylabel('Best Cost');
        save sonuc_pso.mat    
    case 3,
        %% Firefly Results
        fprintf('\nFirefly Results\n');
        fa
        figure(11);
        plot(BestCost,'LineWidth',2, 'Color', 'b');
        xlabel('Iteration');
        ylabel('Best Cost');
        save sonuc_fa.mat  
case 4,
        %% ALO Results
        fprintf('\nAntlion Results\n');
        alo
        figure(11);
        plot(BestCost,'LineWidth',2, 'Color', 'g');
        xlabel('Iteration');
        ylabel('Best Cost');
        save sonuc_alo.mat
case 5,
        %% IALO Results
        fprintf('\nImproved Antlion Results\n');
        ialo
        figure(11);
        plot(BestCost,'LineWidth',2, 'Color', 'k');
        xlabel('Iteration');
        ylabel('Best Cost');
        save sonuc_ialo.mat
case 6,
        %% TALO Results
        fprintf('\nTournament based Antlion Results\n');
        talo
        figure(11);
        plot(BestCost,'LineWidth',2, 'Color', 'y');
        xlabel('Iteration');
        ylabel('Best Cost');
        save sonuc_talo.mat
    case 7, %%% ALL RESULTS...
        %% Firefly Results
        fprintf('\nFirefly Results\n');
        fa
        figure(11);
        plot(BestCost,'LineWidth',2, 'Color', 'b');
        save sonuc_fa.mat

        %% GA Results
        fprintf('\nGenetic Algorithm Results\n');
        ga
        figure(11);
        hold on
        plot(BestCost,'LineWidth',2, 'Color', 'm');
        save sonuc_ga.mat

        %% ALO Results
        fprintf('\nAntlion Results\n');
        alo
        figure(11);
        hold on
        plot(BestCost,'LineWidth',2, 'Color', 'g');
        save sonuc_alo.mat

        %% IALO Results
        fprintf('\nImproved Antlion Results\n');
        ialo
        figure(11);
        hold on
        plot(BestCost,'LineWidth',2, 'Color', 'k');
        save sonuc_ialo.mat

        %% TALO Results
        fprintf('\nTournament based Antlion Results\n');
        talo
        figure(11);
        hold on
        plot(BestCost,'LineWidth',2, 'Color', 'y');
        save sonuc_talo.mat

        %% PSO Results
        fprintf('\nPSO Results\n');
        pso
        figure(11);
        hold on
        plot(BestCost,'LineWidth',2, 'Color', 'r');
        xlabel('Iteration');
        ylabel('Best Cost');
        grid on;
        legend('FA','GA','ALO','IALO','TALO','PSO');
        title('Quadratic Assignment Problem with Heuristic Algorithms');
        save sonuc_pso.mat

end

