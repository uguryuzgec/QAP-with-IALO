function [GBestmem,GBestval,nfeval1,Gmin,Giter,GSure]=ALO_Updated_Tournament(fname,VTR,D,XVmin,XVmax,NP,Max_Iter,refresh,tekrar,solution)
% AOL Updated by Haydar Kilic
%E_AL_pos = Elite Antlion Position
%E_AL_fit = Elite Antlion Fitness
%S_AL_pos = Sorted Antlion Position
%S_AL_fit = Sorted Antlion Fitness
%AL_pos = Antlion Position
%A_pos = Ant Position
%AL_fit = Antlion Fitness
%A_fit = Ant Fitness
%NP = Population Size (Search Agent Number)
%XVmin = Lower Boundary of Population Size
%XVmax = Upper Boundary of Population Size
%D = Dension of Variables
%fname = Objective Function
%RA = Random Walk Around the Selected Antlion by RWh
%RE = Random Walk Around the Elite Antlion
%RWh_AL_Index = R. Wheel Selected Antlion Index
%All_Pop_pos = All Population Position so far
%All_Pop_fit = All Population Fitness so far
%S_All_Pop_pos = Sorted All Population Position so far
%S_All_Pop_fit = Sorted All Population Fitness so far
%function [E_AL_fit,E_AL_pos,Conv_Curve]=ALO_Updated_Rwh(NP, Max_Iter, XVmin, XVmax, D, fname)
%-----Check input variables---------------------------------------------
err=[];
if nargin<1, error('ALO 1st argument must be function name'); else 
  if exist(fname)<1; err(1,length(err)+1)=1; end; end;
if nargin<2, VTR = 1.e-6; else 
  if length(VTR)~=1; err(1,length(err)+1)=2; end; end;
if nargin<3, D = 2; else
  if length(D)~=1; err(1,length(err)+1)=3; end; end; 
if nargin<4, XVmin = [-2 -2];else
  if length(XVmin)~=D; err(1,length(err)+1)=4; end; end; 
if nargin<5, XVmax = [2 2]; else
  if length(XVmax)~=D; err(1,length(err)+1)=5; end; end; 
if nargin<6, NP = 10*D; else
  if length(NP)~=1; err(1,length(err)+1)=6; end; end; 
if nargin<7, Max_Iter = 200; else
  if length(Max_Iter)~=1; err(1,length(err)+1)=7; end; end; 
if nargin<8, refresh = 10; else
  if length(refresh)~=1; err(1,length(err)+1)=8; end; end; 
if nargin<9, tekrar = 10; else
  if length(tekrar)~=1; err(1,length(err)+1)=9; end; end;
if ~isempty(err)% if length(err)>0
  fprintf(stdout,'error in parameter %d\n', err);
  usage('alo (string,scalar,scalar,vector,vector,any,integer,integer,scalar,scalar,integer,integer)');    	
end

if (NP < 5)
   NP=5;
   fprintf(1,' NP increased to minimal value 5\n');
end

if (Max_Iter <= 0)
   Max_Iter = 200;
   fprintf(1,'itermax should be > 0; set to default value 200\n');
end
refresh = floor(refresh);
tournSize=2; %%%%%%%%%%%%%%%%%%%%%%%%%%% Turnuva boyutu
%Initialization Positions of Ants and AntLions
%Initialize the Variables
S_AL_pos=zeros(NP, D);
E_AL_pos=zeros(1,D);
E_AL_fit=inf;
% % % Conv_Curve=zeros(1,Max_Iter);
A_fit=zeros(1,NP);
AL_fit=zeros(1,NP);
GlobalMins=zeros(tekrar,Max_Iter-1); % herbir tekrar icin bulunan tum en iyi minimumlar
Globalbest=zeros(tekrar,D); % herbir tekrar icin bulunan en iyi cozumler (x,y)
GBestval=zeros(tekrar,1); % herbir tekrar icin bulunan en iyi minimum deger
GIter=zeros(tekrar,1); % herbir tekrar icin bulunan iterasyon sayisi
GSure=zeros(tekrar,1);
nfeval = 0;

for r=1:tekrar
fprintf(1,'\n%d. tekrar\n',r);

AL_pos=initialization(NP, D, XVmax, XVmin);
A_pos=initialization(NP, D, XVmax, XVmin);
S_AL_pos=zeros(NP, D);
E_AL_pos=zeros(1,D);
E_AL_fit=inf;
A_fit=zeros(1,NP);
AL_fit=zeros(1,NP);
%Calculation of Initial AntLions' Fitness and Sort Them
for i=1:NP
    AL_fit(1,i)=feval(fname,AL_pos(i,:));
    nfeval=nfeval+1;
end

[S_AL_fit,S_Index]=sort(AL_fit);
S_AL_pos=AL_pos(S_Index,:); % You can use this for (1)
E_AL_pos=S_AL_pos(1,:);
E_AL_fit=S_AL_fit(1);

% Iterations or pseudo time marching
%-------------------------------------------------------------------------------------------------------
% % %   x1 = linspace(XVmin(1), XVmax(1), 101);
% % %   x2 = linspace(XVmin(2), XVmax(2), 101);
% % %   x3 = zeros(length(x1), length(x2));
% % %         
% % %         % simply loop through the function (most functions expect 
% % %         % [N x 2] vectors as input, so meshgrid does not work)
% % %         for i = 1:length(x1)
% % %             for j = 1:length(x2)
% % %                 x3(i, j) = feval(fname,[x1(i), x2(j)]);
% % %             end
% % %         end
% % % 		figure;view(-40, 30);contour(x1', x2', x3); hold on; 
% % %         plot(S_AL_pos(:,1),S_AL_pos(:,2),'b.','MarkerSize',20);
% % %         plot(solution(:,1), solution(:,2), 'r.', 'MarkerSize', 20);
% % %         drawnow; pause(.1);
% % %         hold off
        
%-------------------------------------------------------------------------------------------------------

% Now Main Loop Start from Second Iteration Since First Iteration was Dedicated to Calculating the Fitness of Antlions
delta_val=S_AL_fit;
durdurma_katsayisi=delta_val(NP)-delta_val(1);
Current_Iter = 1;
NNN=Max_Iter*0.2;
saat=tic; %zamanlayici baslat....
while ((Current_Iter < Max_Iter) && (abs(durdurma_katsayisi) > VTR))
    
   % This Loop Simulates Random Walks
   % X ve Y ANTLION RANDWALK BÝR KERE DONGU DISINDA HESAPLANIYOR...
	X = [0 cumsum(2*(rand(NNN,1)>0.5)-1)'];
    Y = [0 cumsum(2*(rand(NNN,1)>0.5)-1)'];
     ax=min(X);  bx=max(X);
     ay=min(Y);  by=max(Y);
   %Tournament Selection
    
    % Repeat the process "tournSize" times
    winners = [];
    n=length(AL_fit);
for i = 1:tournSize
    %Create a random set of competitors
    shuffleOrder = randperm(n);
    competitors = reshape(shuffleOrder, tournSize, n/tournSize)';
    %The winner is the competitor with best fitness
    [winFit, winID] = min(AL_fit(competitors),[],2);
    idMap = (0:tournSize-1)*n/tournSize;
    idMap1 = idMap(winID) + (1:size(competitors,1));
    winners = [winners; competitors(idMap1)];
end
    Select_Index=winners;
       
    for k=1:NP

	 antlion = AL_pos(Select_Index(k),:);
		
		% I is the ratio in Equations (2.10) and (2.11)
		I=1; 
		% Her seferinde tüm if koþullarýný sorgulamasýna gerek yok, bir tek if koþulu gerçekleþir. 
		% Dolayýsýyla I güncellemesi elseif merdiveni haline getirildi.
		if Current_Iter>Max_Iter*(0.95)
			I=1+1000000*(Current_Iter/Max_Iter);
            elseif Current_Iter>Max_Iter*(0.9)
				I=1+100000*(Current_Iter/Max_Iter);
				elseif Current_Iter>Max_Iter*(0.75)
					I=1+10000*(Current_Iter/Max_Iter);
					elseif Current_Iter>Max_Iter*(0.5) 
						I=1+1000*(Current_Iter/Max_Iter); 
						elseif Current_Iter>Max_Iter*(0.1)
							I=1+100*(Current_Iter/Max_Iter);
		end
		% Dicrease boundaries to converge towards antlion
		
        lb_gecici=XVmin/(I); % Equation (2.10) in the paper
		ub_gecici=XVmax/(I); % Equation (2.11) in the paper

		% Move the interval of [lb ub] around the antlion [lb+anlion ub+antlion]
        option = rand;
		if option>0.75
			lb_gecici = antlion + lb_gecici; % Equation (2.8) in the paper
            ub_gecici = antlion + ub_gecici;% Equation (2.9) in the paper
        elseif option>0.5
			ub_gecici = antlion - ub_gecici; % Equation (2.9) in the paper
			lb_gecici = antlion - lb_gecici; % Equation (2.8) in the paper
        elseif option>0.25
            ub_gecici = -antlion + ub_gecici; % Equation (2.9) in the paper
			lb_gecici = -antlion + lb_gecici; % Equation (2.8) in the paper
        else
            ub_gecici = -antlion - ub_gecici; % Equation (2.9) in the paper
			lb_gecici = -antlion - lb_gecici; % Equation (2.8) in the paper            
		end

		% Creates n random walks and normalize accroding to lb and ub vectors 
for i=1:D
% % %             X=[0 cumsum(2*(rand(NNN,1)>0.5)-1)'];
		c=lb_gecici(i); 
        d=ub_gecici(i);
        RA(:,i) = ((X-ax).*(d-c))./(bx-ax)+c;  % Equation (2.7) in the paper
end
		
		% Random Walk around the Pseudo Elite Antlion
		antlion = E_AL_pos; % Elit Antlion icin randomwalk ile RE bulunmasi
		
        % Dicrease boundaries to converge towards antlion
		lb_gecici=XVmin/(I); % Equation (2.10) in the paper
		ub_gecici=XVmax/(I); % Equation (2.11) in the paper
		
        % Move the interval of [lb ub] around the antlion [lb+anlion ub+antlion]
        option = rand;
		if option>0.75
			lb_gecici = antlion + lb_gecici; % Equation (2.8) in the paper
            ub_gecici = antlion + ub_gecici;% Equation (2.9) in the paper
        elseif option>0.5
			ub_gecici = antlion - ub_gecici; % Equation (2.9) in the paper
			lb_gecici = antlion - lb_gecici; % Equation (2.8) in the paper
        elseif option>0.25
            ub_gecici = -antlion + ub_gecici; % Equation (2.9) in the paper
			lb_gecici = -antlion + lb_gecici; % Equation (2.8) in the paper
        else
            ub_gecici = -antlion - ub_gecici; % Equation (2.9) in the paper
			lb_gecici = -antlion - lb_gecici; % Equation (2.8) in the paper            
        end
        
    % Creates n random walks and normalize accroding to lb and ub vectors
for i=1:D
% % %     Y=[0 cumsum(2*(rand(NNN,1)>0.5)-1)'];
		c=lb_gecici(i);
		d=ub_gecici(i);
        RE(:,i) = ((Y-ay).*(d-c))./(by-ay)+c;  
end	
			   
		% Ant Positions around the RWh Selected and Pseudo Elite Antlions
            A_pos(k,:) = (RA(randi(NNN,1,1),:)+RE(randi(NNN,1,1),:))/2;
    end % for k=1:NP
	
	%If the ants with respect to antlions(elite and roulette wheel selection) go beyond the boundries,
    %the antlions' positions and fitnesses has to be updated. Because, if an ant becomes fitter than an antlion,
    % we assume it was caught by the antlion and the antlion updates its position to build the new trap. The following code shows how to do that.
    %We have to check boundries first, if the ant position greater than the upper boundary or smaller than lower boundary, its position is zero :
    for i=1:NP,
         for j=1:D,
            if (A_pos(i,j)>XVmax(j)) || (A_pos(i,j)<XVmin(j))
                A_pos(i,j)=XVmin(j) + rand()*(XVmax(j) - XVmin(j));
            end
        end
	end
    %we have pseudo elite antlion, roulette wheel selected antlion, and ants that will hunt for them.
    %In this condition, one of antlions has to be elite, we use following code to find it.

for i=1:NP
       A_fit(1,i)=feval(fname, A_pos(i,:));
       nfeval=nfeval+1;
       if (A_fit(i) < AL_fit(i))  % if competitor is better than value in "cost array"
         AL_pos(i,:) =A_pos(i,:);  % replace old vector with new one (for new iteration)
         AL_fit(i)   = A_fit(i);  % save value in "cost array"
            if A_fit(i)<E_AL_fit
                E_AL_pos = A_pos(i,:);
                E_AL_fit = A_fit(i);
            end
       end
end      

% % %    if (refresh > 0)
% % %     if (rem(Current_Iter,refresh) == 0)
% % %         figure(2);contour(x1', x2', x3); 
% % %         hold on; 
% % %         plot(solution(:,1), solution(:,2), 'r.', 'MarkerSize', 20);
% % %         plot(S_AL_pos(:,1),S_AL_pos(:,2),'b.','MarkerSize',20);
% % %         plot(A_pos(:,1),A_pos(:,2),'r*','MarkerSize',8);
% % %         drawnow; pause(.1);
% % %         hold off
% % %     end
% % %    end
    % Updating Convergence Curve
% % %     Conv_Curve(Current_Iter)=E_AL_fit;
    %----Output section----------------------------------------------------------
  if (refresh > 0)
    if (rem(Current_Iter,refresh) == 0)
       fprintf(1,'\nIteration: %d,  Best: %f,  NP: %d\n',Current_Iter,E_AL_fit,NP);
       for n=1:D
         fprintf(1,'best(%d) = %f, ',n,E_AL_pos(n));
       end
    end
  end
  
    GlobalMins(r,Current_Iter)=E_AL_fit;
    [delta_val]=sort(AL_fit);
    durdurma_katsayisi=delta_val(NP)-delta_val(1);
    Current_Iter=Current_Iter+1;
% % %     durdurma_katsayisi=1;
end % while ...

sonsaat=toc(saat);
Globalbest(r,:)=E_AL_pos;
GlobalBestval(r,1)=E_AL_fit;
GlobalIter(r,1)=Current_Iter-1;
GlobalSure(r,1)=sonsaat;
end %end for (r=1:tekrar) ...
fprintf(1,'\nOrtalama Sonuc= %e\n',mean(GlobalBestval));
fprintf(1,'\nStandart Sapma= %e\n',std(GlobalBestval));
fprintf(1,'\nOrtalama Cozum= %f,%f\n',sum(Globalbest(:,1))/tekrar,sum(Globalbest(:,2))/tekrar);
fprintf(1,'\nOrtalama iterasyon= %f\n',mean(GlobalIter));
fprintf(1,'\nOrtalama Sure= %fsn\n',mean(GlobalSure));
nfeval1=nfeval/tekrar;
GBestmem = Globalbest;
GBestval = GlobalBestval;
Gmin=GlobalMins;
Giter=GlobalIter;
GSure=GlobalSure;

    