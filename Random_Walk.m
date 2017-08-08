%___________________________________________________________________%
%  Ant Lion Optimizer (ALO) source codes demo version 1.0           %
%                                                                   %
%  Developed in MATLAB R2011b(7.13)                                 %
%                                                                   %
%  Author and programmer: Seyedali Mirjalili                        %
%                                                                   %
%         e-Mail: ali.mirjalili@gmail.com                           %
%                 seyedali.mirjalili@griffithuni.edu.au             %
%                                                                   %
%       Homepage: http://www.alimirjalili.com                       %
%                                                                   %
%   Main paper:                                                     %
%                                                                   %
%   S. Mirjalili, The Ant Lion Optimizer                            %
%   Advances in Engineering Software , in press,2015                %
%   DOI: http://dx.doi.org/10.1016/j.advengsoft.2015.01.010         %
%                                                                   %
%___________________________________________________________________%

% This function creates random walks

function [Rws]=Random_Walk(D,Max_Iter,XVmin, XVmax,antlion,Current_Iter,NNN)
% % if size(lb,1) ==1 && size(lb,2)==1 %Check if the bounds are scalar
% %     lb=ones(1,Dim)*lb;
% %     ub=ones(1,Dim)*ub;
% % end
% % 
% % if size(lb,1) > size(lb,2) %Check if boundary vectors are horizontal or vertical
% %     lb=lb';
% %     ub=ub';
% % end
% % 
% % I=1; % I is the ratio in Equations (2.10) and (2.11)
% % 
% % if current_iter>max_iter/10
% %     I=1+100*(current_iter/max_iter);
% % end
% % 
% % if current_iter>max_iter/2
% %     I=1+1000*(current_iter/max_iter);
% % end
% % 
% % if current_iter>max_iter*(3/4)
% %     I=1+10000*(current_iter/max_iter);
% % end
% % 
% % if current_iter>max_iter*(0.9)
% %     I=1+100000*(current_iter/max_iter);
% % end
% % 
% % if current_iter>max_iter*(0.95)
% %     I=1+1000000*(current_iter/max_iter);
% % end
% % 
% % 
% % % Dicrease boundaries to converge towards antlion
% % lb=lb/(I); % Equation (2.10) in the paper
% % ub=ub/(I); % Equation (2.11) in the paper
% % 
% % % Move the interval of [lb ub] around the antlion [lb+anlion ub+antlion]
% % if rand<0.5
% %     lb=lb+antlion; % Equation (2.8) in the paper
% % else
% %     lb=-lb+antlion;
% % end
% % 
% % if rand>=0.5
% %     ub=ub+antlion; % Equation (2.9) in the paper
% % else
% %     ub=-ub+antlion;
% % end
% % 
% % % This function creates n random walks and normalize accroding to lb and ub
% % % vectors 
% % for i=1:Dim
% %     X = [0 cumsum(2*(rand(max_iter,1)>0.5)-1)']; % Equation (2.1) in the paper
% %     %[a b]--->[c d]
% %     a=min(X);
% %     b=max(X);
% %     c=lb(i);
% %     d=ub(i);      
% %     X_norm=((X-a).*(d-c))./(b-a)+c; % Equation (2.7) in the paper
% %     RWs(:,i)=X_norm;
% % end
I=1; 
		% Her seferinde tüm if koþullarýný sorgulamasýna gerek yok, bir tek if koþulu gerçekleþir. 
		% Dolayýsýyla I güncellemesi elseif merdiveni haline getirildi.
		if Current_Iter>Max_Iter*(0.95)
			I=1+1000000*(Current_Iter/Max_Iter);
            elseif Current_Iter>Max_Iter*(0.9)
				I=1+100000*(Current_Iter/Max_Iter);
				elseif Current_Iter>Max_Iter*(3/4)
					I=1+10000*(Current_Iter/Max_Iter);
					elseif Current_Iter>Max_Iter/2
						I=1+1000*(Current_Iter/Max_Iter); 
						elseif Current_Iter>Max_Iter/10;
							I=1+100*(Current_Iter/Max_Iter);
% % % 							else
% % % 								I=1+10*(Current_Iter/Max_Iter);
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
			ub_gecici = antlion - ub_gecici;
            
			lb_gecici = antlion - lb_gecici; % Equation (2.8) in the paper
        elseif option>0.25
            ub_gecici = -antlion + ub_gecici; % Equation (2.9) in the paper
			lb_gecici = -antlion +lb_gecici; % Equation (2.8) in the paper
        else
            ub_gecici = -antlion - ub_gecici; % Equation (2.9) in the paper
			lb_gecici = -antlion - lb_gecici; % Equation (2.8) in the paper
		end

		% Creates n random walks and normalize accroding to lb and ub vectors 
		for i=1:D
			X = [0 cumsum(2*(rand(NNN,1)>0.5)-1)']; % Equation (2.1) in the paper
			%[a b]--->[c d]
            %----
% 			X=[0 cumsum(2*(trnd(1,Max_Iter,1)>0.5)-1)'];
            a=min(X);
			b=max(X);
			c=lb_gecici(i);
			d=ub_gecici(i);
% % %             eps=rand;
% % %             Norm_Rwalk = (X-a)/(b-a);
% % %             th = 2*atan(((1-eps)/(1+eps))*tan(rand(length(Norm_Rwalk),1)-0.5*(ones(length(Norm_Rwalk),1)))*pi);            
% % %             %lambda = abs(min(lb_gecici)-max(ub_gecici))/2;
% % %             if Walk_Type==1
% % %                 Rws(:,i)=Norm_Rwalk'.*(cos(th)+mean(Norm_Rwalk).*rand(length(th),1)).*(d-c)+c.*ones(length(th),1); %9
% % %             else               
% % %                 Rws(:,i)=Norm_Rwalk'.*(sin(th)+mean(Norm_Rwalk).*rand(length(th),1)).*(d-c)+c.*ones(length(th),1); %9
% % %             end
            Rws(:,i)=((X-a).*(d-c))./(b-a)+c;             
        end
end
