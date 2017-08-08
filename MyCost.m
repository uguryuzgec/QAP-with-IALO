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

function [z, p]=MyCost(s, model)

    n=model.n;
    w=model.w;
    d=model.d;

    [s, p] = sort(s);
    p = p(1:n);
    
    z=0;
    for i=1:n
        for j=i+1:n
            z=z+w(i,j)*d(p(i),p(j));
        end
    end

end