% Script for calling the gateway function AREBENCH.f, for generating
% the continuous and discrete-time AREs from the benchmark collection.
%

% RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
% Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
% D. Sima, Sept. 2001.
%
% Revisions:
% V. Sima, March 2009.

dico = input('Continuous-time example: 1 (default) or discrete-time example: 2   '); 
disp(' ')
if isempty(dico), dico = 1; end
while ((dico<1) || (dico>2))
        dico = input('Example type (1:2)  ');
end
    
disp('Select the group number:')
disp('1 - parameter-free problems of fixed size (default)')
disp('2 - parameter-dependent problems of fixed size')
disp('3 - examples of scalable size without parameters')
disp('4 - parameter-dependent examples of scalable size')
nr1 = input('Group number (1:4)  ');
disp(' ')
if isempty(nr1)
    nr1 = 1; 
end
while ((nr1<1) || (nr1>4))
        nr1 = input('Group number (1:4)  ');
end

if (dico==1)
    % CARE
    ex = [6 9 2 4];
    nr2 = input(['Select the example number (1:',num2str(ex(nr1)),')  ']);
    disp(' ')
    if isempty(nr2) 
        nr2 = 1; 
    end
    while (nr2<1) || (nr2>ex(nr1))
        nr2 = input(['Select the example number (1:',num2str(ex(nr1)),')  ']);
    end
    
    formG = input('Obtain G in unfactored form (1) or as B and R (0-default)?  ');
    disp(' ')
    if isempty(formG) || ((formG~=1) && (formG~=0)) 
        formG = 0; 
    end
    formQ = input('Obtain Q in unfactored form (1) or as C and W (0-default)?  ');
    disp(' ')
    if isempty(formQ) || ((formQ~=1) && (formQ~=0)) 
        formQ = 0; 
    end
   
    switch nr1
    case 1,
        if (formG==1)
            if (formQ==1)
                [X,A,G,Q] = arebench(1,1,nr2);
            else 
                [X,A,G,W,C] = arebench(1,1,nr2,[1 0]);
            end
        else 
            if (formQ==1)
                [X,A,R,Q,B] = arebench(1,1,nr2,[0 1]);
            else 
                [X,A,R,W,B,C] = arebench(1,1,nr2,[0 0]);
            end
        end
    case 2,
        epsilon = input('Epsilon: (press Enter to use the default value)'  );
        disp(' ')
        if isempty(epsilon)
            if (formG==1)
                if (formQ==1)
                    [X,A,G,Q] = arebench(1,2,nr2);
                else 
                    [X,A,G,W,C] = arebench(1,2,nr2,[1 0]);
                end
            else 
                if (formQ==1)
                    [X,A,R,Q,B] = arebench(1,2,nr2,[0 1]);
                else 
                    [X,A,R,W,B,C] = arebench(1,2,nr2,[0 0]);
                end
            end
        else if (formG==1)
                if (formQ==1)
                    [X,A,G,Q] = arebench(1,2,nr2,[1 1],epsilon);
                else 
                    [X,A,G,W,C] = arebench(1,2,nr2,[1 0],epsilon);
                end
            else 
                if (formQ==1)
                    [X,A,R,Q,B] = arebench(1,2,nr2,[0 1],epsilon);
                else 
                    [X,A,R,W,B,C] = arebench(1,2,nr2,[0 0],epsilon);
                end
            end
        end
    case 3,
        sizeN = input('The scalable size N: (press Enter to use the default value)  ');
        disp(' ')
        if isempty(sizeN)
            if (formG==1)
                if (formQ==1)
                    [X,A,G,Q] = arebench(1,3,nr2);
                else 
                    [X,A,G,W,C] = arebench(1,3,nr2,[1 0]);
                end
            else 
                if (formQ==1)
                    [X,A,R,Q,B] = arebench(1,3,nr2,[0 1]);
                else 
                    [X,A,R,W,B,C] = arebench(1,3,nr2,[0 0]);
                end
            end
        else 
            if (formG==1)
                if (formQ==1)
                    [X,A,G,Q] = arebench(1,3,nr2,[1 1],[],sizeN);
                else 
                    [X,A,G,W,C] = arebench(1,3,nr2,[1 0],[],sizeN);
                end
            else 
                if (formQ==1)
                    [X,A,R,Q,B] = arebench(1,3,nr2,[0 1],[],sizeN);
                else 
                    [X,A,R,W,B,C] = arebench(1,3,nr2,[0 0],[],sizeN);
                end
            end
        end
    otherwise, 
        switch nr2
        case 1,
            sizeN = input('The scalable size N: (press Enter to use the default value)  ');
            disp(' ')
            if ~isempty(sizeN)
                q = input('Parameter q: ');
                disp(' ')
                r = input('Parameter r: ');
                disp(' ')
                param = [q r];
            else param = [];    
            end
        case 2,
            sizeN = input('The scalable size N: (press Enter to use the default value)  ');
            disp(' ')
            if ~isempty(sizeN)
                a = input('Parameter a: ');
                disp(' ')
                b = input('Parameter b: ');
                disp(' ')
                c = input('Parameter c: ');
                disp(' ')
                beta1 = input('Parameter beta1: ');
                disp(' ')
                beta2 = input('Parameter beta2: ');
                disp(' ')
                gam1 = input('Parameter gamma1: ');
                disp(' ')
                gam2 = input('Parameter gamma2: ');
                disp(' ')
                param = [a b c beta1 beta2 gam1 gam2];
            else param = [];    
            end
        case 3,
            sizeN = input('The scalable size N: (press Enter to use the default value)  ');
            disp(' ')
            if ~isempty(sizeN)
                mu = input('Parameter mu: ');
                disp(' ')
                delta = input('Parameter delta: ');
                disp(' ')
                ka = input('Parameter kappa: ');
                disp(' ')
                param = [mu delta ka];
            else param = [];    
            end  
        otherwise,    
            sizeN = input('The scalable size N: (press Enter to use the default value)  ');
            disp(' ')
        end
        
        if isempty(sizeN)
            if (formG==1)
                if (formQ==1)
                    [X,A,G,Q] = arebench(1,4,nr2,[1 1]);
                else 
                    [X,A,G,W,C] = arebench(1,4,nr2,[1 0]);
                end
            else if (formQ==1)
                    [X,A,R,Q,B] = arebench(1,4,nr2,[0 1]);
                 else 
                    [X,A,R,W,B,C] = arebench(1,4,nr2,[0 0]);
                 end
            end
        else if (formG==1)
                 if (formQ==1)
                     [X,A,G,Q] = arebench(1,4,nr2,[1 1],param,sizeN);
                 else 
                     [X,A,G,W,C] = arebench(1,4,nr2,[1 0],param,sizeN);
                 end
             else if (formQ==1)
                     [X,A,R,Q,B] = arebench(1,4,nr2,[0 1],param,sizeN);
                 else 
                     [X,A,R,W,B,C] = arebench(1,4,nr2,[0 0],param,sizeN);
                 end
             end
        end
    end
%
else
    % DAREs
    ex = [13 5 0 1];
    if (nr1==3) 
        disp('Section empty')
        break
    end
    nr2 = input(['Select the example number (1:',num2str(ex(nr1)),')  ']);
    disp(' ')
    if isempty(nr2) 
        nr2 = 1;
    end
    while (nr2<1) || (nr2>ex(nr1))
        nr2 = input(['Select the example number (1:',num2str(ex(nr1)),')  ']);
    end
    
    formG = input('Obtain G in unfactored form (1) or as B and R (0-default)?  ');
    disp(' ')
    if isempty(formG) || ((formG~=1) && (formG~=0)) 
        formG = 0;
    end
    formQ = input('Obtain Q in unfactored form (1) or as C and W (0-default)?  ');
    disp(' ')
    if isempty(formQ) || ((formQ~=1) && (formQ~=0)) 
        formQ = 0;
    end
    formS = input('Return S (1=yes, 0=no (default))?'  );
    disp(' ')
    if isempty(formS) || ((formS~=1) && (formS~=0)) 
        formS = 0;
    end

    switch nr1
    case 1,
        if (formG==1)
            if (formQ==1)
                if (formS==1)
                    [X,A,G,Q,S] = arebench(2,1,nr2);
                else 
                    [X,A,G,Q] = arebench(2,1,nr2,[1 1 0]);
                end
            else
                if (formS==1)
                    [X,A,G,W,C,S] = arebench(2,1,nr2,[1 0 1]);
                else 
                    [X,A,G,W,C] = arebench(2,1,nr2,[1 0 0]);
                end
            end
        else 
            if (formQ==1)
                if (formS==1)
                    [X,A,G,Q,S] = arebench(2,1,nr2,[0 1 1]);
                else 
                    [X,A,G,Q] = arebench(2,1,nr2,[0 1 0]);
                end
            else
                if (formS==1)
                    [X,A,G,W,C,S] = arebench(2,1,nr2,[0 0 1]);
                else 
                    [X,A,G,W,C] = arebench(2,1,nr2,[0 0 0]);
                end
            end
        end
   case 2,
       if (nr2<5)
        delta = input('delta: (press Enter to use the default value)  ');
        disp(' ')
        if isempty(delta)
            if (formG==1)
                if (formQ==1)
                    if (formS==1)
                        [X,A,G,Q,S] = arebench(2,2,nr2);
                    else 
                        [X,A,G,Q] = arebench(2,2,nr2,[1 1 0]);
                    end
                else
                    if (formS==1)
                        [X,A,G,W,C,S] = arebench(2,2,nr2,[1 0 1]);
                    else 
                        [X,A,G,W,C] = arebench(2,2,nr2,[1 0 0]);
                    end
                end
            else 
                if (formQ==1)
                    if (formS==1)
                        [X,A,G,Q,S] = arebench(2,2,nr2,[0 1 1]);
                    else 
                        [X,A,G,Q] = arebench(2,2,nr2,[0 1 0]);
                    end
                else
                    if (formS==1)
                        [X,A,G,W,C,S] = arebench(2,2,nr2,[0 0 1]);
                    else 
                        [X,A,G,W,C] = arebench(2,2,nr2,[0 0 0]);
                    end
                end
            end
        else
            if (formG==1)
                if (formQ==1)
                    if (formS==1)
                        [X,A,G,Q,S] = arebench(2,2,nr2,[1 1 1],delta);
                    else 
                        [X,A,G,Q] = arebench(2,2,nr2,[1 1 0],delta);
                    end
                else
                    if (formS==1)
                        [X,A,G,W,C,S] = arebench(2,2,nr2,[1 0 1],delta);
                    else 
                        [X,A,G,W,C] = arebench(2,2,nr2,[1 0 0],delta);
                    end
                end
            else 
                if (formQ==1)
                    if (formS==1)
                        [X,A,G,Q,S] = arebench(2,2,nr2,[0 1 1],delta);
                    else 
                        [X,A,G,Q] = arebench(2,2,nr2,[0 1 0],delta);
                    end
                else
                    if (formS==1)
                        [X,A,G,W,C,S] = arebench(2,2,nr2,[0 0 1],delta);
                    else 
                        [X,A,G,W,C] = arebench(2,2,nr2,[0 0 0],delta);
                    end
                end
            end
        end
       else
        tau = input('Parameter tau: (press Enter to use the default value)  ');
        disp(' ')
        if isempty(tau)
            if (formG==1)
                if (formQ==1)
                    if (formS==1)
                        [X,A,G,Q,S] = arebench(2,2,nr2);
                    else 
                        [X,A,G,Q] = arebench(2,2,nr2,[1 1 0]);
                    end
                else
                    if (formS==1)
                        [X,A,G,W,C,S] = arebench(2,2,nr2,[1 0 1]);
                    else 
                        [X,A,G,W,C] = arebench(2,2,nr2,[1 0 0]);
                    end
                end
            else 
                if (formQ==1)
                    if (formS==1)
                        [X,A,G,Q,S] = arebench(2,2,nr2,[0 1 1]);
                    else 
                        [X,A,G,Q] = arebench(2,2,nr2,[0 1 0]);
                    end
                else
                    if (formS==1)
                        [X,A,G,W,C,S] = arebench(2,2,nr2,[0 0 1]);
                    else 
                        [X,A,G,W,C] = arebench(2,2,nr2,[0 0 0]);
                    end
                end
            end
        else
            D = input('Parameter D: ');
            disp(' ')
            K = input('Parameter K: ');
            disp(' ')
            r = input('Parameter r: ');
            disp(' ')
            param = [tau D K r];
            if (formG==1)
                if (formQ==1)
                    if (formS==1)
                        [X,A,G,Q,S] = arebench(2,2,nr2,[1 1 1],param);
                    else 
                        [X,A,G,Q] = arebench(2,2,nr2,[1 1 0],param);
                    end
                else
                    if (formS==1)
                        [X,A,G,W,C,S] = arebench(2,2,nr2,[1 0 1],param);
                    else 
                        [X,A,G,W,C] = arebench(2,2,nr2,[1 0 0],param);
                    end
                end
            else 
                if (formQ==1)
                    if (formS==1)
                        [X,A,G,Q,S] = arebench(2,2,nr2,[0 1 1],param);
                    else 
                        [X,A,G,Q] = arebench(2,2,nr2,[0 1 0],param);
                    end
                else
                    if (formS==1)
                        [X,A,G,W,C,S] = arebench(2,2,nr2,[0 0 1],param);
                    else 
                        [X,A,G,W,C] = arebench(2,2,nr2,[0 0 0],param);
                    end
                end
            end
        end
       end
    case 4,
        sizeN = input('The scalable size N: (press Enter to use the default value)  ');
        disp(' ')
        if isempty(sizeN)
            if (formG==1)
                if (formQ==1)
                    if (formS==1)
                        [X,A,G,Q,S] = arebench(2,4,nr2);
                    else 
                        [X,A,G,Q] = arebench(2,4,nr2,[1 1 0]);
                    end
                else
                    if (formS==1)
                        [X,A,G,W,C,S] = arebench(2,4,nr2,[1 0 1]);
                    else 
                        [X,A,G,W,C] = arebench(2,4,nr2,[1 0 0]);
                    end
                end
            else                     
                if (formQ==1)
                    if (formS==1)
                        [X,A,G,Q,S] = arebench(2,4,nr2,[0 1 1]);
                    else 
                        [X,A,G,Q] = arebench(2,4,nr2,[0 1 0]);
                    end
                else
                    if (formS==1)
                        [X,A,G,W,C,S] = arebench(2,4,nr2,[0 0 1]);
                    else 
                        [X,A,G,W,C] = arebench(2,4,nr2,[0 0 0]);
                    end
                end
            end
        else
            r = input('Parameter r: ');
            disp(' ')
            if (formG==1)
                if (formQ==1)
                    if (formS==1)
                        [X,A,G,Q,S] = arebench(2,4,nr2,[1 1 1],r,sizeN);
                    else 
                        [X,A,G,Q] = arebench(2,4,nr2,[1 1 0],r,sizeN);
                    end
                else
                    if (formS==1)
                        [X,A,G,W,C,S] = arebench(2,4,nr2,[1 0 1],r,sizeN);
                    else 
                        [X,A,G,W,C] = arebench(2,4,nr2,[1 0 0],r,sizeN);
                    end
                end
            else 
                if (formQ==1)
                    if (formS==1)
                        [X,A,G,Q,S] = arebench(2,4,nr2,[0 1 1],r,sizeN);
                    else 
                        [X,A,G,Q] = arebench(2,4,nr2,[0 1 0],r,sizeN);
                    end
                else
                    if (formS==1)
                        [X,A,G,W,C,S] = arebench(2,4,nr2,[0 0 1],r,sizeN);
                    else 
                        [X,A,G,W,C] = arebench(2,4,nr2,[0 0 0],r,sizeN);
                    end
                end
            end
        end
    end
end
