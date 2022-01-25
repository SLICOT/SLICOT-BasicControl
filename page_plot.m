function page_plot(res,name,number)
%Makes the three subplots in the current 'page'

%   RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%   Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
%   Contributor:
%   D. Sima, Univ. of Bucharest, Romania, June 1999.
%
%   Revisions:
%   V. Sima,  28-05-2005, 04-03-2009, 06-01-2018.
%   
pause off

[m,n] = size(res);

set(gcf,'Name',strcat(name,' performances'),...
        'NumberTitle','off');

retu = findobj(gcf,'String','Return');
next = findobj(gcf,'String','Next>>');

if number == 1
     set(retu,'Enable','off');
else, set(retu,'Enable','on');
end

if number == 2
     set(next,'Enable','off');
else, set(next,'Enable','on');
end

switch n
case 7
 if strcmp(name,'slgesg')==0  
   sub = findobj(gcf,'Tag','Axes1');
   axes(sub)
   plot(res(:,1),res(:,2:3));
   legend('SLICOT','MATLAB','Location','Best')
   xlabel('Dimension n')
   ylabel('Time')
   set(sub,'Tag','Axes1')

   sub = findobj(gcf,'Tag','Axes2');
   axes(sub)
   plot(res(:,1),res(:,4:5));
   legend('SLICOT','MATLAB','Location','Best')
   xlabel('Dimension n')
   ylabel('Relative Residuals')
   set(sub,'Tag','Axes2')

   sub = findobj(gcf,'Tag','Axes3');
   axes(sub)
   plot(res(:,1),res(:,6:7));
   legend('SLICOT','MATLAB','Location','Best')
   xlabel('Dimension n')
   ylabel('Relative Errors')
   set(sub,'Tag','Axes3')
   
 else
   sub = findobj(gcf,'Tag','Axes1');
   axes(sub)
   stem3(res(:,1),res(:,2),res(:,3));
   xlabel('Dimension n')
   ylabel('Dimension m')
   zlabel('Time')
   set(sub,'Tag','Axes1')

   sub = findobj(gcf,'Tag','Axes2');
   axes(sub)
   bar(res(:,4:5));
   legend('in X','in Y','Location','Best')
   xlabel('Example #')
   ylabel('Relative Residuals')
   set(sub,'Tag','Axes2')

   sub = findobj(gcf,'Tag','Axes3');
   axes(sub)
   bar(res(:,6:7));
   legend('in X','in Y','Location','Best')
   xlabel('Example #')
   ylabel('Relative Errors')
   set(sub,'Tag','Axes3')
 end  
case 8
 if strcmp(name,'slgesg')==0 
   sub = findobj(gcf,'Tag','Axes1');
   axes(sub)
   bar(res(:,3:4));
   if m>10 
      set(sub,'XTick',[1 5 10 15]) 
   end
   legend('SLICOT','MATLAB','Location','Best')
   xlabel('Example #')
   ylabel('Time')
   set(sub,'Tag','Axes1')

   sub = findobj(gcf,'Tag','Axes2');
   axes(sub)
   bar(res(:,5:6));
   if m>10 
      set(sub,'XTick',[1 5 10 15]) 
   end
   legend('SLICOT','MATLAB','Location','Best')
   xlabel('Example #')
   ylabel('Relative Residuals')
   set(sub,'Tag','Axes2')

   sub = findobj(gcf,'Tag','Axes3');
   axes(sub)
   bar(res(:,7:8));
   if m>10 
      set(sub,'XTick',[1 5 10 15]) 
   end
   legend('SLICOT','MATLAB','Location','Best')
   xlabel('Example #')
   ylabel('Relative Errors')
   set(sub,'Tag','Axes3')
 else
   sub = findobj(gcf,'Tag','Axes1');
   axes(sub)
   stem3(res(:,1),res(:,2),res(:,3));
   xlabel('Dimension n')
   ylabel('Dimension m')
   zlabel('Time')
   set(sub,'Tag','Axes1')

   sub = findobj(gcf,'Tag','Axes2');
   axes(sub)
   bar(res(:,4:5));
   legend('in X','in Y','Location','Best')
   xlabel('Example #')
   ylabel('Relative Residuals')
   set(sub,'Tag','Axes2')

   sub = findobj(gcf,'Tag','Axes3');
   axes(sub)
   bar(res(:,6:7));
   legend('in X','in Y','Location','Best')
   xlabel('Example #')
   ylabel('Relative Errors')
   set(sub,'Tag','Axes3')
 end
case {4 , 5}
 if (strcmp(name,'slgely')==1) || ...
    (strcmp(name,'slgest')==1)
      
   sub = findobj(gcf,'Tag','Axes1');
   axes(sub)
   plot(res(:,1),res(:,2));
   %legend('SLICOT',2)
   xlabel('Dimension n')
   ylabel('Time')
   set(sub,'Tag','Axes1')

   sub = findobj(gcf,'Tag','Axes2');
   axes(sub)
   plot(res(:,1),res(:,3));
   %legend('SLICOT',2)
   xlabel('Dimension n')
   ylabel('Relative Residuals')
   set(sub,'Tag','Axes2')

   sub = findobj(gcf,'Tag','Axes3');
   axes(sub)
   plot(res(:,1),res(:,4));
   %legend('SLICOT',2)
   xlabel('Dimension n')
   ylabel('Relative Errors')
   set(sub,'Tag','Axes3')
   
 elseif (strcmp(name,'sllyap')==1) || ...
        (strcmp(name,'slstei')==1) || ...
        (strcmp(name,'slstly')==1) || ...
        (strcmp(name,'slstst')==1)
     
   sub = findobj(gcf,'Tag','Axes1');
   axes(sub)
   plot(res(:,1),res(:,2:3));
   legend('SLICOT','MATLAB','Location','Best')
   xlabel('Dimension n')
   ylabel('Time')
   set(sub,'Tag','Axes1')

   sub = findobj(gcf,'Tag','Axes2');
   axes(sub)
   plot(res(:,1),res(:,4:5));
   legend('SLICOT','MATLAB','Location','Best')
   xlabel('Dimension n')
   ylabel('Relative Residuals')
   set(sub,'Tag','Axes2')
 
   sub = findobj(gcf,'Tag','Axes3');
   chil = get(sub,'Children');
   legend(sub,'off')
   set(sub,'Tag','Axes3','Visible','off')
   set(chil,'Visible','off')

 else 
   
   sub = findobj(gcf,'Tag','Axes1');
   axes(sub)
   stem3(res(:,1),res(:,2),res(:,3));
   xlabel('Dimension n')
   ylabel('Dimension m')
   zlabel('Time')
   set(sub,'Tag','Axes1')

   sub = findobj(gcf,'Tag','Axes2');
   axes(sub)
   stem3(res(:,1),res(:,2),res(:,4));
   xlabel('Dimension n')
   ylabel('Dimension m')
   zlabel('Relative Residuals')
   set(sub,'Tag','Axes2')

   sub = findobj(gcf,'Tag','Axes3');
   chil = get(sub,'Children');
   legend(sub,'off')
   set(sub,'Tag','Axes3','Visible','off')
   set(chil,'Visible','off')

 end
case 3
         
   sub = findobj(gcf,'Tag','Axes1');
   axes(sub)
   plot(res(:,1),res(:,2));
   %legend('SLICOT',2)
   xlabel('Dimension n')
   ylabel('Time')
   set(sub,'Tag','Axes1')

   sub = findobj(gcf,'Tag','Axes2');
   axes(sub)
   plot(res(:,1),res(:,3));
   %legend('SLICOT',2)
   xlabel('Dimension n')
   ylabel('Relative Residuals')
   set(sub,'Tag','Axes2')

   sub = findobj(gcf,'Tag','Axes3');
   chil = get(sub,'Children');
   legend(sub,'off')
   set(sub,'Tag','Axes3','Visible','off')
   set(chil,'Visible','off')

end
   
   