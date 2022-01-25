function callbk(action)
%Callback functions for the examples controller

%   RELEASE 2.0 of SLICOT Basic Systems and Control Toolbox.
%   Based on SLICOT RELEASE 5.7, Copyright (c) 2002-2020 NICONET e.V.
%
text1 = {...
         'Examples for testing the performance of'
         'Sylvester and Lyapunov solvers.'};
text2 = {...
         'Testing the performance of SLICOT '
         'Sylvester solver.'};
text3 = {...
         'Testing the performance of SLICOT '
         'continuous-time Lyapunov solver.'};
text4 = {...
         'Testing the performance of SLICOT '
         'discrete-time Lyapunov solver.'};
text5 = {... 
         'Testing the performance of SLICOT '
         'stable continous-time Lyapunov solver.'};   
text6 = {...
         'Testing the performance of SLICOT '
         'stable discrete-time Lyapunov solver.'};
text7 = {...
         'Examples for testing the performance of'
         'generalized Sylvester solvers.'};
text8 = {...
         'Testing the performance of SLICOT '
         'generalized  Sylvester solver.'};
text9 = {...
         'Testing the performance of SLICOT '
         'generalized continuous-time Lyapunov'
         'solver.'};
text10 = {...
         'Testing the performance of SLICOT '
         'generalized discrete-time Lyapunov'
         'solver.'};
text11 = {...
          'Testing the performance of SLICOT '
          'stable generalized continuous-time '
          'Lyapunov solver.'};
text12 = {...
          'Testing the performance of SLICOT '
          'stable generalized discrete-time '
          'Lyapunov solver.'};
texts = {text1 text2 text3 text4 text5 text6 text7...
         text8 text9 text10 text11 text12};
textinfo = {...
            'This demo was conceived to show how '
            'the Slicot routines for solving various '
            'linear or generalized matrix equations, '
            'work, when they are grouped in Matlab '
            'MEX-files, compared to the corresponding '
            'Matlab functions.'
            ' '
            'The test programs are written by the '
            'Slicot''s librarian, Vasile Sima, and'
            'the graphical part is due to Diana Sima.'
            };
            
switch action
case 'list'
   Val = get(gcbo,'Value');
   Str = get(gcbo,'String');
   beginBtn = findobj(gcbf,'String','Begin example');
   explic = findobj(gcbf,'Tag','explic');
   if Str{Val}(1)==' '
      set(beginBtn,'Enable','on');
   else
      set(beginBtn,'Enable','off');
   end;
   set(explic,'String',texts{Val});
case 'begin'
   List = findobj(gcbf,'Tag','Listbox');
   Val = get(List,'Value');
   Str = get(List,'String');
   if Str{Val}(1)==' '
      aux = Str{Val}(4:length(Str{Val}));
      eval(strcat('perf_',aux)) 
   end;
case 'next'
   uiresume
case 'info' 
   explic = findobj(gcbf,'Tag','explic');
   set(explic,'String',textinfo);
   beginBtn = findobj(gcbf,'String','Begin example');
   set(beginBtn,'Enable','off');

end;
