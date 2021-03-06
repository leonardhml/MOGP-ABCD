function[]=okop(akce1);
%OKOP  calls up a demonstration window to show quality of estimation
%      regression function dependence on smoothing parameter
%
%     Use only in procedure KERN
%
% (C) Jan Kolacek, Masaryk University (Czech Republic)

if nargin<1
 akce1='inic';
end;

if strcmp(akce1,'inic')
CCfig=figure( ...
	'Visible','on', ...
	'Name','"Eye" method', ...
   'Tag','kclos',...
   'NumberTitle','off');
kr1=.1;kr2=.1;kr3=.1;h1=0;h2=0;h3=0;kon=0;
CC = get(CCfig,'Number');
if ischar(CC)
    CC = CCfig;
end

save v kr1 kr2 kr3 h1 h2 h3 CC kon -append;
axes( ...
        'Units','normalized', ...
        'Position',[0.08 0.08 0.75 0.87]);
load w;load v;
     plot(x,y,'x');%nx=min(x);mx=max(x);rr='b';
     if ff==1 hold on;plot(t,fx,rr);end;
     

% The CONSOLE frame
    uicontrol( ...
        'Style','frame', ...
        'Units','normalized', ...
        'FontUnits','normalized',...
        'Position',[0.86,0,.15,1], ...
        'BackgroundColor',[0.50 0.50 0.50]);


% The NW label
    
    uicontrol( ...
	     'Style','text', ...
        'Units','normalized', ...
        'Position',[.875,.95,.1,.05], ...
        'BackgroundColor',[0.50 0.50 0.50], ...
        'ForegroundColor','red', ...
        'FontSize',16,...
        'FontUnits','normalized',...
        'String','NW');

% The play button
    labelStr='->>';
    callbackStr1='load v;load w;set(gcbo,''Enable'',''off'');';
    callbackStr2=['hndlList=get(gcf,''UserData'');',...
       'hHndl1=hndlList(4);h1=get(hHndl1,''String'');h1=str2num(h1);',...
       'krHndl1=hndlList(5);kr1=get(krHndl1,''String'');kr1=str2num(kr1);'];
    callbackStr3=['while h1<1-kr1 h1=h1+kr1;stopHndl1=hndlList(2);',...
       'st=get(stopHndl1,''UserData'');set(stopHndl1,''Enable'',''on'');',...
       'if st==-1 break;end;nw(x,y,n,k,m,h1,CC);set(hHndl1,''String'',num2str(h1));',...
       'if ff==1 hold on;plot(t,fx,rr);hold off;end;pause(1);end;',...
       'set(stopHndl1,''UserData'',0,''Enable'',''off'');okop(''prac'');'];
 pl1=uicontrol( ...
        'Style','push', ...
        'Units','normalized', ...
        'FontUnits','normalized',...
        'Position',[0.87,.89,.06,0.05], ...
        'String',labelStr, ...
        'Callback',[callbackStr1,callbackStr2,callbackStr3]);

% The stop button
    labelStr='Stop';
    callbackStr1=['set(gcbo,''Userdata'',-1);',...
        'hndlList=get(gcf,''UserData'');pl1=hndlList(1);',...
        'set(pl1,''Enable'',''on'');'];
stopHndl1=uicontrol( ...
        'Style','push', ...
        'Units','normalized', ...
        'FontUnits','normalized',...
        'Position',[0.93,.89,.06,0.05], ...
        'String',labelStr, ...
        'UserData',0,...
        'Enable','off',...
        'Callback',[callbackStr1]);
% The <- button
    labelStr='<-';
    callbackStr1=[''''';load v;load w;nw(x,y,n,k,m,h1-kr1,CC);',...
          'hndlList=get(gcf,''UserData'');hHndl1=hndlList(4);',...
          'set(hHndl1,''String'',num2str(h1-kr1));pl1=hndlList(1);',...
          'set(pl1,''Enable'',''on'');okop(''prac'');'];
    callbackStr2='if ff==1 hold on;plot(t,fx,rr);hold off;end;';
    if h1==0 EnableStr='off';else EnableStr='on';end;
    leftHndl1=uicontrol( ...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.88,.83,.04,0.05], ...
        'String',labelStr, ...
        'FontUnits','normalized',...
        'Enable',EnableStr,...
        'Callback',[callbackStr1,callbackStr2]);
% The -> button
    labelStr='->';
    callbackStr1=[''''';load v;load w;nw(x,y,n,k,m,h1+kr1,CC);',...
          'hndlList=get(gcf,''UserData'');hHndl1=hndlList(4);',...
          'set(hHndl1,''String'',num2str(h1+kr1));okop(''prac'');'];
    callbackStr2='if ff==1 hold on;plot(t,fx,rr);hold off;end;';
    uicontrol( ...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.94,.83,.04,0.05], ...
        'FontUnits','normalized',...
        'String',labelStr, ...
        'Callback',[callbackStr1,callbackStr2]);
     
% Then the text label
    
    uicontrol( ...
	     'Style','text', ...
        'Units','normalized', ...
        'Position',[.875,.78,.04,.04], ...
        'BackgroundColor',[0.50 0.50 0.50], ...
        'FontUnits','normalized',...
	     'ForegroundColor',[1 1 1], ...
        'String','h =');
%pause;
% Then the editable text field
   hHndl1=uicontrol( ...
	'Style','edit', ...
        'Units','normalized', ...
        'Max',1, ...
        'Min',1, ...
        'String',num2str(h1),...%
        'BackgroundColor',[1 1 1], ...
 	     'Callback','okop(''prac'')' , ...
        'Position',[.91,.78,.08,.045], ...
        'HorizontalAlignment','left');
     
% Then the text label
    
    uicontrol( ...
	'Style','text', ...
        'Units','normalized', ...
        'Position',[.862,.71,.05,.05], ...
        'BackgroundColor',[0.50 0.50 0.50], ...
        'ForegroundColor',[1 1 1], ...
        'String','Step');
     
% Then the editable text field
   krHndl1=uicontrol( ...
	'Style','edit', ...
        'Units','normalized', ...
        'Max',1, ...
        'Min',1, ...
        'String',num2str(kr1),...%
        'BackgroundColor',[1 1 1], ...
 	     'Callback','okop(''prac'')' , ...
        'Position',[.91,.725,.08,.045], ...
        'HorizontalAlignment','left');
     
% The LL label
    
    uicontrol( ...
	'Style','text', ...
        'Units','normalized', ...
        'Position',[.875,.65,.1,.07], ...
        'BackgroundColor',[0.50 0.50 0.50], ...
        'ForegroundColor','red', ...
        'FontSize',16,...
        'String','LL');
    
% The play button
    labelStr='->>';
    callbackStr1='load v;load w;set(gcbo,''Enable'',''off'');';
    callbackStr2=['hndlList=get(gcf,''UserData'');',...
       'hHndl2=hndlList(9);h2=get(hHndl2,''String'');h2=str2num(h2);',...
       'krHndl2=hndlList(10);kr2=get(krHndl2,''String'');kr2=str2num(kr2);'];
 callbackStr3=['while h2<1-kr2 h2=h2+kr2;stopHndl2=hndlList(7);',...
  'st=get(stopHndl2,''UserData'');set(stopHndl2,''Enable'',''on'');',...
  'if st==-1 break;end;ll(x,y,n,k,m,h2,CC);set(hHndl2,''String'',num2str(h2));',...
  'if ff==1 hold on;plot(t,fx,rr);hold off;end;pause(1);end;',...
  'set(stopHndl2,''UserData'',0,''Enable'',''off'');okop(''prac'');'];

pl2=uicontrol( ...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.87,.6,.06,0.05], ...
        'String',labelStr, ...
        'Callback',[callbackStr1,callbackStr2,callbackStr3]);

% The stop button
    labelStr='Stop';
    callbackStr1=['set(gcbo,''Userdata'',-1);',...
          'hndlList=get(gcf,''UserData'');pl2=hndlList(6);',...
       'set(pl2,''Enable'',''on'');'];
stopHndl2=uicontrol( ...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.93,.6,.06,0.05], ...
        'String',labelStr, ...
        'UserData',0,...
        'Enable','off',...
        'Callback',[callbackStr1]);
% The <- button
    labelStr='<-';
    callbackStr1=[''''';load v;load w;ll(x,y,n,k,m,h2-kr2,CC);',...
          'hndlList=get(gcf,''UserData'');hHndl2=hndlList(9);',...
          'set(hHndl2,''String'',num2str(h2-kr2));pl2=hndlList(6);',...
       'set(pl2,''Enable'',''on'');okop(''prac'');'];
    callbackStr2='if ff==1 hold on;plot(t,fx,rr);hold off;end;';
    if h2==0 EnableStr='off';else EnableStr='on';end;
    leftHndl2=uicontrol( ...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.88,.54,.04,0.05], ...
        'String',labelStr, ...
        'Enable',EnableStr,...
        'Callback',[callbackStr1,callbackStr2]);
% The -> button
    labelStr='->';
    callbackStr1=[''''';load v;load w;ll(x,y,n,k,m,h2+kr2,CC);',...
          'hndlList=get(gcf,''UserData'');hHndl2=hndlList(9);',...
          'set(hHndl2,''String'',num2str(h2+kr2));okop(''prac'');'];
    callbackStr2='if ff==1 hold on;plot(t,fx,rr);hold off;end;';
    uicontrol( ...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.94,.54,.04,0.05], ...
        'String',labelStr, ...
        'Callback',[callbackStr1,callbackStr2]);
     
% Then the text label
    
    uicontrol( ...
	'Style','text', ...
        'Units','normalized', ...
        'Position',[.875,.485,.04,.04], ...
        'BackgroundColor',[0.50 0.50 0.50], ...
	'ForegroundColor',[1 1 1], ...
        'String','h =');
% Then the editable text field
   hHndl2=uicontrol( ...
	'Style','edit', ...
        'Units','normalized', ...
        'Max',1, ...
        'Min',1, ...
        'String',num2str(h2),...%
        'BackgroundColor',[1 1 1], ...
 	     'Callback','okop(''prac'')' , ...
        'Position',[.91,.49,.08,.045], ...
        'HorizontalAlignment','left');
     
% Then the text label
    
    uicontrol( ...
	'Style','text', ...
        'Units','normalized', ...
        'Position',[.862,.415,.05,.05], ...
        'BackgroundColor',[0.50 0.50 0.50], ...
        'ForegroundColor',[1 1 1], ...
        'String','Step');
% Then the editable text field
   krHndl2=uicontrol( ...
	'Style','edit', ...
        'Units','normalized', ...
        'Max',1, ...
        'Min',1, ...
        'String',num2str(kr2),...%
        'BackgroundColor',[1 1 1], ...
 	     'Callback','okop(''prac'')' , ...
        'Position',[.91,.43,.08,.045], ...
        'HorizontalAlignment','left');
     
 % The GM label
    
    uicontrol( ...
	'Style','text', ...
        'Units','normalized', ...
        'Position',[.875,.35,.1,.07], ...
        'BackgroundColor',[0.50 0.50 0.50], ...
        'ForegroundColor','red', ...
        'FontSize',16,...
        'String','GM');
    
% The play button
    labelStr='->>';
    callbackStr1='load v;load w;set(gcbo,''Enable'',''off'');';
    callbackStr2=['hndlList=get(gcf,''UserData'');',...
       'hHndl3=hndlList(14);h3=get(hHndl3,''String'');h3=str2num(h3);',...
       'krHndl3=hndlList(15);kr3=get(krHndl3,''String'');kr3=str2num(kr3);'];
 callbackStr3=['while h3<1-kr3 h3=h3+kr3;stopHndl3=hndlList(12);',...
  'st=get(stopHndl3,''UserData'');set(stopHndl3,''Enable'',''on'');',...
  'if st==-1 break;end;gm(x,y,n,k,m,h3,CC);set(hHndl3,''String'',num2str(h3));',...
  'if ff==1 hold on;plot(t,fx,rr);hold off;end;pause(1);end;',...
  'set(stopHndl3,''UserData'',0,''Enable'',''off'');okop(''prac'');'];

pl3=uicontrol( ...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.87,.3,.06,0.05], ...
        'String',labelStr, ...
        'Callback',[callbackStr1,callbackStr2,callbackStr3]);

% The stop button
    labelStr='Stop';
    callbackStr1=['set(gcbo,''Userdata'',-1);',...
          'hndlList=get(gcf,''UserData'');pl3=hndlList(11);',...
       'set(pl3,''Enable'',''on'');'];
stopHndl3=uicontrol( ...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.93,.3,.06,0.05], ...
        'String',labelStr, ...
        'UserData',0,...
        'Enable','off',...
        'Callback',[callbackStr1]);
% The <- button
    labelStr='<-';
    callbackStr1=[''''';load v;load w;gm(x,y,n,k,m,h3-kr3,CC);',...
          'hndlList=get(gcf,''UserData'');hHndl3=hndlList(14);',...
          'set(hHndl3,''String'',num2str(h3-kr3));pl3=hndlList(11);',...
       'set(pl3,''Enable'',''on'');okop(''prac'');'];
    callbackStr2='if ff==1 hold on;plot(t,fx,rr);hold off;end;';
    if h3==0 EnableStr='off';else EnableStr='on';end;
    leftHndl3=uicontrol( ...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.88,.24,.04,0.05], ...
        'String',labelStr, ...
        'Enable',EnableStr,...
        'Callback',[callbackStr1,callbackStr2]);
% The -> button
    labelStr='->';
    callbackStr1=[''''';load v;load w;gm(x,y,n,k,m,h3+kr3,CC);',...
          'hndlList=get(gcf,''UserData'');hHndl3=hndlList(14);',...
          'set(hHndl3,''String'',num2str(h3+kr3));okop(''prac'');'];
    callbackStr2='if ff==1 hold on;plot(t,fx,rr);hold off;end;';
    uicontrol( ...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.94,.24,.04,0.05], ...
        'String',labelStr, ...
        'Callback',[callbackStr1,callbackStr2]);
     
% Then the text label
    
    uicontrol( ...
	'Style','text', ...
        'Units','normalized', ...
        'Position',[.875,.185,.04,.04], ...
        'BackgroundColor',[0.50 0.50 0.50], ...
	'ForegroundColor',[1 1 1], ...
        'String','h =');
% Then the editable text field
   hHndl3=uicontrol( ...
	'Style','edit', ...
        'Units','normalized', ...
        'Max',1, ...
        'Min',1, ...
        'String',num2str(h3),...%
        'BackgroundColor',[1 1 1], ...
 	     'Callback','okop(''prac'')' , ...
        'Position',[.91,.19,.08,.045], ...
        'HorizontalAlignment','left');
     
% Then the text label
    
    uicontrol( ...
	'Style','text', ...
        'Units','normalized', ...
        'Position',[.862,.115,.05,.05], ...
        'BackgroundColor',[0.50 0.50 0.50], ...
        'ForegroundColor',[1 1 1], ...
        'String','Step');
% Then the editable text field
   krHndl3=uicontrol( ...
	'Style','edit', ...
        'Units','normalized', ...
        'Max',1, ...
        'Min',1, ...
        'String',num2str(kr3),...%
        'BackgroundColor',[1 1 1], ...
 	     'Callback','okop(''prac'')' , ...
        'Position',[.91,.13,.08,.045], ...
        'HorizontalAlignment','left');

 % The CLOSE button
    labelStr='Close';
    callbackStr='close;';
    closeHndl=uicontrol( ...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.88,.03,.1,0.07], ...
        'String',labelStr, ...
        'Callback',callbackStr);
     
     hndlList=[pl1 stopHndl1 leftHndl1 hHndl1 krHndl1,...
               pl2 stopHndl2 leftHndl2 hHndl2 krHndl2,...
               pl3 stopHndl3 leftHndl3 hHndl3 krHndl3];
    set(gcf,'Visible','on', ...
	'UserData',hndlList);
objects=get(gcf,'children');
set(objects,'FontUnits','normalized');
set(gcf,'Units','normalized','Position',[0.1059 0.1655 0.7700 0.6898]);

  
elseif strcmp(akce1,'prac')
hndlList=get(gcf,'UserData');   
hHndl1=hndlList(4);
hHndl2=hndlList(9);
hHndl3=hndlList(14);
krHndl1=hndlList(5);
krHndl2=hndlList(10);
krHndl3=hndlList(15);
leftHndl1=hndlList(3);
leftHndl2=hndlList(8);
leftHndl3=hndlList(13);
h1=get(hHndl1,'String');
kr1=get(krHndl1,'String');
h1=str2num(h1);
kr1=str2num(kr1);
h2=get(hHndl2,'String');
kr2=get(krHndl2,'String');
h2=str2num(h2);
kr2=str2num(kr2);
h3=get(hHndl3,'String');
kr3=get(krHndl3,'String');
h3=str2num(h3);
kr3=str2num(kr3);
l=length(h1);
if l>1
   h1=h1(l);
end;
if h1>0 
   set(leftHndl1,'Enable','on');
end;
if h2>0
   set(leftHndl2,'Enable','on');
end;
if h3>0
   set(leftHndl3,'Enable','on');
end;
save v.mat h1 h2 h3 kr1 kr2 kr3 -append;
end;
