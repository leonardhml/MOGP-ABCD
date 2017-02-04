function[]=simulquality()
%SIMULQUALITY  simulation of data for quality measures
%
%  SIMULQUALITY starts a window with menu for simulation of data for
%  quality measures
%
% (C) Jan Kolacek, Masaryk University (Czech Republic)

%if nargin<1
    
if ~exist('w.mat','file')
    inic=[];save w.mat inic;
end
    
callbackStr=['set(gcf,''CloseRequestFcn'',''closereq'');',...
      'clear;load puvprom;',...
      'delete(gcf);'];

crStr= 'save puvprom;';
% unts = get(0,'Units');
% set(0,'Units','Normalized');
% scrsz = get(0,'ScreenSize');
% set(0,'Units',unts);
% 
simfig=figure( ...
	'Visible','on', ...
   'Name','SIMULATION OF MODEL VARIABLES', ...
   'CloseRequestFcn',callbackStr,...
   'Units','Normalized',...
   'Tag','kclos',...
   'CreateFcn',crStr,...
   'NumberTitle','off');

axes('Units','normalized', ...
     'Position',[0.08 0.08 0.72 0.72]);

top=.85;
prom=char(who('-file','puvprom.mat'));
fsz=0.31;  %FontSize

% if size(prom,1)<1
%    rozh='off';
%    prom=char('[]','[]');
% else
%     rozh='on';
% end

%     call=['G=get(gcf,''userdata'');m_def=get(G(1),''String'');',...
%         'if ~isempty(m_def) set(G(6),''Enable'',''on'');else set(G(6),''Enable'',''off'');end;'];
call='G=get(gcf,''userdata'');G(1)=0;set(gcf,''userdata'',G);myrandtool;';
sim0Hndl=uicontrol( ...
        'Style','push', ...
        'Units','normalized', ...
        'FontUnits','normalized',...
        'Position',[.08,top+.03,.35,.07], ...
        'String','Generate a random sample for G0', ...
        'HorizontalAlignment','left',...
        'FontSize',fsz*1.3,...
        'Callback',call);


%     call=['G=get(gcf,''userdata'');m_def=get(G(1),''String'');',...
%         'if ~isempty(m_def) set(G(6),''Enable'',''on'');else set(G(6),''Enable'',''off'');end;'];
call='G=get(gcf,''userdata'');G(1)=1;set(gcf,''userdata'',G);myrandtool;';
sim1Hndl=uicontrol( ...
        'Style','push', ...
        'Units','normalized', ...
        'FontUnits','normalized',...
        'Position',[.45,top+.03,.35,.07], ...
        'String','Generate a random sample for G1', ...
        'HorizontalAlignment','left',...
        'FontSize',fsz*1.3,...
        'Callback',call);
    

% The SAVE text     
uicontrol( ...
	     'Style','text', ...
        'Units','normalized', ...
        'FontUnits','normalized',...
        'HorizontalAlignment','left',...
        'Position',[.84,.75,.2,.05], ...
        'BackgroundColor',[0.8 0.8 0.8], ...
        'FontSize',2*fsz,...
        'String','Save Data:');

% The VAR button
call=['checkLabels = {''Save design points for G0 to variable named:'','...
               '''Save design points for G1 to variable named:''};',... 
               'varNames = {''x_def'',''y_def''};',... 
'items = {x_def,y_def};save ppne x_def y_def; clear x_def y_def;',...
'save ppp;save pp checkLabels varNames items;clear;load pp;delete pp.mat;',...
'[hdialog,ok_pressed]=export2wsdlg(checkLabels,varNames,items,''Save Data'');'...
'if ok_pressed G=get(gcf,''userdata'');set(G(5),''enable'',''on'');',...
'clear G hdialog ok_pressed checkLabels varNames items;varNames=who;',...
'save w varNames -append;save(''puvprom'',varNames{:},''-append'');',...
'else load ppne;end;load ppp;delete ppp.mat;delete ppne.mat;'];
    labelStr='var';
    varHndl=uicontrol( ...
        'Style','push', ...
        'Enable','off',...
        'Units','normalized', ...
        'FontUnits','normalized',...
        'Position',[0.85,.7,.05,0.05], ...
        'String',labelStr, ...
        'FontSize',2*fsz,...
        'Callback',call);

% The FILE button
    labelStr='file';
    call='load w; uisave(varNames,''data_sim'')';
    fileHndl=uicontrol( ...
        'Style','push', ...
        'Enable','off',...
        'Units','normalized', ...
        'FontUnits','normalized',...
        'Position',[0.9,.7,.05,0.05], ...
        'String',labelStr, ...
        'FontSize',2*fsz,...
        'Callback',call);

% The CLOSE button
    labelStr='Close';
 uicontrol( ...
        'Style','push', ...
        'Units','normalized', ...
        'FontUnits','normalized',...
        'Position',[0.85,.07,.1,0.1], ...
        'String',labelStr, ...
        'FontSize',fsz,...
        'Callback',callbackStr);
     
osy=[0 1 0 1];
ktera=0;
% lgnd.x=0;
% lgnd.y=0;
hndlList=[ktera,sim0Hndl,sim1Hndl,varHndl,fileHndl,osy];
     set(gcf,'UserData',hndlList);    
     %set(C,'Position',[0 0.0347 1.0000 0.9190]);
     set(simfig,'Position',[0.1059 0.1655 0.7700 0.6898]);
     %set(C,'Position',[1 1 .97 .9].*scrsz+[.01 .012 0 0]);

% else
%     %rozhaz(fce,a,b,n,r);
% end