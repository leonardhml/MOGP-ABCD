function[]=HAZARD_setting()
% HAZARD_setting  - setting smoothing parameters
%
%     Use only in procedure ksdesn
%
% (C) Jiri Zelinka, Jan Kolacek, Masaryk University (Czech Republic)

   
HAZARD_sethndl=figure(...
   'Visible','on', ...
   'Units','normalized',...
   'Name','HAZARD_setting', ...
   'Color',[.3 .3 .3],...
   'Tag','HAZARD_setting',...
   'NumberTitle','off');

setting_data=struct('K',[],'h',[],'bhandles',[]);
mf=findobj('Tag','HAZARD_MAIN');
main_data=get(mf,'UserData');

K=main_data.K;
h=main_data.h;
xx=main_data.xx;

setting_data.K=K;
setting_data.h=h;
setting_data.X=main_data.X;
setting_data.d=main_data.d;

set(HAZARD_sethndl,'UserData',setting_data);

lb=min(xx);
ub=max(xx);
nx=length(xx);

%
% help
%

hlpstr=char(['First choose the kernel which will be used for kernel estimation. You can select a predefined or optimal kernel. ',...
	'For optimal kernel set the smoothness of the kernel. The kernel will be ',...
        'of order nu=0. The parameter mu is the ',...
        'smoothness of this kernel, mu=-1 is for non-continuous kernel. Finally set the bandwidth h ',...
        'using some of the given methods. You can also set the points for drawing the estimate.']);

% Then the editable text field
    eHndl=uicontrol( ...
	'Style','text', ...
        'Units','normalized', ...
        'Max',10, ...
        'BackgroundColor',[.9 .9 .9], ...
        'Position',[.05,.75,.9,.2], ...
        'FontSize',10, ...
        'FontUnits','normalized',...
        'HorizontalAlignment','left', ...
        'String',hlpstr);
% Then the text label

%
% kernel choice
%
    
% predefined kernels
    uicontrol( ...
        'Style','text', ...
        'Units','normalized', ...
        'FontUnits','normalized',...
        'Position',[.05,.68,.2,.05], ...
        'BackgroundColor',[0.30 0.30 0.30], ...
        'ForegroundColor',[0 0 0], ...
        'String','Predefined kernels:');

kchoice=uicontrol( ...
        'Style','popup', ...
        'Units','normalized', ...
        'FontUnits','normalized',...
        'Position',[0.03,.62,.2,0.05], ...
        'String',' Epanechnikov| Quartic| Gaussian| Rectangular| Triangular', ...
        'Tag','HAZARD_KPREDEF');

%  button OK
predef1='kdb=findobj(''Tag'',''Draw_k_button'');set(kdb,''Enable'',''on'');kch=findobj(''Tag'',''HAZARD_KPREDEF'');preval=get(kch,''Value'');';
predef2='switch preval, case 1,  K=K_def(''epan''); case 2, K=K_def(''quar''); case 3, K=K_def(''gaus''); case 4, K=K_def(''rect''); case 5, K=K_def(''tria'');end; ';
predef3='setfig=findobj(''Tag'',''HAZARD_setting'');udata=get(setfig,''UserData'');udata.K=K;set(setfig,''UserData'',udata);HAZARD_countbands;';
predef4='bh=udata.bhandles; nofbh=length(bh); for ii=1:nofbh, set(bh(1,ii),''Enable'',''on'');end, eyeh=findobj(''Tag'',''HAZARD_eyebut'');set(eyeh,''Enable'',''on'');bokh=findobj(''Tag'',''HAZARD_setbandOK'');set(bokh,''Enable'',''on'');';
predef5='ki=findobj(''Tag'',''HAZARD_kerinfo'');ks=get(kch,''String'');ks1=deblank(ks(preval,:));kis=[''Selected kernel: '',ks1,''.''];set(ki,''String'',kis);';
uicontrol( ...
        'Style','push', ...
        'Units','normalized', ...
        'FontUnits','normalized',...
        'Position',[0.25,0.62,.05,0.05], ...
        'String','OK', ...
        'Callback',[predef1,predef2,predef3,predef4,predef5]);

    uicontrol( ...
        'Style','text', ...
        'Units','normalized', ...
        'FontUnits','normalized',...
        'Position',[.05,.54,.2,.05], ...
        'BackgroundColor',[0.30 0.30 0.30], ...
        'ForegroundColor',[0 0 0], ...
        'String','Optimal kernel:');

% frame
    uicontrol( ...
        'Style','frame', ...
        'Units','normalized', ...
        'FontUnits','normalized',...
        'Position',[.03,0.47,.2,.08], ...
        'BackgroundColor',[0.50 0.50 0.50]);
% string
    uicontrol( ...
        'Style','text', ...
        'Units','normalized', ...
        'FontUnits','normalized',...
        'Position',[.05,.48,.08,.05], ...
        'BackgroundColor',[0.50 0.50 0.50], ...
        'ForegroundColor',[1 1 1], ...
        'String','Set nu');

% popup menu
    nHndl=uicontrol( ...
        'Style','popup', ...
        'Units','normalized', ...
        'FontUnits','normalized',...
        'Position',[.15,.495,.06,.04], ...
        'Tag','HAZARD_KOPT_nu',...
        'String',num2str((0)'));


% frame
    
    uicontrol( ...
        'Style','frame', ...
        'Units','normalized', ...
        'FontUnits','normalized',...
        'Position',[.03,0.39,.2,.08], ...
        'BackgroundColor',[0.50 0.50 0.50]);

% string

uicontrol( ...
   	'Style','text', ...
      'Units','normalized', ...
      'FontUnits','normalized',...
      'Position',[.05,.4,.08,.05], ...
      'BackgroundColor',[0.50 0.50 0.50], ...
      'ForegroundColor',[1 1 1], ...
      'String','Set k');
% popup menu
kHndl=uicontrol( ...
        'Style','popup', ...
        'Units','normalized', ...
        'FontUnits','normalized',...
        'Position',[.15,.415,.06,.04], ...
        'Tag','HAZARD_KOPT_k',...
        'String',num2str([2,4,6,8,]'));

% frame
    
    uicontrol( ...
        'Style','frame', ...
        'Units','normalized', ...
        'FontUnits','normalized',...
        'Position',[.03,0.31,.2,.08], ...
	'BackgroundColor',[0.50 0.50 0.50]);

% string
    uicontrol( ...
        'Style','text', ...
        'Units','normalized', ...
        'FontUnits','normalized',...
        'Position',[.05,.32,.08,.05], ...
        'BackgroundColor',[0.50 0.50 0.50], ...
        'ForegroundColor',[1 1 1], ...
        'String','Set mu');
% popup menu
    mHndl=uicontrol( ...
	'Style','popup', ...
        'Units','normalized', ...
        'FontUnits','normalized',...
        'Position',[.15,.335,.06,.04], ...
        'Tag','HAZARD_KOPT_mu',...
        'Value',2,...
        'String',num2str((-1:6)'));

%  button OK
optker1='kdb=findobj(''Tag'',''Draw_k_button'');set(kdb,''Enable'',''on'');KOPT_k=findobj(''Tag'',''HAZARD_KOPT_k'');k=2*get(KOPT_k,''Value'');';
optker2='KOPT_nu=findobj(''Tag'',''HAZARD_KOPT_nu'');nu=get(KOPT_nu,''Value'')-1;KOPT_mu=findobj(''Tag'',''HAZARD_KOPT_mu'');mu=get(KOPT_mu,''Value'')-2;';
optker3='if k<=nu, HAZARD_warning(''Parameters k has to be greater than nu. Correction was made, please check the parameters again.'');k=nu+2;set(KOPT_k,''Value'',k-1);end;';
optker4='K=K_def(''opt'',nu,k,mu);setfig=findobj(''Tag'',''HAZARD_setting'');udata=get(setfig,''UserData'');udata.K=K;set(setfig,''UserData'',udata);HAZARD_countbands;';
optker5='bh=udata.bhandles; nofbh=length(bh); for ii=1:nofbh, set(bh(1,ii),''Enable'',''on'');end, eyeh=findobj(''Tag'',''HAZARD_eyebut'');set(eyeh,''Enable'',''on'');bokh=findobj(''Tag'',''HAZARD_setbandOK'');set(bokh,''Enable'',''on'');';
optker6='ki=findobj(''Tag'',''HAZARD_kerinfo'');kis=[''Selected kernel: Optimal with nu='',num2str(nu),'', k='',num2str(k),'', mu='',num2str(mu),''.''];set(ki,''String'',kis);';
uicontrol( ...
        'Style','push', ...
        'Units','normalized', ...
        'FontUnits','normalized',...
        'Position',[0.25,0.4,.05,0.05], ...
        'String','OK', ...
        'Callback',[optker1,optker2,optker3,optker4,optker5,optker6]);

if isempty(K), enab_but='Off'; else, enab_but='On'; end
     DKH=uicontrol('Style','push',...
          'Units','normalized', ...
          'Enable',enab_but,...
          'Tag','Draw_k_button',...
          'FontUnits','normalized',...
          'Position',[.03,.18,.2,.08], ...
          'Callback','HAZARD_kerdraw', ...
          'BackgroundColor',[.5,0.5,0.5], ...
          'String','Draw the kernel');

 kselstr='Selected kernel:';
 K=setting_data.K;
 if ~isempty(K)
  switch K.type
   case 'tri', kselstr=[kselstr,' triangular'];
   case 'gau', kselstr=[kselstr,' gaussian'];
   case 'opt',
    k=K.k;mu=K.mu;nu=K.nu;
    if k==2 && nu==0
     if mu==-1, kselstr=[kselstr,' rectangular'];
     elseif mu==0, kselstr=[kselstr,' Epanechnikov'];
     elseif mu==1, kselstr=[kselstr,' Quartic'];
     else kselstr=[kselstr,' Optimal with nu=',num2str(nu),', k=',num2str(k),', mu=',num2str(mu),'.'];
     end
    else kselstr=[kselstr,' Optimal with nu=',num2str(nu),', k=',num2str(k),', mu=',num2str(mu),'.'];
    end
  end
 end
    kerinfohndl=uicontrol( ...
        'Style','text', ...
        'Units','normalized', ...
        'FontUnits','normalized',...
        'Tag','HAZARD_kerinfo',...
        'Position',[.03,.05,.25,.1], ...
        'BackgroundColor',[0.30 0.30 0.30], ...
        'ForegroundColor',[0 0 0], ...
        'String',kselstr, ...
        'HorizontalAlignment','left');

%
% bandwidth choice
%

    uicontrol( ...
        'Style','text', ...
        'Units','normalized', ...
        'FontUnits','normalized',...
        'Position',[.35,.68,.25,.05], ...
        'BackgroundColor',[0.30 0.30 0.30], ...
        'ForegroundColor',[0 0 0], ...
        'String','Choose the bandwidth:');

   banchoiHndl = uibuttongroup('visible','off','Position',[.35,.3,.275,.375],'Tag','HAZARD_BandChoi');
% bandwidth for iterative method
if isempty(K), enab_but='Off'; else, enab_but='On'; end
   Ithndl = uicontrol('Style','Radio',...
            'Enable',enab_but,...
	    'String','Iterative',...
	    'Units','normalized', ...
	    'Tag','HAZARD_BCIt', ...
	    'FontUnits','normalized',...
	    'FontSize',0.45,...
	    'Position',[.01,.82,.7,.16], ...
	    'parent',banchoiHndl,...
	    'HandleVisibility','on');
% bandwidth for cross-validation method
   CVhndl = uicontrol('Style','Radio',...
            'Enable',enab_but,...
	    'String','Cross-validation',...
	    'Units','normalized', ...
	    'Tag','HAZARD_BCCV', ...
	    'FontUnits','normalized',...
	    'FontSize',0.45,...
	    'Position',[.01,.57,.7,.16], ...
	    'parent',banchoiHndl,...
	    'HandleVisibility','on');
% bandwidth for maximal likelihood method
   MLhndl = uicontrol('Style','Radio',...
            'Enable',enab_but,...
	    'Tag','HAZARD_BCML', ...
	    'String','Max. likelihood',...
	    'Units','normalized', ...
	    'FontUnits','normalized',...
	    'FontSize',0.45,...
	    'Position',[.01,.32,.7,.16], ...
	    'parent',banchoiHndl,...
	    'HandleVisibility','on');
% bandwidth given manually
   Manhndl = uicontrol('Style','Radio',...
            'Enable',enab_but,...
	    'String','Set manualy',...
	    'Tag','HAZARD_BCMan', ...
	    'Units','normalized', ...
	    'FontUnits','normalized',...
	    'FontSize',0.45,...
	    'Position',[.01,.07,.7,.16], ...
	    'parent',banchoiHndl,...
	    'HandleVisibility','on');
% values of bandwidths for particular methods
   Ithndl2= uicontrol('Style','Text',...
	    'String',' ',...
	    'Units','normalized', ...
	    'Tag','HAZARD_BCIt_val', ...
	    'FontUnits','normalized',...
	    'FontSize',0.4,...
	    'Position',[.75,.80,.25,.15], ...
	    'parent',banchoiHndl,...
	    'HandleVisibility','on', ...
	    'HorizontalAlignment','left');
   CVhndl2= uicontrol('Style','Text',...
	    'String',' ',...
	    'Units','normalized', ...
	    'Tag','HAZARD_BCCV_val', ...
	    'FontUnits','normalized',...
	    'FontSize',0.4,...
	    'Position',[.75,.55,.25,.15], ...
	    'parent',banchoiHndl,...
	    'HandleVisibility','on', ...
	    'HorizontalAlignment','left');
   MLhndl2= uicontrol('Style','Text',...
	    'Tag','HAZARD_BCML_val', ...
	    'String',' ',...
	    'Units','normalized', ...
	    'FontUnits','normalized',...
	    'FontSize',0.4,...
	    'Position',[.75,.30,.25,.15], ...
	    'parent',banchoiHndl,...
	    'HandleVisibility','on', ...
	    'HorizontalAlignment','left');
   Manhndl2 = uicontrol('Style','Edit',...
	    'String',' ',...
	    'Tag','HAZARD_BCMan_val', ...
	    'Units','normalized', ...
	    'FontUnits','normalized',...
	    'FontSize',0.4,...
            'BackgroundColor',[0.90 0.90 0.90], ...
	    'Position',[.75,.050,.25,.15], ...
	    'parent',banchoiHndl,...
	    'HandleVisibility','on', ...
	    'HorizontalAlignment','left');

  bandhandles=[Ithndl, CVhndl, MLhndl, Manhndl;Ithndl2, CVhndl2, MLhndl2, Manhndl2];
  udata=get(HAZARD_sethndl,'UserData');
  udata.bhandles=bandhandles;
  set(HAZARD_sethndl,'UserData',udata);

  set(banchoiHndl,'SelectedObject',Ithndl);  % No selection
  set(banchoiHndl,'Visible','on');

  chband1='bc=findobj(''Tag'',''HAZARD_BandChoi'');bc0=get(bc,''SelectedObject'');bc1=findobj(''Tag'',''HAZARD_BCMan'');';
  chband2='if bc0==bc1, bcval=findobj(''Tag'',''HAZARD_BCMan_val'');h=eval(get(bcval,''String''));else, h=get(bc0,''UserData'');end,';
  chband3='sf=findobj(''Tag'',''HAZARD_setting'');udata=get(sf,''UserData'');udata.h=h;sb=findobj(''Tag'',''SelBand'');set(sb,''String'',num2str(h));set(sb,''UserData'',h);';
  chband4='mainok=findobj(''Tag'',''HAZARD_setMainOK'');set(mainok,''Enable'',''on'');';
  uicontrol( ...
        'Style','push', ...
        'Enable',enab_but,...
        'Tag','HAZARD_setbandOK', ...
        'Units','normalized', ...
        'FontUnits','normalized',...
        'Position',[0.65,0.4,.05,0.05], ...
        'String','OK', ...
        'Callback',[chband1,chband2,chband3,chband4]);

  eyehndl = uicontrol('Style','push','Units','normalized', ...
            'Enable',enab_but,...
            'Tag','HAZARD_eyebut', ...
            'FontUnits','normalized',...
            'Position',[.35,.18,.2,.08], ...
            'Callback','HAZARD_eye', ...
            'BackgroundColor',[.5,0.5,0.5], ...
            'String','"Eye" method');


  selband1= uicontrol( ...
        'Style','text', ...
        'Units','normalized', ...
        'FontUnits','normalized',...
        'Position',[.35,.075,.2,.075], ...
        'BackgroundColor',[0.30 0.30 0.30], ...
        'ForegroundColor',[0 0 0], ...
        'String','Selected bandwidth:', ...
	'HorizontalAlignment','left');

  if isempty(setting_data.h), selhstr='none';else selhstr=num2str(setting_data.h); end
  selband2= uicontrol( ...
        'Style','text', ...
        'Tag','SelBand', ...
        'Units','normalized', ...
        'UserData',setting_data.h, ...
        'FontUnits','normalized',...
        'Position',[.5,.1,.1,.05], ...
        'BackgroundColor',[0.30 0.30 0.30], ...
        'ForegroundColor',[0 0 0], ...
        'String',selhstr);

 if ~isempty(K)
  HAZARD_countbands;
 end

% selection of drawing points

    uicontrol( ...
        'Style','text', ...
        'Units','normalized', ...
        'FontUnits','normalized',...
        'Position',[.725,.68,.2,.05], ...
        'BackgroundColor',[0.30 0.30 0.30], ...
        'ForegroundColor',[0 0 0], ...
        'String','Points for drawing:');

% l-band
uicontrol( ...
        'Style','text', ...
        'Units','normalized', ...
        'FontUnits','normalized',...
        'FontSize',0.5,...
        'Position',[.725,.59,.14,.05], ...
        'HorizontalAlignment','left',...
        'BackgroundColor',[0.30 0.30 0.30], ...
        'ForegroundColor',[0 0 0], ...
        'String','Lower bound');

% edit field
uicontrol( ...
        'Style','edit', ...
        'Units','normalized', ...
        'FontUnits','normalized',...
        'Position',[.875,.60,.075,.05], ...
        'HorizontalAlignment','left',...
        'BackgroundColor',[0.9 0.9 0.9], ...
        'ForegroundColor',[0 0 0], ...
        'Tag','HAZARD_LBOUND',...
        'String',[' ',num2str(lb)]);

% u-band
uicontrol( ...
        'Style','text', ...
        'Units','normalized', ...
        'FontUnits','normalized',...
        'FontSize',0.5,...
        'Position',[.725,.53,.14,.05], ...
        'HorizontalAlignment','left',...
        'BackgroundColor',[0.30 0.30 0.30], ...
        'ForegroundColor',[0 0 0], ...
        'String','Upper bound');
% edit field
uicontrol( ...
        'Style','edit', ...
        'Units','normalized', ...
        'FontUnits','normalized',...
        'Position',[.875,.54,.075,.05], ...
        'HorizontalAlignment','left',...
        'BackgroundColor',[0.9 0.9 0.9], ...
        'ForegroundColor',[0 0 0], ...
        'Tag','HAZARD_UBOUND',...
        'String',[' ',num2str(ub)]);

% nx
uicontrol( ...
        'Style','text', ...
        'Units','normalized', ...
        'FontUnits','normalized',...
        'FontSize',0.5,...
        'Position',[.725,.47,.14,.05], ...
        'HorizontalAlignment','left',...
        'BackgroundColor',[0.30 0.30 0.30], ...
        'ForegroundColor',[0 0 0], ...
        'String','Num. of points');
% edit field
uicontrol( ...
        'Style','edit', ...
        'Units','normalized', ...
        'FontUnits','normalized',...
        'Position',[.875,.48,.075,.05], ...
        'HorizontalAlignment','left',...
        'BackgroundColor',[0.9 0.9 0.9], ...
        'ForegroundColor',[0 0 0], ...
        'Tag','HAZARD_NX',...
        'String',[' ',num2str(nx)]);

%%  button OK
% acc1='lb_hndl=findobj(''Tag'',''HAZARD_LBOUND'');lb=eval(get(lb_hndl,''String''));ub_hndl=findobj(''Tag'',''HAZARD_UBOUND'');ub=eval(get(ub_hndl,''String''));';
% acc2='nx_hndl=findobj(''Tag'',''HAZARD_NX'');nx=eval(get(nx_hndl,''String''));xx=linspace(lb,ub,nx);sf=findobj(''Tag'',''HAZARD_setting'');udata=get(sf,''UserData'');';
% acc3='udata.xx=xx;set(sf,''UserData'',udata);';
% uicontrol( ...
%        'Style','push', ...
%        'Units','normalized', ...
%        'FontUnits','normalized',...
%        'Position',[0.9,0.4,.05,0.05], ...
%        'String','OK', ...
%        'Callback',[acc1,acc2,acc3]);

     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% button Accept
  acc1='lb_hndl=findobj(''Tag'',''HAZARD_LBOUND'');lb=eval(get(lb_hndl,''String''));ub_hndl=findobj(''Tag'',''HAZARD_UBOUND'');ub=eval(get(ub_hndl,''String''));';
  acc2='nx_hndl=findobj(''Tag'',''HAZARD_NX'');nx=eval(get(nx_hndl,''String''));xx=linspace(lb,ub,nx);';
  acc3='mf=findobj(''Tag'',''HAZARD_MAIN'');sf=findobj(''Tag'',''HAZARD_setting'');selh=findobj(''Tag'',''SelBand'');mdata=get(mf,''UserData'');sdata=get(sf,''UserData'');hdata=get(selh,''UserData'');';
  acc4='LoH=mdata.LoH;nH=length(LoH);for ii=7:nH, set(LoH(ii),''Enable'',''On''); end; ';
  acc5='mdata.K=sdata.K;mdata.xx=xx;mdata.h=hdata;X=mdata.X;d=mdata.d;lam_est=K_hafest(X,d,K,xx,h);mdata.lam_est=lam_est;set(mf,''UserData'',mdata);set(gcf,''CloseRequestFcn'',''closereq'');delete(gcf);HAZARD_estimdraw;';
       uicontrol('Style','push','Units','normalized', ...
          'FontUnits','normalized',...
          'FontSize',0.3,...
          'Tag','HAZARD_setMainOK', ...
          'Enable',enab_but,...
          'Position',[.8,.15,.14,.08], ...
          'Callback',[acc1,acc2,acc3,acc4,acc5], ...
          'BackgroundColor',[.5,0.5,0.5], ...
          'String','Accept values');

% button Cancel
       closeStr='set(gcf,''CloseRequestFcn'',''closereq'');delete(gcf);';
       uicontrol('Style','push','Units','normalized', ...
          'FontUnits','normalized',...
          'Position',[.8,.02,.14,.08], ...
          'Callback',closeStr, ...
          'BackgroundColor',[.5,0.5,0.5], ...
          'String','Cancel');
       
set(gcf,'CloseRequestFcn',closeStr);
     
 
        
 
set(gcf,'Position',[0.1059 0.1655 0.7700 0.6898]);

% uloz data



