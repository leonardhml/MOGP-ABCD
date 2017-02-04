function HAZ_Fempdraw

mf=findobj('Tag','HAZ_MAIN');
udata=get(mf,'UserData');
X=udata.X;
xx=udata.xx;
%hist_bits=udata.hist_bits;
yy=distf_emp(X,xx);
figure(mf);
plot(xx,yy);
tit=title('Empirical distribution function');
set(tit,'FontUnits','Normalized');
set(tit,'FontSize',0.05);


