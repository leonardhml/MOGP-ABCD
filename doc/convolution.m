covmatrix;
mu = [0 0];
Sigma = [0.02 0.0; 0.0 0.4];
Sigma2 = [0.02 0.0; 0.0 0.02];

x1 = -3:.2:3; x2 = -3:.2:3;
[X1,X2] = meshgrid(x1,x2);
F1 = mvnpdf([X1(:) X2(:)],mu,Sigma);
F1 = reshape(F1,length(x2),length(x1));
F12 = mvnpdf([X1(:) X2(:)],mu,Sigma2);
F12 = reshape(F12,length(x2),length(x1));

%
cov=cri; hyp=hypri;
%
cov=cmi; hyp=hypmi;
%
cov=cgi; hyp=hypgi;
%
cov=cpe; hyp=hyppe;
%
cov = {'covSum',{cgi,cpe}}; hyp = [hypgi; hyppe];      % sum
%
cov = {@covProd,{cgi,cpe}};   hyp = [hypgi; hyppe];       % product

figure;

[F2,dK] = feval(cov{:},hyp,x1')
%

%F2 = mvnpdf([X1(:) X2(:)],mu,Sigma2);
%F2 = reshape(F2,length(x2),length(x1));

%F3=conv_fft2(F1,F2)
%F3=F3(16:46,16:46);
F4=conv2(F12,F2)
F4=F4(16:46,16:46);

%subplot(2,2,1)
%surf(x1,x2,F1)
%imagesc(F1)
%subplot(2,2,2)
%imagesc(F2)
%subplot(2,2,3)
%imagesc(F3)
subplot(2,2,4)
imagesc(F4)