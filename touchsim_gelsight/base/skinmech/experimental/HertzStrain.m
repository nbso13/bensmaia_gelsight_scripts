function [strain,names]=HertzStrain(P,PLoc,PRad,AffLoc,AffDepth,StrainComp)

% Cylindric
% Introduction to Contact Mechanics Second Edition
% Fischer-Cripps, Anthony C
% pp 88-89
%
% Louapre, Breder 2015
%

[nsamp,npin] = size(P); %(nsamp,npin)
nrec=size(AffLoc,1);

if(npin>1)
    error('this is Hertzian contact, only one probe should be used')
end

nu = 0.3;
E=.05;


x = bsxfun(@minus,AffLoc(:,1)',PLoc(:,1));    % (npin,nrec)
y = bsxfun(@minus,AffLoc(:,2)',PLoc(:,2));    % (npin,nrec)
z = bsxfun(@times,ones(npin,1),AffDepth(:)'); % (npin,nrec)

tt = atan2(y,x); % (npin,nrec)
s = sin(tt);     % (npin,nrec)
c = cos(tt);     % (npin,nrec)

r = hypot(x,y);

k=9/16*(1-nu^2);
a=(4/3*k*P*PRad/E).^(1/3)
sig0=3*P/2/pi/a.^2;

R=r/a;
Z=z/a;

L=sqrt(.5*(R.^2+Z.^2-1+sqrt((R.^2+Z.^2-1).^2+4*Z.^2)));
truc=Z./L.*(L*(1+nu).*atan(1./L)-(1-nu)*L.^2./(1+L.^2)-2*nu);

% stresses from Louapre&Breder, Int J Appl Ceram Technol (2015)
s_r=-L.^3.*R.^2.*Z./((L.^4+Z.^2).*(1+L.^2).^2)...
    -(1-2*nu).*(Z./(L.*(1+L.^2))-(1-(Z./L).^3)./(3*R.^2))+truc;
s_t=-(1-2*nu)./(3*R.^2).*(1-(Z./L).^3)+truc;
s_rz=-L.*R.*Z.^2./((L.^4+Z.^2).*(1+L.^2));
s_z=-Z.^3./(L.*(L.^4+Z.^2));

% on symmetry axis z==0
idx=Z==0&R<=1;
s_r(idx)=(1-2*nu)./(3*R(idx).^2).*(1-(1-R(idx).^2).^(3/2))-...
    (1-R(idx).^2).^(1/2);
s_t(idx)=-(1-2*nu)./(3*R(idx).^2).*(1-(1-R(idx).^2).^(3/2))-...
    2*nu*(1-R(idx).^2).^(1/2);
s_z(idx)=-(1-R(idx).^2).^(1/2);

idx=Z==0&R>1;
s_r(idx)=(1-2*nu)./(3*R(idx).^2);
s_t(idx)=-s_r(idx);
s_z(idx)=0;

% on symmetry axis r==0
s_r(r==0)=(1+nu)*Z(r==0).*atan(1./Z(r==0))+...
    (1-2*(1+nu)*(1+Z(r==0).^2))./(2*(1+Z(r==0).^2));
s_t(r==0)=s_r(r==0);
s_z(r==0)=-1./(1+Z(r==0).^2);

% Pressure rotated stress matrix (x,y,z)
s_x    = s_r.*c.^2 + s_t.*s.^2;
s_y    = s_r.*s.^2 + s_t.*c.^2;
s_xy   = (s_r-s_t).*s.*c;
s_xz   = s_rz.*c;
s_yz   = s_rz.*s;



sigma=zeros(nsamp,nrec,6);
sigma(:,:,1)    = sig0*s_x;  % (nsamp,nrec,6)
sigma(:,:,2)    = sig0*s_y;
sigma(:,:,3)    = sig0*s_z;
sigma(:,:,4)    = sig0*s_xy;
sigma(:,:,5)    = sig0*s_xz;
sigma(:,:,6)    = sig0*s_yz;

[strain,names]=stressmat2singlecomp(sigma,StrainComp);