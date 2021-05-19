function [strain,names]=HertzStrainFull(P,Q,PLoc,PRad,AffLoc,AffDepth,StrainComp)

% Cartesian
% Hamilton 1983

[nsamp,npin] = size(P); %(nsamp,npin)
nrec=size(AffLoc,1);

if(npin>1)
    error('this is Hertzian contact, only one probe should be used')
end

nu = 0.4;
E=.05;

x = bsxfun(@minus,AffLoc(:,1)',PLoc(:,1));    % (npin,nrec)
y = bsxfun(@minus,AffLoc(:,2)',PLoc(:,2));    % (npin,nrec)
z = bsxfun(@times,ones(npin,1),AffDepth(:)'); % (npin,nrec)

x2=x.^2; y2=y.^2; z2=z.^2; r2=x2+y2;
a=(3/4*(1-nu^2)*P*PRad/E).^(1/3);

A=r2+z2-a^2;
S=sqrt(A.^2+4.*a^2.*z2);
M=sqrt((S+A)/2);
N=sqrt((S-A)/2);

phi=atan(a./M);
G=M.^2-N.^2+z.*M-a*N;
H=2*M.*N+a*M+z.*N;

% normal load

nxx = (1+nu).*z.*phi+1./r2.*((y2-x2)./r2.*((1-nu).*N.*z2-(1-2.*nu)./3.*(N.*S+2.*A.*N+1)-nu.*M.*z)-N.*(x2+2.*nu.*y2)-M.*x2.*z./S);
nyy = (1+nu).*z.*phi+1./r2.*((x2-y2)./r2.*((1-nu).*N.*z2-(1-2.*nu)./3.*(N.*S+2.*A.*N+1)-nu.*M.*z)-N.*(y2+2.*nu.*x2)-M.*y2.*z./S);
nzz = -N+a*z.*M./S;
nxy = x.*y.*(1-2.*nu)./r2.^2.*(-N.*r2+2./3.*N.*(S+2.*A)-z.*(z.*N+M)+2./3)+x.*y.*z./r2.^2.*(-M.*r2./S-z.*N+M);
nyz = -z.*(y.*N./S-y.*z.*H./(G.^2+H.^2));
nxz = -z.*(x.*N./S-x.*z.*H./(G.^2+H.^2));

% tangential load

txx = -x.*(nu./4+1).*phi+x.*M./r2.^2.*((3./2-2.*x2./r2).*(S.*nu-2.*A.*nu+z2)+x2.*z2./S+7.*nu.*r2./4-2.*nu.*x2+r2)+x.*z.*N./r2.^2.*((3./2-2.*x2./r2).*(-S./6.*(1-2.*nu)-A./3.*(1-2.*nu)-1./2.*(z2+3))+x2./S-nu.*r2./4-7.*r2./4)+4.*x.*z./(3.*r2.^2).*(3./2-2.*x2./r2).*(1-2.*nu);
tyy = -3.*nu.*x.*phi./4+x.*M./r2.^2.*((1./2-2.*y2./r2).*(nu.*(S-2.*A+r2)+z2)+y2.*z2./S+3./4.*nu.*r2)+z.*x.*N./r2.^2.*((1./2-2.*y2./r2).*(-S./6.*(1-2.*nu)-A./3.*(1-2.*nu)-z2./2-3./2)+y2./S-3./4.*nu.*r2-r2./4)+4./3.*z.*x./r2.^2.*(1./2-2.*y2./r2).*(1-2.*nu);
tzz = z.*x.*N./(2.*r2).*(1-(r2+z2+1)./S);
txy = y./2.*(nu./2-1).*phi+y.*M./r2.^2.*(x2.*z2./S+nu.*((S-2.*A).*(1./2-2.*x2./r2)-2.*x2+r2./4)+r2./2+z2.*(1./2-2.*x2./r2))+y.*z.*N./r2.^2.*((1./2-2.*x2./r2).*((2.*nu-1).*(S./6+A./3)-z2./2-3./2-r2./2)+r2.*nu./4+x2./S-y2./2-3.*x2./2)+4.*y.*z./(3.*r2.^2).*(1./2-2.*x2./r2).*(1-2.*nu);
tyz = x.*y.*z./(2.*r2.^2).*(M.*(1./2+1./S.*(z2./2-3./2-r2./2))+z.*N./2.*(-3+1./S.*(5+z2+r2)));
txz = 3.*z.*phi./2+z.*M./r2.*(1+x2./r2-x2./S)+N./r2.*(-3./4.*(S+2.*A)+z2-3./4-1./4.*r2+z2./2.*(1./2-2.*x2./r2));

P0=3*P/2/pi/a.^3;
Q0=3*Q/2/pi/a.^3;

sigma=zeros(nsamp,nrec,6);
sigma(:,:,1)    = P0*nxx + Q0*txx;  % (nsamp,nrec,6)
sigma(:,:,2)    = P0*nyy + Q0*tyy;
sigma(:,:,3)    = P0*nzz + Q0*tzz;
sigma(:,:,4)    = P0*nxy + Q0*txy;
sigma(:,:,5)    = P0*nxz + Q0*txz;
sigma(:,:,6)    = P0*nyz + Q0*tyz;

[strain,names]=stressmat2singlecomp(sigma,StrainComp);