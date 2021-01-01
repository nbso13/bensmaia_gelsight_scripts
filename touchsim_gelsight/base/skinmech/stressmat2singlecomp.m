function [strain,names]=stressmat2singlecomp(stressmat,Comp)

E=.05;
nu = 0.4;

sigma_x=stressmat(:,:,1);
sigma_y=stressmat(:,:,2);
sigma_z=stressmat(:,:,3);
sigma_xy=stressmat(:,:,4);
sigma_zx=stressmat(:,:,5);
sigma_zy=stressmat(:,:,6);

StrainTrace=[];StressTrace=[];StrainDet=[];StressDet=[];
% strain matrix
if((Comp<9||Comp>15)&&Comp<19)
    % hooke's law
    eps_x  = (sigma_x - nu*(sigma_y+sigma_z))/E;
    eps_y  = (sigma_y - nu*(sigma_x+sigma_z))/E;
    eps_z  = (sigma_z - nu*(sigma_x+sigma_y))/E;
    eps_xy = (1+nu)/E*sigma_xy;
    eps_yz = (1+nu)/E*sigma_zy;
    eps_zx = (1+nu)/E*sigma_zx;
    
    if(Comp>2)
        [a,b,c]=eig33sym(eps_x,eps_xy,eps_zx,eps_y,eps_yz,eps_z);
        eigval=cat(3,a,b,c);
        PrincipalStrains = sort(eigval,3);
        StrainTrace = sum(PrincipalStrains,3) ;
        StrainDet   = prod(PrincipalStrains,3);
    end
else
    eps_x  = 0;    eps_y  = 0;    eps_z  = 0;
    eps_xy = 0;    eps_yz = 0;    eps_zx = 0;
    if(Comp>10)
        [a,b,c]=eig33sym(sigma_x,sigma_xy,sigma_zx,sigma_y,sigma_zy,sigma_z);
        eigval=cat(3,a,b,c);
        PrincipalStresses = sort(eigval,3);
        StressTrace = sum(PrincipalStresses,3);
        StressDet   = prod(PrincipalStresses,3);
    end
end

if(Comp==2), eps_h = abs(0.5*(eps_x+eps_y + sqrt( (eps_x-eps_y).^2 + 4*eps_xy.^2) ));else eps_h=0; end
if(Comp==10),sigma_h = abs(0.5*(sigma_x+sigma_y + sqrt( (sigma_x-sigma_y).^2 + 4*sigma_xy.^2) ));else sigma_h=[]; end
if(Comp==3||Comp==5), eps_c = abs(PrincipalStrains(:,:,1)); else eps_c=[];end
if(Comp==4||Comp==5), eps_t = abs(PrincipalStrains(:,:,3)); else eps_t=[];end
if(Comp==5), eps_m = max(eps_c,eps_t); else eps_m=[]; end
if(Comp==11||Comp==13), sigma_c = abs(PrincipalStresses(:,:,1)); else sigma_c=[];end
if(Comp==12||Comp==13), sigma_t = abs(PrincipalStresses(:,:,3)); else sigma_t=[];end
if(Comp==13), sigma_m = max(sigma_c,sigma_t); else sigma_m=[];end
if(Comp==16)
    StressTensor = cat(3,sigma_x, sigma_y, sigma_z, sigma_xy, sigma_zy, sigma_zx);
    StrainTensor = cat(3,eps_x  , eps_y  , eps_z  , eps_xy, eps_yz, eps_zx);
    SED = sum(StressTensor.*StrainTensor,3)/2;
else    SED=[]; end
if(Comp==8)
    a = abs(1-PrincipalStrains(:,:,1));
    b = abs(1-PrincipalStrains(:,:,2));
    c = abs(1-PrincipalStrains(:,:,3));
    p = 1.6;
    AreaChange = abs(( ((a.*b).^p + (b.*c).^p + (c.*a).^p)/3 ).^1/p - 1);
else    AreaChange = [];end
if(Comp==17),    VolumeChange = abs(prod(1-PrincipalStrains,3)-1);
else    VolumeChange = [];end
if(Comp==18),eps_s=.5*abs(PrincipalStrains(:,:,1)-PrincipalStrains(:,:,3));
else eps_s=[]; end
if(Comp==19),sigma_s=.5*abs(PrincipalStresses(:,:,1)-PrincipalStresses(:,:,3));
else sigma_s=[]; end

%           1     2     3     4     5       6           7         8        9       10      11      12      13        14         15      16       17
strain = {eps_z eps_h eps_c eps_t eps_m StrainTrace StrainDet AreaChange sigma_z sigma_h sigma_c sigma_t sigma_m StressTrace StressDet SED VolumeChange eps_s sigma_s};
names = {'eps_z','eps_h','eps_c','eps_t','eps_m','StrainTrace','StrainDet','AreaChange', 'sigma_z','sigma_h','sigma_c','sigma_t','sigma_m','StressTrace','StressDet','SED','VolumeChange','max_shear_strain','max_shear_stress'};

strain=strain{Comp};
names=names(Comp);

end

% function [ev1 ev2 ev3]=eig33sym(a11,a12,a13,a22,a23,a33)
%
%  Computes the eigenvalues of a 3x3 symetric matrix

function [ev1, ev2, ev3]=eig33sym(a11,a12,a13,a22,a23,a33)

c2 = - a11 - a22 - a33;
c1 = a11.*a22 + a11.*a33 + a22.*a33 - a12.^2 - a13.^2 - a23.^2;
c0 = a11.*a23.^2 + a22.*a13.^2 + a33.*a12.^2 - a11.*a22.*a33 - 2*a13.*a12.*a23;

p= c2.^2 - 3*c1;
q=-13.5*c0 - c2.^3 + 4.5*c1.*c2;

l=6.75*(c1.^2).*(p-c1)+27*c0.*(q+6.75*c0);
l(l<0)=0;

phi=atan2(sqrt(l),q)/3;

ev1=sqrt(p)/3.*(2*cos(phi))-c2/3;
ev2=sqrt(p)/3.*(2*cos(phi+2*pi/3))-c2/3;
ev3=sqrt(p)/3.*(2*cos(phi-2*pi/3))-c2/3;

end