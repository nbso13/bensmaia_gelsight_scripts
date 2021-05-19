function [P, Pdyn,S1] = HertzPressure(S0,xy,samp_freq,ProbeRad)

s=size(S0);

if(s(2)~=size(xy,1))
    disp(['Trace size : ' num2str(size(S0))])
    disp(['Pin number : ' num2str(size(xy,1))])
    error('incoherent number of pins ');
end

if(s(2)>1)
    error('Hertz Contact, single pins allowed');
end

E=0.050000; % 50kPa = 50,000 N/m^2 = 0.05 N/mm^2;
nu=0.4;

P=4/3*S0.^(3/2)*sqrt(ProbeRad)*E/(1-nu^2);

S1=S0;

S1p=([S1(2:end,:); nan(1,size(S1,2))] - [nan(1,size(S1,2)) ; S1(1:end-1,:)])/2*samp_freq;
S1p(1,:)=S1p(2,:); S1p(end,:)=S1p(end-1,:);

Pdyn=4/3*sign(S1p).*abs(S1p).^(3/2)*sqrt(ProbeRad)*E/(1-nu^2);