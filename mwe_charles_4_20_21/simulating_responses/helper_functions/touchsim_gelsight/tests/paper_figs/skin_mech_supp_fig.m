% figures for skin mech suppl figure
clear all; close all; clc

%% HAND
h1=newfig('hand','unit','pixel','pos',[400 260 350 450]);
plot_hand(gca,'names',1,'axes',1,'centers',1);
set(gca,'pos',[0.01 0.01 .98 .98]); axis equal

screenshot(h1,'figs/skinmechsuppl_hand','pdf')

%% PIN
h2=figure;
set(h2,'pos',[680   743   560   235])
% deflection profile
[x,y]=meshgrid(-10:.25:10);
r=hypot(x,y); ProbeRad=2; E=.05; nu=.4;P=1;
D = (1-nu^2)/pi/ProbeRad*asin(ProbeRad./r)/E*P;
D(r<=ProbeRad)=(1-nu^2)/2/ProbeRad/E*P;
surfc(x,y,-D)
hold on
%cylinder
theta = (0:100)/100*2*pi;
sintheta = sin(theta); sintheta(101) = 0;
r=[2 2]';
x = r * cos(theta);
y = r * sintheta;
z = [-3 1]'*ones(1,101);
surf(x,y,z,zeros(size(z)))
% cylinder contour
x = 2 * cos(theta);
y = 2 * sintheta;
z = ones(1,101);
plot3(x,y,z,'k')
plot3([sqrt(2) sqrt(2)],-[sqrt(2) sqrt(2)],[-2 1],'k')
plot3(-[sqrt(2) sqrt(2)],[sqrt(2) sqrt(2)],[-2 1],'k')
xlabel('X [mm]');ylabel('Y [mm]');
% fig prop
shading interp;axis equal; view(-45,7.5);colormap([parula;0.7 0.7 0.7])

screenshot(h2,'figs/skinmechsuppl_skindeflection','pdf')

%% STRAIN AND WAVE PROFILES

f=200; sf=20000;
t=(1/sf:1/sf:.05)';

prad=.05;
rad=2;
[x,y]=meshgrid(-rad:.1:rad);
x=x(:); y=y(:); r=hypot(x,y);
x(r>rad)=[]; y(r>rad)=[]; r(r>rad)=[];

s0=.25+.25*sin(2*pi*f*t);
tic;stim1=Stimulus(s0,[0 0],sf,rad);toc
tic;stim2=Stimulus(repmat(s0,1,length(x)),[x y],sf,prad);toc

ap=AfferentPopulation;
xx=0:.0125:100;
for ii=1:length(xx)
    ap.add_afferents('RA',[xx(ii) 0],'depth',1);
end

tic
rc=ap.response(stim1);
r1=[rc.responses(:).propagated_struct];
toc
tic
rc=ap.response(stim2);
r2=[rc.responses(:).propagated_struct];
toc

s1=cat(2,r1(:).stat_comp);
u1=cat(2,r1(:).dyn_comp);

s2=cat(2,r2(:).stat_comp);
u2=cat(2,r2(:).dyn_comp);
%%
h3=figure;
cols=lines;
t=linspace(0,2*pi,200)';
set(h3,'pos',[680   280   400   700])
ax=subplot_ax(2,1,'nextplot','add');
plot(ax(1),xx,s1(end,:),xx,s2(end,:))
yl=ylim(ax(1));
plot(ax(1),[rad rad],yl,'k--','linew',1)
plot(ax(2),xx/10,u1(end,:),xx/10,u2(end,:))
yl=ylim(ax(2));
plot(ax(2),[rad rad]/10,yl,'k--','linew',1)
xlim(ax(1),[0 8]);xlim(ax(2),[0 8])
xlabel(ax(1),'Distance from probe center [mm]')
xlabel(ax(2),'Distance from probe center [cm]')
ylabel(ax(1),'Quasi-static')
ylabel(ax(2),'Dynamic')
lgd(ax(1),'Single probe','Probe grid','location','N')
set(ax(1:2),'box','off','ytick',0,'ygrid','on')

ax(3)=axes('pos',[.65 .62 .25 .25]);
ax(4)=axes('pos',[.65 .2 .25 .25]);
plot(ax(3),rad*cos(t),rad*sin(t),'color',cols(1,:)); axis equal
plot(ax(4),bsxfun(@plus,prad*cos(t),x'),bsxfun(@plus,prad*sin(t),y'),'color',cols(2,:),'linew',.5)
hold(ax(3),'on')
plot(ax(3),[0 rad*cosd(45)],[0 rad*sind(45)],'k')
text(.05*rad,.3*rad,'2 mm','par',ax(3),'rot',45)
hold(ax(4),'on')
plot(ax(4),[0 rad*cosd(45)],[0 rad*sind(45)],'k')
text(.05*rad,.3*rad,'2 mm','par',ax(4),'rot',45)
axis(ax(3),'equal'); axis(ax(4),'equal'); 
set(ax(3:4),'visible','off')

screenshot(h3,'figs/skinmechsuppl_singlevsmutli','pdf')

