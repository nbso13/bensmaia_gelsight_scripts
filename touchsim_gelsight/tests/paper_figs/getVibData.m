function [SINE,NOISE,DIHARM]=getVibData()

naff=17;
nrep=5;
ntype=[repmat({'SA'},4,1);repmat({'RA'},9,1);repmat({'PC'},4,1)];

%%
load mat/SINE.mat

%      SA         RA                         PC
nidx = [2 6 7 10   1 3 13 14 15 16 22 23 25   8 9 11 27];

stim_len = SINE.state_v(:,8);
idx_sine = 1:length(stim_len);
[~,stimidx]=sortrows(SINE.state_v(:,[5 2]));

amp=SINE.state_v(stimidx,2);
freq=SINE.state_v(stimidx,5);
len=SINE.state_v(stimidx,8);
stim=table(amp,freq,len);

spikes=cell(naff,1);
for ii=1:naff
    spikes{ii}=SINE.ObservedSpikes{nidx(ii)}(stimidx,:);
end

clear SINE

SINE.spikes=spikes;
SINE.stim=stim;
SINE.ncond=height(stim);
SINE.nrep=nrep;
SINE.naff=naff;
SINE.affclass=ntype;


%% NOISE
load mat/NOISE
chub=load('mat/pma_chubNoise1');
msha=load('mat/pma_msNoise1Atten0');
%      SA        RA                         PC
nidx = [2 4 5 8   1 3 11 12 13 14 18 19 20   6 7 9 21];


spikes=cell(naff,1);
for ii=1:naff
    spikes{ii}=NOISE.ObservedSpikes{nidx(ii)};
end

state_v=[msha.L2_str.state_v{1};chub.L2_str.state_v{1}];

amp=state_v(:,3);
freq1=state_v(:,4);
freq2=state_v(:,5);
len=state_v(:,1);

stim=table(amp,freq1,freq2,len);

clear NOISE

NOISE.spikes=spikes;
NOISE.stim=stim;
NOISE.ncond=height(stim);
NOISE.nrep=nrep;
NOISE.naff=naff;
NOISE.affclass=ntype;

%% DI HARM
load mat/DI
%      SA         RA                         PC
nidx = [2 6 7 10   1 3 13 14 15 16 21 22 24   8 9 11 27];
[~,stimidx]=sortrows(DI.state_v(:,[5 6 2]));

spikes=cell(naff,1);
for ii=1:naff
    spikes{ii}=DI.ObservedSpikes{nidx(ii)}(stimidx,:);
end

amp1=DI.state_v(stimidx,2);
amp2=DI.state_v(stimidx,3);
freq1=DI.state_v(stimidx,5);
freq2=DI.state_v(stimidx,6);
len=DI.state_v(stimidx,8);

stim=table(amp1,amp2,freq1,freq2,len);
clear DI

DIHARM.spikes=spikes;
DIHARM.stim=stim;
DIHARM.ncond=height(stim);
DIHARM.nrep=nrep;
DIHARM.naff=naff;
DIHARM.affclass=ntype;
