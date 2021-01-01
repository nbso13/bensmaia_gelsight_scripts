% stimuli as in muniak et al.
% amp : zero to peak sine amplitude
% frq : sine frq
% rad : probe radius (def:.5)
% pre_indent : pre-indentation (def: 0)

function s = stim_probe(amp,frq,rad,pre_indent,dur)

if(nargin<3)
    rad=.5;
end

if(nargin<4)
    pre_indent=0;
end

if(nargin<5)
    dur=max(.1,5/frq);
end

sf=20000;
t=(1/sf:1/sf:dur)';

s=Stimulus(pre_indent+amp*sin(2*pi*frq*t),[0 0],sf,rad);
