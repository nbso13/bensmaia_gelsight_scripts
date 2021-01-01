function [spikes,state] = MN_neuron_stream_wrapper(parameters,state,trace)

displacement = reshape(trace,1,[]);

d = displacement;
vt = diff(d);
v = [vt 0];
at = diff(v);
a = [at 0];

d = resample(d,1,4);
v = resample(v,1,4);
a = resample(a,1,4);
dt = 1/5000;

% preparation of transduction current
vp =  v; vp(vp<0)=0;
vn = -v; vn(vn<0)=0;
dp =  d; dp(dp<0)=0;
dn = -d; dn(dn<0)=0;
ap =  a; ap(ap<0)=0;
an = -a; an(an<0)=0;

transduction_current = ...
    vp*parameters(5)/100 + vn*parameters(6)/100 + ...
    dp*parameters(7) + dn*parameters(8) + ...
    ap*parameters(9)/1000 + an*parameters(10)/1000;

I_0 = parameters(11)/1e9;

saturation = transduction_current;

tt = saturation(saturation>0);
saturation(saturation>0) = tt*I_0./(I_0+tt);
tt = saturation(saturation<0);
saturation(saturation<0) = tt*I_0./(I_0-tt);

input_current = saturation;

si_current = zeros(2,length(input_current));

[spikes,state] = ...
    MN_neuron_stream(dt, parameters(1), 150e-12, -.07, -.07, parameters(4), 10, -.03, -.03,...
    [.005; .05], parameters(2:3)'/1e9, 0, input_current, si_current, state);
end
