function [spikes,state] = MN_neuron_stream(dt,tau,C,Vr,V0,alpha,beta,Tr,T0,...
    sic_param_time_constant,sic_param_amplitude,arp,input_current,...
    si_current,state)
%#codegen
% MEX version of MN_neuron


TC = zeros(size(input_current));    % total current
SIC = si_current;       % spike induced current
V = NaN*TC;             % membrane potential
Theta = NaN*TC;         % threshold
Sp = NaN*TC;            % spikes
Time = NaN*TC;          % time


Time(1) = state(1);
V(1) = state(2);
Theta(1) = state(3);
Sp(1) = 0;
TC(1) = input_current(1);
for sici=1:size(SIC,1)
    TC(1) = TC(1)+SIC(sici,1);
end

latest_spike_time = state(4);
for ti=2:length(Time)
    
    Time(ti) = Time(ti-1) + dt;
    V(ti) =  V(ti-1)+( -(V(ti-1)-Vr)/tau + TC(ti-1)/C )*dt;
    Theta(ti) =  Theta(ti-1)+(  alpha*(V(ti-1)-Vr) - beta*(Theta(ti-1)-Tr)  )*dt;
    for sici=1:size(SIC,1)
        SIC(sici,ti) = SIC(sici,ti-1) - SIC(sici,ti-1)/sic_param_time_constant(sici)*dt;
    end
    Sp(ti) = 0;
    
    if arp/dt > (Time(ti)-latest_spike_time)/dt
        V(ti) = V(ti-1);
        Theta(ti) = Theta(ti-1);
    elseif V(ti)>=Theta(ti)
        Sp(ti) = 1;
        latest_spike_time = Time(ti);
        V(ti) = V0;
        Theta(ti) = max(Theta(ti),T0);
        for sici=1:size(SIC,1)
            SIC(sici,ti) = SIC(sici,ti)+sic_param_amplitude(sici);
        end
    end
    
    TC(ti) = input_current(ti);
    for sici=1:size(SIC,1)
        TC(ti) = TC(ti)+SIC(sici,ti);
    end
    
end

spikes = Time(Sp==1);
state = [Time(end); V(end); Theta(end); latest_spike_time];
