% compile MN_neuron_stream.m

sample_freq = 30;

maxlen = 1/sample_freq*5000*4;
codegen MN_neuron_stream_wrapper -o MN_neuron_stream_wrapper...
    -args {coder.typeof(double(0),[1 13]),coder.typeof(double(0),[4 1]),coder.typeof(double(0),[1 maxlen],[false false])}
