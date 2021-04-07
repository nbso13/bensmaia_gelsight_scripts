function [] = plotWelch(pxx, f)
%plotWelch plots the normalized power spectrum as passed.
plot(f, 10*log10(pxx))
ylabel("PSD (db/Hz)")
xlabel("Frequency (1/mm)")
end

