function trace = createProperTrace(x, F, D)
%   Standard fourier-domain delay of temporal signal.
N = length(D);
trace = zeros(length(F), N);

for i = 1:length(D)
    P = ifft( fft(x).*exp(1i*2*pi*D(i)*F) );
    trace(:, i) = fftshift( abs(fft(P.*x)/N) ).^2;
end
end