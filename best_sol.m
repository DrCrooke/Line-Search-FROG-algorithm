function [pulse, ref] = best_sol(obj, ref)
% code from https://doi.org/10.1364/OPTICA.3.001320
% It has not been modified in any way except this comment. / Christoffer O K
c1 = abs( xcorr(obj,ref) );
c2 = abs( xcorr(conj( obj(end:-1:1) ),ref) );
c3 = abs( xcorr(-obj,ref) );
c4 = abs( xcorr(conj( -obj(end:-1:1) ),ref)  );
[~, ind] = max([max(c1) max(c2) max(c3) max(c4)]);
switch ind
    case 1
    [~, indM] = max(abs(c1));
    pulseT = circshift( obj, numel(obj)-indM);
    case 2
    [~, indM] = max(abs(c2));
    pulseT = circshift( conj(obj(end:-1:1)), numel(obj)-indM);
    case 3
    [~, indM] = max(abs(c3));
    pulseT = circshift( -obj, numel(obj)-indM);
    case 4
    [~, indM] = max(abs(c4));
    pulseT = circshift( conj(-obj(end:-1:1)), numel(obj)-indM);
end
error = inf;
N = numel(obj);
F = ifftshift(-N/2:N/2-1);
for ii=1:200
    d = -1+ii/100;
    pulseR = ifft( fft(pulseT).*exp(1i*2*pi*d*F'/N) );
    if sqrt(sum(sum(abs( abs(pulseR)-abs(ref) ).^2)))<error
        error = sqrt(sum(sum(abs( abs(pulseR)-abs(ref) ).^2)));
        pulse = pulseR;
    end
end
end