function [noiseseq, nofadein, nnonzero] = velvetnsmkx(len, density, delta, fs)
%VELVETNS Extended MK velvet noise sequence
% This function generates a velvet noise sequence 
% of length LEN samples, which contains DENSITY 
% ones and minus ones per second while all other 
% numbers in the sequence are zero. This sparse
% random noise sequence is useful in artificial
% reverberation algorithms.
%
% This version extends Matti's description by allowing less variation
% for the random jitter than average impulse spacing (Td). When delta = 1,
% you get Matti's original velvet noise sequence. When delta < 1, you get
% a sparse randomized sequence with less jitter. In the extreme case, 
% delta = 0, you get a pseudo-regular pulse sequence with random sample 
% values (-1 or +1). The pulse4 locations may still not be quite
% equidistant with delta = 0, because rounding of the pulse location may
% cause offsets. 
%
% Created by Vesa Valimaki, Olari, Espoo 18.5.2011
% Modified by Vesa Valimaki, Otaniemi, Espoo 29.1.2013
% Modified by Gloria Dal Santo, Otaniemi, Espoo 24.5.2022 to support
% arbitrary sampling rates

Td = fs/density;  % Average distance (in samples) between impulses 
noiseseq = zeros(1,len);  % Initialize velvet noise sequence
M = round(1.0001*len/Td);  % Estimated max number of pulses (0.01% extra)

% seed for replicabilty (seed = 13 for Arni)
seed = 18; 
rng(seed);
% Two random sequences are needed
rnd = rand(1,M);  % Random numbers between {0,1}
a = 2*round(rand(1,M))-1;  % Random sequence of -1's and +1's

for m = 1:M
    %k = round(Td*((m-1) + delta*rnd(m)));  % Index of next pulse (w/ delta)
    k = round(Td*(m-1) + (Td-1)*delta*rnd(m));%+1; % Next pulse location (w/ delta)
    if k > 0, noiseseq(k) = a(m); end;  % If added Jan. 24, 2013 (Vesa)
                                 % to avoid a rare error caused by k=0. 
end
noiseseq = noiseseq(1:len);  % Select the desired length
nnonzero = sum(abs(noiseseq));  % Count non-zero samples
nofadein = noiseseq;  % Save the original sequence

% Define envelope function
if length(noiseseq) >= 88
    env = [0.5*(1-cos(pi*(1:44)/44)) ones(1,len-88) 0.5*(1+cos(pi*(1:44)/44))];
    noiseseq = env .* noiseseq;  % Fade in and fade out to avoid clicks 
else
    fprintf('Cannot apply envelope function. Sequence too short \n');
end

end