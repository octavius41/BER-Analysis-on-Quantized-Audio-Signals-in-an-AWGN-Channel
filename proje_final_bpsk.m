close all; clc;

afp = "C:\Users\Grachus\OneDrive\Desktop\kom2\proje\Pulp Fiction - Apartment Scene Complete-[AudioTrimmer.com] (1).wav";
[y , fs] = audioread(afp);


% code modulates and starts to demodulate and decode immediately in 5
% seconds after that for demodulation of each SNR presented with bit error
% ratio in command window, performing each step approximately in 5 seconds

%% downsampling to decrease data points etc.
dsrate = 4;
yp = mean(y,2);
yp = reshape(yp,[1,length(yp)]);
t = (0:length(yp)-1) / fs;
N = length(yp);
w = linspace(-fs/2,fs/2,N);
fspec = fftshift(fft(yp,N));                            % frequency spectrum of the message signal in case we want to examine
ypds = zeros(1,round(length(yp)/dsrate));
tp = zeros(1,round(length(yp)/dsrate));
j = 1;
%assigning every one of 4 elements of the original signal in order to lower
%computational cost
for i = 1:dsrate:length(yp)
    tp(j) = t(i);
    ypds(j) = yp(i);                                    
    j =  j+1;

end

ypds = ypds/(max(abs(ypds)));                           % applying normalization to perform better distrubition to integer levels

%% quantization
n = 8;                                                  
% we performed shift scale operation in order to directly encode from integer
% levels, for this reason we choose n = 2^k-1, where k is the number of
% bits 
                                                        
ypds_qt = round(ypds * (n-1));              % rounded integer levels
ypds_qt = ypds_qt - min(ypds_qt);           % shifting minimum level to 0 rearranging integer levels as 0 - 15 for 4 bit encoding

qt_temp = (ypds_qt - mean(ypds_qt));           % re-configuring to initial levels 
qt_temp = qt_temp/max(abs(qt_temp));           % back normalization

q_err = sqrt(mean((ypds- ypds_qt).^2));        % quantization error

%%binary encoding

ypds_bi = dec2bin(ypds_qt, 4) - '0';        
ypds_bi = ypds_bi';                        % 4 bit binary encoding of 15 distinct levels 
ypds_bi = ypds_bi(:)';

%% bpsk mod.

bd = 0.001;fc = 10^5;A = 1;dt = 1/fs;bs = bd/dt;          % bit duration, carrier frequency, amplitude etc.
tc = 0:dt:(bd-dt);

fc1 = 10*10^3;
dp = pi/2;
md0 = A*cos(2*pi*fc1*tc);                                 % modulator signals
md1 = A*cos(2*pi*fc1*tc + pi/2);                          

bi_d = ypds_bi;          % a window taken out from binary encoded representation to keep computational cost low (don't care this)

bpsk=zeros(1,length(tc)*length(bi_d));            % pre assigning the space for modulated signal

%% modulation
for i = 1:length(bi_d)
    bpsk((1:bs)+bs*(i-1)) = A*cos(2*pi*fc1.*tc + bi_d(i)*pi/2);
end

%%

% GEÇERSİZ MODULASYON ŞEMASI HESAPLAMASI ÇOK UZUN
%{
for i = 1:length(bi_d);
    if bi_d(i) == 0
        bi_bpsk = [bi_bpsk md0];
    else
        bi_bpsk = [bi_bpsk md1];
    end
end
%}
%% awgn channel with different snrs' & demodulation and normalization

snr = [-30:30];
dcd_norm = zeros(length(snr),length(ypds_qt));
bers = zeros(1,length(snr));
for i = 1:length(snr)
    snr(i)
    ps = sum(abs(bi_d).^2)/length(bi_d);
    snrlin = 10^(0.1*snr(i));               %manually calculating the noise vector for given snr
    varn = ps / snrlin;     
    nt = sqrt(varn).*randn(1,length(bpsk));
    
    bi_nbpsk = nt+bpsk;                     % adding noise to the modulated signal
    
    bi_dmd = length(bi_d);
    
% correlator demodulation, window by window multiplication of whole
% received signal with the modulated signal to check similarity
    for j = 1:length(bi_d);
        th = xcorr(bi_nbpsk((1:bs)+bs*(j-1)),md1,0)/bs; 
        if th<0.42
            bi_dmd(j) = 0;
        else
            bi_dmd(j) = 1;
        end
    end
    
    % binary to decimal (digital to analog conversion)
    dcd = reshape(bi_dmd,[4,length(bi_dmd)/4])';                % 4 represents the 4 bit encoded levels as a known value
    
    dcd_out = [];
    for j=1:length(dcd)
        dcd_out(j) = bin2dec(num2str(dcd(j,:)));
    end
     
    temp = dcd_out - mean(dcd_out);
    dcd_norm(i,:) = temp/max(abs(temp));
    [m,n] = biterr(bi_dmd,bi_d);
    n
    bers(i) = n;
end
%%

bi_nbpsk = awgn(bpsk,-10);

figure(1)
plot(tp,ypds);title("downsampled and normalized signal"); xlabel("time");ylabel("amplitude") ; legend("ypds(t)");

figure(2)
plot(tp,ypds_qt);title("shift-scaled 16-level (4-bit) quantized signal"); xlabel("time");ylabel("amplitude") ; legend("ypds_qt(t)");

figure(3)
plot(bpsk);title("bpsk modulated signal"); xlabel("time");ylabel("amplitude") ; legend("bpsk(t)");xlim([26000 26500])

figure(4)
plot(bi_nbpsk);title("noisy bpsk signal, snr = -10"); xlabel("time");ylabel("amplitude") ; legend("nbfsk(t)");xlim([26000 26500])

figure(5)
plot(tp,qt_temp); title("4 bit quantized signal"); xlabel("time(s)");ylabel("amplitude") ; legend("quantized(t)");

figure(6)
plot(tp,dcd_norm(21,:)); title("demodulated part of the signal at snr = -10"); xlabel("time(s)");ylabel("amplitude") ; legend("dcd_norm(t)");

figure(7)
plot(tc,md0); title("modulator signal for bit 0"); xlabel("time(s)");ylabel("amplitude") ; legend("md0(t)");

figure(8)
plot(tc,md1); title("modulator signal for bit 1"); xlabel("time(s)");ylabel("amplitude") ; legend("md1(t)");

%% BER

load('bpsk_ber.mat');
load('bfsk_ber.mat');
snr = [-30:30];
EbNo = (-30:30)';
M = 2; % Modulation order

berPtheor = berawgn(EbNo,'psk',M,'nondiff');
berFtheor = berawgn(EbNo,'fsk',M,'coherent');

figure(9);
semilogy(snr,bpsk_ber);
hold on;
semilogy(snr,bfsk_ber);
hold on;
semilogy(snr,berFtheor);
hold on;
semilogy(snr,berPtheor);title("BER ANALYSIS"); xlabel("SNR"); ylabel("BER");legend("bpsk","bfsk","theoryF","theoryP");xlim([-30 12]);ylim([10^-5 1]);
grid on;
hold off;

%%


