close all; clc;

afp = "C:\Users\Grachus\OneDrive\Desktop\kom2\proje\Pulp Fiction - Apartment Scene Complete-[AudioTrimmer.com] (1).wav";
[y , fs] = audioread(afp);

%t_iv = 12500:16100;
%b_iv = 50001:64400;

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
fspec = fftshift(fft(yp,N));
ypds = zeros(1,round(length(yp)/dsrate));
tp = zeros(1,round(length(yp)/dsrate));
j = 1;
for i = 1:dsrate:length(yp)
    tp(j) = t(i);
    ypds(j) = yp(i);
    j =  j+1;

end

%% normalization
ypds = ypds/(max(abs(ypds)));
%% quantization
n = 4;



ypds_qt = round(ypds * (n-1));              % integers levels between -8 and +7
ypds_qt = ypds_qt - min(ypds_qt);           % shifting minimum level to 0 rearranging integer levels as 0 - 15

qt_temp = (ypds_qt - mean(ypds_qt));           % re-configuring to initial levels 
qt_temp = qt_temp/max(abs(qt_temp));

q_err = sqrt(mean((ypds- qt_temp).^2));     

%%binary encoding

ypds_bi = dec2bin(ypds_qt, 3) - '0';        
ypds_bi = ypds_bi';                        % 4 bit binary encoding of 15 distinct levels 
ypds_bi = ypds_bi(:)';

%% bfsk mod.

bd = 0.001;fc = 10^5;A = 1;dt = 1/fs;bs = bd/dt;          % bit duration, carrier frequency, amplitude etc.
tc = 0:dt:(bd-dt);

fc1 = 10*10^3;
df = 10^3;
fc2 = fc1+df;
md0 = A*cos(2*pi*fc1*tc);
md1 = A*cos(2*pi*fc2*tc);

bi_d = ypds_bi;                     % a window taken out from binary encoded representation to keep computational cost low

bfsk=zeros(1,length(tc)*length(bi_d));

for i = 1:length(bi_d)
    bfsk((1:bs)+bs*(i-1)) = A*cos(2*pi*(fc1 + bi_d(i)*df).*tc);
end


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
dcd_norm3 = zeros(length(snr),length(ypds_qt));
bers = zeros(1,length(snr));
for i = 1:length(snr)
    snr(i)
    ps = sum(abs(bi_d).^2)/length(bi_d);
    snrlin = 10^(0.1*snr(i));
    varn = ps / snrlin;     
    nt = sqrt(varn).*randn(1,length(bfsk));
    
    bi_nbfsk = nt+bfsk;
    
    bi_dmd = length(bi_d);
    
% correlator demodulation
    for j = 1:length(bi_d);
        th = xcorr(bi_nbfsk((1:bs)+bs*(j-1)),md1,0)/bs;
        if th<0.42
            bi_dmd(j) = 0;
        else
            bi_dmd(j) = 1;
        end
    end
    
    % binary to decimal (digital to analog conversion)
    dcd = reshape(bi_dmd,[3,length(bi_dmd)/3])';                % 4 represents the 4 bit encoded levels as a known value
    
    dcd_out = [];
    for j=1:length(dcd)
        dcd_out(j) = bin2dec(num2str(dcd(j,:)));
    end
        
    temp = dcd_out - mean(dcd_out);
    dcd_norm3(i,:) = temp/max(abs(temp));
    [m,n] = biterr(bi_dmd,bi_d);
    n
    bers(i) = n;
 
end
%%

%upload proje_final_bfsk_ws.mat and run section

figure(1)
plot(tp,ypds);title("downsampled and normalized signal"); xlabel("time");ylabel("amplitude") ; legend("ypds(t)");

figure(2)
plot(tp,ypds_qt);title("shift-scaled 16-level quantized signal"); xlabel("time");ylabel("amplitude") ; legend("ypds_qt(t)");

figure(3)
plot(bfsk);title("bfsk modulated signal"); xlabel("time");ylabel("amplitude") ; legend("bfsk(t)");xlim([26000 26500])
%%
figure(4)
plot(bi_nbfsk);title("noisy bfsk signal"); xlabel("time");ylabel("amplitude") ; legend("nbfsk(t)");xlim([26000 26500])
%%
figure(5)
plot(tp,qt_temp); title("window of quantized signal"); xlabel("time(s)");ylabel("amplitude") ; legend("qt_window(t)");

figure(6)
plot(tp,dcd_norm(41,:)); title("demodulated part of the signal at snr = 10"); xlabel("time(s)");ylabel("amplitude") ; legend("dcd_norm(t)");

figure(7)
plot(tc,md0); title("modulator signal for bit 0"); xlabel("time(s)");ylabel("amplitude") ; legend("md0(t)");

figure(8)
plot(tc,md1); title("modulator signal for bit 1"); xlabel("time(s)");ylabel("amplitude") ; legend("md1(t)");

%% BER

load('bpsk_ber_final.mat');                                                 
load('bfsk_ber_final.mat');
snr = [-30:30];
EbNo = (-30:30)';
M = 2; % Modulation order

berPtheor = berawgn(EbNo,'psk',M,'nodiff');
berFtheor = berawgn(EbNo,'fsk',M,'coherent');

figure(9);
semilogy(snr,bpsk_ber_up);
hold on;
semilogy(snr,bfsk_ber_up);
hold on;
semilogy(snr,berFtheor);
hold on;
semilogy(snr,berPtheor);title("BER ANALYSIS"); xlabel("SNR"); ylabel("BER");legend("bpsk","bfsk","theoryF","theoryP");xlim([-30 12]);ylim([10^-5 1]);
grid on;
hold off;
%%


bfsk_ber_up = bers;






