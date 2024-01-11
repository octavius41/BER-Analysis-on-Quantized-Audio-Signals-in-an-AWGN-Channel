close all; clc;

% FOR SINGLE RUN WITH SPECIFIED PARAMETERS LIKE NUMBER OF BITS, SNR, MOD.
% TYPE ETC

afp ="C:\Users\Grachus\OneDrive\Desktop\kom2\proje\Pulp Fiction - Apartment Scene Complete-[AudioTrimmer.com].wav";
[y , fs] = audioread(afp);   % reading the audio file
bd = 0.001;fc = 10^5;A = 1;dt = 1/fs;bs = bd/dt;          % bit duration, carrier frequency, amplitude etc.
tc = 0:dt:(bd-dt);

%% downsampling to decrease data points etc. (mimicing a sampling operation asssuming that the first read of the audio file
dsrate = 4;                             % CHANGE DOWNSAMLPING RATE IF YOU ASK                                                             % is a continuos reading
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

figure(1)
plot(tp,ypds);title("downsampled and normalized signal"); xlabel("time");ylabel("amplitude") ; legend("ypds(t)");
%% quantization
n = 32;                                      %CHANGE QUANTIZATION LEVELS/2

ypds_qt = round(ypds * (n-1));              % integers levels between -8 and +7
ypds_qt_err = ypds_qt/(n-1);
ypds_qt = ypds_qt - min(ypds_qt);           % shifting minimum level to 0 rearranging integer levels as 0 - 15 
q_err = sqrt(mean((ypds - ypds_qt_err).^2));  

%%binary encoding

ypds_bi = dec2bin(ypds_qt, 6) - '0';          %CHANGE NUMBER OF BITS
ypds_bi = ypds_bi';                        % 4 bit binary encoding of 15 distinct levels 
ypds_bi = ypds_bi(:)';

binary_window = ypds_bi;
%%
figure(2)
plot(tp,ypds_qt);title("shift-scaled 8-level quantized signal"); xlabel("time");ylabel("amplitude") ; legend("ypds_qt(t)");xlim([0.02 0.03])



%% bpsk mod.
fc1 = 10^4;
df = 10^3;
fc2 = fc1+df;
dp = pi/2;
md0 = A*cos(2*pi*fc1*tc);
md1 = A*cos(2*pi*fc2*tc);                           %CHANGE MODULATOR SIGNAL FOR ASKED MODULATION TYPE

bfsk=zeros(1,length(tc)*length(binary_window));
bpsk=zeros(1,length(tc)*length(binary_window));

for i = 1:length(binary_window)
    bfsk((1:bs)+bs*(i-1)) = A*cos(2*pi*(fc1 + df*binary_window(i)).*tc);
end

for i = 1:length(binary_window)
    bpsk((1:bs)+bs*(i-1)) = A*cos(2*pi*fc1.*tc + binary_window(i)*pi/2);
end

%%

figure(3)
plot(bpsk);title("bpsk modulated signal"); xlabel("time");ylabel("amplitude") ; legend("bpsk(t)");xlim([26000 26500])

figure(4)
plot(bfsk);title("bfsk modulated signal"); xlabel("time");ylabel("amplitude") ; legend("bfsk(t)");xlim([26000 26500])


%%

%{
bfsk =  A*cos(2*pi.*fses.*tm);   
bpsk =  A*cos(2*pi*fc1.*tm + pses);
%}


ps = sum(abs(binary_window).^2)/length(binary_window);
snr = 10;                                               %CHANGE SNR VALUE
snrlin = 10^(0.1*snr);
varn = ps / snrlin;     
nt = sqrt(varn).*randn(1,length(bfsk));

bi_nbfsk = awgn(bfsk,10);


bi_nbpsk = nt+bpsk;

%%

figure(5)
plot(bi_nbpsk);title("noisy bfsk signal"); xlabel("time");ylabel("amplitude") ; legend("nbfsk(t)");xlim([26000 26500])

figure(6)
plot(bi_nbfsk);title("noisy bpsk signal"); xlabel("time");ylabel("amplitude") ; legend("nbpsk(t)");xlim([26000 26500])


%%
bi_dmd = zeros(1,length(binary_window));


    for j = 1:length(binary_window);
        th = xcorr(bi_nbfsk((1:bs)+bs*(j-1)),md1,0)/bs;         %CHANGE THE DEMODULATOR INPUT
        if th<0.42
            bi_dmd(j) = 0;
        else
            bi_dmd(j) = 1;
        end
    end

%figure(5)
%plot(tm,bpsk);
%{
%% correlator demodulation
for i = 1:length(binary_window)
    l0 = xcorr(bi_nbpsk((1:bs)+bs*(i-1)),A*cos(2*pi*fc1.*(tm(1:bs)+bs*(i-1))),0);
    l1 = xcorr(bi_nbpsk((1:bs)+bs*(i-1)),A*cos(2*pi*fc2.*(tm(1:bs)+bs*(i-1))),0);
    th = l0-l1
    if th<0
        bi_dmd(i) = 0;
    else
        bi_dmd(i) = 1;
    end
end
%}
[m,n] = biterr(binary_window,bi_dmd)
%%
qt_temp = (ypds_qt - mean(ypds_qt));           % re-configuring to initial levels 
qt_temp = qt_temp/max(abs(qt_temp));


    dcdr = reshape(bi_dmd,[6,length(bi_dmd)/6])';                % 4 represents the 4 bit encoded levels as a known value
                           % CHANGE NUMBER OF BITS
    dcd = [];
    for j=1:length(dcdr)
        dcd(j) = bin2dec(num2str(dcdr(j,:)));
    end

%%
    temp = dcd - mean(dcd);
    out = temp/max(abs(temp));

figure(7)
plot(tp,qt_temp); title("window of quantized signal"); xlabel("time(s)");ylabel("amplitude") ; legend("qt_window(t)");

figure(8)
plot(tp,out); title("demodulated part of the signal"); xlabel("time(s)");ylabel("amplitude") ; legend("out(t)");

%%

%%
figure(9)
plot(md0)

figure(10)
plot(md1)


%%

