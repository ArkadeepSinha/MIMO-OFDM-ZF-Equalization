
close all;
clear all;
rng('shuffle');

Nsub = 256;
Ncp = round(Nsub/4);
nBlocks = 100;
L=2;
r=2;
t=2;
EbdB = [1:4:45];
Eb=10.^(EbdB/10);
No=1;
SNR=2*Eb/No;
SNRdB=10*log10(SNR);
BER = zeros(size(SNRdB));
BERt = zeros(size(SNRdB));

for bx = 1:nBlocks
    ChNoise = sqrt(No/2)*(randn(r,L+Nsub+Ncp-1)+ 1j*randn(r,L+Nsub+Ncp-1));
    H = 1/sqrt(2)*(randn(r,L,t) + 1j*randn(r,L,t));
    BitsI = randi([0,1], [t,Nsub]);
    BitsQ = randi([0,1], [t,Nsub]);
    Sym = (2*BitsI-1) + 1j*(2*BitsQ-1);
    Hfft = zeros(r,Nsub,t);
    for tx = 1:t
        Hfft(:,:,tx) = fft([H(:,:,tx),zeros(r,Nsub-L)],[],2);
    end
    for K= 1:length(SNRdB)
        LoadedSym = sqrt(Eb(K))*Sym;
        TxSamples = ifft(LoadedSym,[],2);
        TxSamCP = [TxSamples(:,Nsub-Ncp+1:Nsub),TxSamples];
        RxSamCP = zeros(r,L+Nsub+Ncp-1);
        for rx = 1:r
            for tx = 1:t
                RxSamCP(rx,:) = RxSamCP(rx,:)+conv(H(rx,:,tx),TxSamCP(tx,:));
            end
        end
        RxSamCP = RxSamCP + ChNoise;
        RxSamples = RxSamCP(:,Ncp+1:Ncp+Nsub);
        RxSym = fft(RxSamples,[],2);
        for nx = 1:Nsub
            Hsub = squeeze(Hfft(:,nx,:));
            ZFout = pinv(Hsub)*RxSym(:,nx);
            DecBitsI_ZF = (real(ZFout)>0);
            DecBitsQ_ZF = (imag(ZFout)>0);
            BER(K) = BER(K) + sum(sum(DecBitsI_ZF~=BitsI(:,nx)))+ sum(sum(DecBitsQ_ZF~=BitsQ(:,nx)));
        end
    end
end

BER = BER/nBlocks/Nsub/t/2;

semilogy(SNRdB,BER,'b -s ','linewidth',3.0,'MarkerFaceColor','b','MarkerSize',9.0);
hold on;
axis tight;
grid on;
legend('BER')
xlabel('SNR (dB)');
ylabel('BER');
title('BER vs SNR(dB) for MIMO-OFDM');

