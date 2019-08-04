clear all,
%I = load('I_OFDM1.txt').';
%Q = load('Q_OFDM1.txt').';
I = load('Real.txt').';
Q = load('Imag.txt').';
figure(1),subplot(2,1,1),plot(I),title('I data'),xlabel('Samples'),ylabel('Amplitude');
hold on, subplot(2,1,2),plot(Q),title('Q data'),xlabel('Samples'),ylabel('Amplitude');
Y = I + sqrt(-1)*Q;
fs = 20*10^6;
Ts = 1/fs; %duration
%sync by using ST
STLength = 20*8; %short training
J = zeros(1,length(Y)-STLength+1);
R = zeros(1,length(Y)-STLength+1);
M = zeros(1,length(Y)-STLength+1);
for l = 0:length(Y)-STLength
    r1 = Y(l+1:l+144);
    r2 = Y(l+17:l+160);
    r3 = Y(l+1:l+160);
    J(l+1) = sum(conj(r1).*(r2));
    R(l+1) = sum(abs(r3).^2);
    M(l+1)=(abs(J(l+1)).^2)/R(l+1);
end
[mx,idx] = max(abs(M));
Cor = abs(M);
figure(2),plot(Cor),title('Sync by correlation'),xlabel('Samples'),ylabel('Amplitude');%correlation output
%CFO by using LT
GILength = 20*1.6; %GI
LTLength = 20*(8-1.6); %long training
r4 = Y(idx+192+64:idx+192+64+63);
r5 = Y(idx+192:idx+192+63);
Ylong = sum((r4).*conj(r5));
CFO = angle(Ylong)*fs/(pi*128);
YCFO=exp(-1*sqrt(-1)*2*pi*CFO*Ts*(0:(length(Y)-1))) .*Y; %CFO Correction data
%channel estimation 
ylt1 = YCFO(idx+192:idx+192+63); %rec long training 1 seq
ylt2 = YCFO(idx+192+64:idx+192+64+63); %rec long training 2 seq
%FFT to freq domain
YLT1 = fft(ylt1);
YLT2 = fft(ylt2);
%removing padding
l_YLT1 = YLT1(2:27);
u_YLT1 = YLT1(39:64);
YLT1_Rm = [l_YLT1,u_YLT1];
l_YLT2 = YLT2(2:27);
u_YLT2 = YLT2(39:64);
YLT2_Rm = [l_YLT2,u_YLT2];
Xlt_minus = [1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,1,-1,1,-1,1,1,1,1];
Xlt_plus = [1,-1,-1,1,1,-1,1,-1,1,-1,-1,-1,-1,-1,1,1,-1,-1,1,-1,1,-1,1,1,1,1]; %long traing seq
Xlt = [Xlt_plus,Xlt_minus];
Hls = 0.5*1./Xlt.*(YLT1_Rm+YLT2_Rm);
%signal field
Y_Sig = YCFO(idx+160*2+16:idx+160*2+16+63);
FFT_YSig = fft(Y_Sig);
%removing padding
l_FFT_YSig = FFT_YSig(2:27);
u_FFT_YSig = FFT_YSig(39:64);
FFT_YSig_RmPd = [l_FFT_YSig,u_FFT_YSig];
%channel equalization
YSig_uncomp = FFT_YSig_RmPd./Hls;
%pilots (+-7,+-21)
Pilot_Seq = [1,-1,1,1];
idx_Pilots = [7,21,32,46];
Y_OFDM_comp = YSig_uncomp(idx_Pilots)./Pilot_Seq;
Y_OFDM_comp = Y_OFDM_comp./abs(Y_OFDM_comp);
YSig_comp_avg = mean(Y_OFDM_comp);
YSig_res_comp = YSig_uncomp*conj(YSig_comp_avg);
figure(3),scatter(real(YSig_res_comp),imag(YSig_res_comp)),title('Signal field constellation'),xlabel('I Amplitude'),ylabel('Q Amplitude'); 
axis([-1,1,-0.8,0.8]);%constellation
%demodulate signal field with BPSK
%remove pilots
YSig_res_comp(idx_Pilots)=[];
idx_bpsk1 = find(real(YSig_res_comp)>=0);
idx_bpsk0 = find(real(YSig_res_comp)<0);
YSig_res_comp(idx_bpsk1)=1;
YSig_res_comp(idx_bpsk0)=0;
YSig_res_comp = [YSig_res_comp(25:end),YSig_res_comp(1:24)];
%deinterleaving
N = length(YSig_res_comp);
s = 1;
j = 0:N-1;
i = s*floor(j/s)+mod((j+floor(16*j/N)),s);
k = 16*i-(N-1)*floor(16*i/N);%index
[B,I] = sort(k);
m = 1:N-1;
YSig_deinterlv = zeros(1,N);
YSig_deinterlv(m) = YSig_res_comp(I(m));
%viterbi decoding
trellis = poly2trellis(7,[133,171]);
YSig_bits = vitdec(YSig_deinterlv,trellis,24,'term','hard');
RATE = YSig_bits(1:4);%rate
LENGTH = YSig_bits(6:17);%length
Signal_Tail = YSig_bits(19:end);%tail 6'0'
if Signal_Tail ~= zeros(1,6)
   disp('signal field does not match');
elseif RATE == [1,1,0,1]
   disp('BPSK with 1/2 coding rate');
   disp(['OFDM length (octets): ',num2str(bi2de(LENGTH))]);
   ModSchm = 1;
   CodRate = 1/2;
   NBPSC = 1;
   NCBPS = 48;
   NDBPS = NCBPS*CodRate;
elseif RATE == [1,1,1,1]
   disp('BPSK with 3/4 coding rate');
   disp(['OFDM length (octets): ',num2str(bi2de(LENGTH))]); 
   ModSchm = 1;
   CodRate = 3/4;
   NBPSC = 1;
   NCBPS = 48;
   NDBPS = NCBPS*CodRate;
elseif RATE == [0,1,0,1]
   disp('QPSK with 1/2 coding rate');
   disp(['OFDM length (octets): ',num2str(bi2de(LENGTH))]); 
   ModSchm = 2;
   CodRate = 1/2;
   NBPSC = 2;
   NCBPS = 96;
   NDBPS = NCBPS*CodRate;
elseif RATE == [0,1,1,1]
   disp('QPSK with 3/4 coding rate');
   disp(['OFDM length (octets): ',num2str(bi2de(LENGTH))]);
   ModSchm = 2;
   CodRate = 3/4;
   NBPSC = 2;
   NCBPS = 96;
   NDBPS = NCBPS*CodRate;
elseif RATE == [1,0,0,1]
   disp('16QAM with 1/2 coding rate');
   disp(['OFDM length (octets): ',num2str(bi2de(LENGTH))]); 
   ModSchm = 3;
   CodRate = 1/2;
   NBPSC = 4;
   NCBPS = 192;
   NDBPS = NCBPS*CodRate;
elseif RATE == [1,0,1,1]
   disp('16QAM with 3/4 coding rate');
   disp(['OFDM length (octets): ',num2str(bi2de(LENGTH))]); 
   ModSchm = 3;
   CodRate = 3/4;
   NBPSC = 4;
   NCBPS = 192;
   NDBPS = NCBPS*CodRate;
elseif RATE == [0,0,0,1]
   disp('64QAM with 2/3 coding rate');
   disp(['OFDM length (octets): ',num2str(bi2de(LENGTH))]); 
   ModSchm = 4;
   CodRate = 2/3;
   NBPSC = 6;
   NCBPS = 288;
   NDBPS = NCBPS*CodRate;
elseif RATE == [0,0,1,1]
   disp('64QAM with 3/4 coding rate');
   disp(['OFDM length (octets): ',num2str(bi2de(LENGTH))]);
   ModSchm = 4;
   CodRate = 3/4;
   NBPSC = 6;
   NCBPS = 288;
   NDBPS = NCBPS*CodRate;
else
    disp('no match scheme');
end
%Data field
N_Oct = bi2de(LENGTH);%no of octets
N_PSDUbits = 8*N_Oct;%no of message bits
N_OFDM = ceil((N_PSDUbits+22)/NDBPS);%no of OFDM blocks
Y_OFDMData = YCFO(idx+160*2+80:idx+160*2+80+80*N_OFDM-1);
%fft 
Y_OFDMData = reshape(Y_OFDMData,[80,N_OFDM]);
%remove GI
Y_OFDMData(1:16,:) = [];
Y_OFDMDATA = fft((Y_OFDMData));
%remove padding
l_YOFDM = Y_OFDMDATA(2:27,:);
u_YOFDM = Y_OFDMDATA(39:64,:);
Y_OFDM_Rm = [l_YOFDM;u_YOFDM];
Y_OFDM_uncomp = Y_OFDM_Rm./Hls.';
figure(4),plot(Y_OFDM_uncomp(:),'o'),title('OFDM Data constellation'),xlabel('I Amplitude'),ylabel('Q Amplitude'); 
%pilots (+-7,+-21)
Pilot_SeqBase = [1,-1,1,1].';
%58 ofdm blocks correspond to 58 polarities
Pilot_polarity = [1,1,1,-1,-1,-1,1,-1,-1,-1,-1,1,1,-1,1,-1,-1,1,1,-1,1,1,-1,1,1,1,1,1,1,-1,1,1,1,-1,1,1,-1,-1,1,1,1,-1,1,-1,-1,-1,1,-1,1,-1,-1,1,-1,-1,1,1,1,1];
Pilot_OFDM = zeros(4,N_OFDM);
for i = 1:N_OFDM
    Pilot_OFDM(:,i) = Pilot_SeqBase.*Pilot_polarity(i);
end
idx_Pilots = [7,21,32,46];
%residual compensation
Y_OFDM_comp = Y_OFDM_uncomp(idx_Pilots,:)./Pilot_OFDM;
Y_OFDM_comp = Y_OFDM_comp./abs(Y_OFDM_comp);
Y_OFDM_comp_avg = mean(Y_OFDM_comp);
Y_OFDM_res_comp = Y_OFDM_uncomp.*repmat(conj(Y_OFDM_comp_avg),52,1);
%remove pilots
Y_OFDM_res_comp(idx_Pilots,:)=[];
Y_OFDM_res_comp = [Y_OFDM_res_comp(25:end,:);Y_OFDM_res_comp(1:24,:)];
figure(5),plot(Y_OFDM_res_comp(:),'o'),title('OFDM Data constellation'),xlabel('I Amplitude'),ylabel('Q Amplitude'); 
%demodulate (1:BPSK,2:QPSK,3:16QAM,4:64QAM)
switch ModSchm
    case 1 %BPSK
        idx_bpsk1 = find(real(Y_OFDM_res_comp)>=0);
        idx_bpsk0 = find(real(Y_OFDM_res_comp)<0);
        Y_OFDM_res_comp(idx_bpsk1)=1;
        Y_OFDM_res_comp(idx_bpsk0)=0;
    case 2 %QPSK
        p01 = [-1,1];
        p11 = [1,1];
        p00 = [-1,-1];
        p10 = [1,-1];
        p_QPSK = [p01;p11;p00;p10];
        bit_QPSK = [0,1;1,1;0,0;1,0];
        Y_OFDM_demo = [];
         for i = 1:N_OFDM %ofdm index
            for k = 1:48 %subcarrier index
                p_OFDM = [real(Y_OFDM_res_comp(k,i)),imag(Y_OFDM_res_comp(k,i))]./(1/sqrt(2));
                for j = 1:4
                    dist(j) = norm(p_OFDM-p_QPSK(j,:));
                end
                [mindist,minidx] = min(dist);
                Y_OFDM_demo1 = bit_QPSK(minidx,:);
                Y_OFDM_demo = [Y_OFDM_demo,Y_OFDM_demo1];
            end
        end
    case 3 %16QAM
        p_16QAM = [-3,3;-1,3;-3,1;-1,1;1,3;3,3;1,1;3,1;-3,-1;-1,-1;-3,-3;-1,-3;1,-1;3,-1;1,-3;3,-3];
        bit_16QAM = [0,0,1,0;0,1,1,0;0,0,1,1;0,1,1,1;1,1,1,0;1,0,1,0;1,1,1,1;1,0,1,1;0,0,0,1;0,1,0,1;0,0,0,0;0,1,0,0;1,1,0,1;1,0,0,1;1,1,0,0;1,0,0,0];
        Y_OFDM_demo = [];
         for i = 1:N_OFDM %ofdm index
            for k = 1:48 %subcarrier index
                p_OFDM = [real(Y_OFDM_res_comp(k,i)),imag(Y_OFDM_res_comp(k,i))]./(1/sqrt(10));
                for j = 1:16
                    dist(j) = norm(p_OFDM-p_16QAM(j,:));
                end
                [mindist,minidx] = min(dist);
                Y_OFDM_demo1 = bit_16QAM(minidx,:);
                Y_OFDM_demo = [Y_OFDM_demo,Y_OFDM_demo1];
            end
         end
    case 4 %64QAM
        p_64QAM = [-7,7;-5,7;-3,7;-1,7;-7,5;-5,5;-3,5;-1,5;-7,3;-5,3;-3,3;-1,3;-7,1;-5,1;-3,1;-1,1;1,7;3,7;5,7;7,7;1,5;3,5;5,5;7,5;1,3;3,3;5,3;7,3;1,1;3,1;5,1;7,1;-7,-1;-5,-1;-3,-1;-1,-1;-7,-3;-5,-3;-3,-3;-1,-3;-7,-5;-5,-5;-3,-5;-1,-5;-7,-7;-5,-7;-3,-7;-1,-7;1,-1;3,-1;5,-1;7,-1;1,-3;3,-3;5,-3;7,-3;1,-5;3,-5;5,-5;7,-5;1,-7;3,-7;5,-7;7,-7];
        bit_64QAMup = [0,0,0,1,0,0;0,0,1,1,0,0;0,1,1,1,0,0;0,1,0,1,0,0;0,0,0,1,0,1;0,0,1,1,0,1;0,1,1,1,0,1;0,1,0,1,0,1;0,0,0,1,1,1;0,0,1,1,1,1;0,1,1,1,1,1;0,1,0,1,1,1;0,0,0,1,1,0;0,0,1,1,1,0;0,1,1,1,1,0;0,1,0,1,1,0;1,1,0,1,0,0;1,1,1,1,0,0;1,0,1,1,0,0;1,0,0,1,0,0;1,1,0,1,0,1;1,1,1,1,0,1;1,0,1,1,0,1;1,0,0,1,0,1;1,1,0,1,1,1;1,1,1,1,1,1;1,0,1,1,1,1;1,0,0,1,1,1;1,1,0,1,1,0;1,1,1,1,1,0;1,0,1,1,1,0;1,0,0,1,1,0;];
        bit_64QAMdown = [0,0,0,0,1,0;0,0,1,0,1,0;0,1,1,0,1,0;0,1,0,0,1,0;0,0,0,0,1,1;0,0,1,0,1,1;0,1,1,0,1,1;0,1,0,0,1,1;0,0,0,0,0,1;0,0,1,0,0,1;0,1,1,0,0,1;0,1,0,0,0,1;0,0,0,0,0,0;0,0,1,0,0,0;0,1,1,0,0,0;0,1,0,0,0,0;1,1,0,0,1,0;1,1,1,0,1,0;1,0,1,0,1,0;1,0,0,0,1,0;1,1,0,0,1,1;1,1,1,0,1,1,;1,0,1,0,1,1;1,0,0,0,1,1;1,1,0,0,0,1;1,1,1,0,0,1;1,0,1,0,0,1;1,0,0,0,0,1;1,1,0,0,0,0;1,1,1,0,0,0;1,0,1,0,0,0;1,0,0,0,0,0];
        bit_64QAM = [bit_64QAMup;bit_64QAMdown];
        Y_OFDM_demo = [];
         for i = 1:N_OFDM %ofdm index
            for k = 1:48 %subcarrier index
                p_OFDM = [real(Y_OFDM_res_comp(k,i)),imag(Y_OFDM_res_comp(k,i))]./(1/sqrt(42));
                for j = 1:64
                    dist(j) = norm(p_OFDM-p_64QAM(j,:));
                end
                [mindist,minidx] = min(dist);
                Y_OFDM_demo1 = bit_64QAM(minidx,:);
                Y_OFDM_demo = [Y_OFDM_demo,Y_OFDM_demo1];
            end
         end
end
Y_OFDM_demo = reshape(Y_OFDM_demo,[NCBPS,N_OFDM]);
%deinterleaving
Y_OFDM_deinterlv = [];
s1 = max(NBPSC/2,1);
j = 0:NCBPS-1;
i = s1*floor(j/s1)+mod((j+floor(16*j/NCBPS)),s1);
k = 16*i-(NCBPS-1)*floor(16*i/NCBPS);%index
[B1,I1] = sort(k);
for a = 1:N_OFDM
    m = 1:NCBPS;
    Y_OFDM_deinterlv1(m) = Y_OFDM_demo(I1(m),a);
    Y_OFDM_deinterlv = [Y_OFDM_deinterlv;Y_OFDM_deinterlv1];
end
%viterbi decoding
Y_OFDM_vit =[];
TBlength = NDBPS;
switch CodRate %puncturing pattern
    case 1/2
        puncpat = [1;1];
    case 2/3
        puncpat = [1;1;1;0];
    case 3/4
        puncpat = [1;1;1;0;0;1];
end
trellis_data = poly2trellis(7,[133,171]);
for b = 1:N_OFDM
    Y_OFDM_bits1 = vitdec(Y_OFDM_deinterlv(b,:),trellis_data,TBlength,'term','hard',puncpat);
    Y_OFDM_vit = [Y_OFDM_vit,Y_OFDM_bits1];
end
%discrambling
%find where to start
initial_state = Y_OFDM_vit(1:7);
scrambler_pattern= [0,0,0,0,1,1,1,0,1,1,1,1,0,0,1,0,1,1,0,0,1,0,0,1,0,0,0,0,0,0,1,0,0,0,1,0,0,1,1,0,0,0,1,0,1,1,1,0,1,0,1,1,0,1,1,0,0,0,0,0,1,1,0,0,1,1,0,1,0,1,0,0,1,1,1,0,0,1,1,1,1,0,1,1,0,1,0,0,0,0,1,0,1,0,1,0,1,1,1,1,1,0,1,0,0,1,0,1,0,0,0,1,1,0,1,1,1,0,0,0,1,1,1,1,1,1,1];
diffSize = length(scrambler_pattern)-length(initial_state);
match = zeros(1,diffSize);
for i = 1:diffSize
    match(i) = all(scrambler_pattern(i:i-1+length(initial_state))==initial_state);
end
match_idx = find(match==1);
%reordered scrambler
if match_idx ~= 1
    scrambler = [scrambler_pattern(match_idx:end),scrambler_pattern(1:match_idx-1)];
else
    scrambler = scrambler_pattern;
end
Len_bits = length(Y_OFDM_vit);
N_scrambler = floor(Len_bits/127);
scrambler_bits = repmat(scrambler,1,N_scrambler);
Left_bits = Len_bits - length(scrambler_bits);
scrambler_bits = [scrambler_bits,scrambler_bits(1:Left_bits)];
OFDM_bits = xor(scrambler_bits,Y_OFDM_vit);
OFDM_bits = double(OFDM_bits);
%string form for later hex conversion
OFDM_bits_char = [];
for i = 1:length(OFDM_bits)
    OFDM_bits_char1 = num2str(OFDM_bits(i));
    OFDM_bits_char = [OFDM_bits_char OFDM_bits_char1];
end
%mac layer
%service bits
service_bits = OFDM_bits(1:16);
%frame control
frame_control1 = OFDM_bits_char(17:17+2*8-1);
TYPE = fliplr(frame_control1(3:4));
subTYPE = fliplr(frame_control1(5:8));
frame_control = bi2hex(frame_control1);
disp(['Frame control: ',frame_control]);
disp(['TYPE: ',TYPE]);
disp(['subTYPE: ',subTYPE]);
%duration
Duration1 = OFDM_bits_char(33:33+8*2-1);
Duration = bi2hex(Duration1);
disp(['Duration: ',Duration]);
%Receiver addr
RecAddr1 = OFDM_bits_char(49:49+6*8-1);
RecAddr = bi2hex(RecAddr1);
disp(['Receiver addr.: ',RecAddr]);
%Transmitter addr
TransAddr1 = OFDM_bits_char(97:97+6*8-1);
TransAddr = bi2hex(TransAddr1);
disp(['Transmitter addr.: ',TransAddr]);
%Destination addr
DestAddr1 = OFDM_bits_char(145:145+6*8-1);
DestAddr = bi2hex(DestAddr1);
disp(['Destination addr.: ',DestAddr]);
%Sequence control
SeqCtl1 = OFDM_bits_char(193:193+2*8-1);
SeqCtl = bi2hex(SeqCtl1);
disp(['Sequence control: ',SeqCtl]);
%=====binary to hex===============
function output = bi2hex(array)
    bit_length = length(array);
    Num_hex = bit_length/4;
    around_hex = floor(Num_hex);
    left_hex = Num_hex-around_hex;
    while left_hex > 0
        array = ['0' array];
        bit_length = length(array);
        Num_hex = bit_length/4;
        around_hex = floor(Num_hex);
        left_hex = Num_hex-around_hex;
    end
    Num_hex = length(array)/4;
    HEX = {'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F'};
    HexStr = [];
    for w = 1:Num_hex
        P = (w*4)-3;
        word = array(P:P+3);
        Dec = bin2dec(word);
        Hex = HEX{Dec+1};
        HexStr = [HexStr Hex];
    end
    Hex_length = length(HexStr);
    
    output = HexStr;
end



