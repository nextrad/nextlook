clear; 
% Set these two correctly if you hope to read your file:
PULSES = 30;
RECORDING_BLOCK = 512;%1024;



IQ_SAMPLES = RECORDING_BLOCK*PULSES; %32 bits per sample
BUF = IQ_SAMPLES; 
LEN = IQ_SAMPLES*2;

%fida=fopen('C:\ReadyFlow\71620\Windows\4.2\win64\examples\ddc_multichan\adc0data.dat','r');
%fidb=fopen('C:\ReadyFlow\71620\Windows\4.2\win64\examples\ddc_multichan\adc1data.dat','r');
%fidc=fopen('C:\ReadyFlow\71620\Windows\4.2\win64\examples\ddc_multichan\adc2data.dat','r');

fida=fopen('/media/darryn/FLASH/2017_11_17_15_22_13_adc0data.dat','r');
fidb=fopen('../adc1data.dat','r');
fidc=fopen('../adc2data.dat','r');

% Read first part

aa=fread(fida,LEN,'int16');
a=aa(1:2:LEN)+j*aa(2:2:LEN);


%aa=fread(fida,LEN,'int16');
%a=aa(1:2:LEN)+j*aa(2:2:LEN);


%bb=fread(fidb,LEN,'int16');
%b=bb(1:2:LEN)+j*bb(2:2:LEN);

%cc=fread(fidc,LEN,'int16');
%c=cc(1:2:LEN)+j*cc(2:2:LEN);


%Plot first five
if 1
hold off;
plot(real(a(1:BUF)),'r');
hold on;
plot(imag(a(1:BUF)),'b');

figure

hold off;
plot(real(a(1:BUF)),'r');
hold on;
plot(imag(a(1:BUF)),'b');

%plot(real(a(1:1000)));
%plot(unwrap(angle(a(1:500))));
%hold on;
%plot(real(b(1:BUF)),'r');
%plot(real(c(1:BUF)),'k');
%plot(unwrap(angle(b(1:500))),'r');
%plot(real(b(1:1000)),'r');
end


if 0
hold off;
figure
plot(real(a));
%plot(real(a(1:1000)));
%plot(unwrap(angle(a(1:500))));
hold on;
plot(real(b),'r');
plot(real(c),'k');
%plot(unwrap(angle(b(1:500))),'r');
%plot(real(b(1:1000)),'r');
end


if 0
	
flen=length(a);

vert = 20*log10(flen)+20*log10(2^16)-21; % adjust vertical scale for any length FFT and 16 bit data

    vert = 167;
    figure
    plot((1/flen)*[0:flen-1],20*log10(abs(fftshift(fft(a.*hanning(flen)))))-vert);
    bottom=-150;
    ax=axis;
    ax(1)=0;
    ax(2)=1;
    ax(3)=bottom;
    ax(4)=10;
    grid on;
    scale=bottom:10:0;
    set(gca,'ytick',scale);
    axis(ax);

end;

fclose(fida);
fclose(fidb);