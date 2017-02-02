signal_name = '/home/darryn/Git/nextrad-processor/reference_signals/matlab_format/Ref_PL_Xband_MPALNA_7_5us';
in_file_name = strcat(signal_name, '.mat');
out_file_name = strcat(signal_name, '.dat');

load(in_file_name);

f_out_id = fopen(out_file_name,'w');

for i = 1 : length(Ref_sig)
    fwrite(f_out_id, real(Ref_sig(i)), 'int16', 'ieee-le');
    fwrite(f_out_id, imag(Ref_sig(i)), 'int16', 'ieee-le');
end;

fclose(f_out_id);

figure(1)
hold on;
plot(real(Ref_sig));
plot(imag(Ref_sig));
legend('real', 'imag');
