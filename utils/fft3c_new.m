function [res] = fft3c_new(x)
% fctr = 1%size(x,1)*size(x,2)*size(x,3);

% 	res = 1/sqrt(fctr)*fftshift(fftshift(fftshift(fft(fft(fft(ifftshift(ifftshift(ifftshift(x,1),2),3),[],1),[],2),[],3),1),2),3);
res = fftc_new(fftc_new(fftc_new(x,1),2),3);
end

