function [res] = ifft3c_new(x)
% function res = ifft3c_new(x)
% fctr = 1%size(x,1)*size(x,2)*size(x,3);

% 	res = sqrt(fctr)*ifftshift(ifftshift(ifftshift(ifft(ifft(ifft(fftshift(fftshift(fftshift(x,1),2),3),[],1),[],2),[],3),1),2),3);
res = ifftc_new(ifftc_new(ifftc_new(x,1),2),3);
end

 