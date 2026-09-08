function out = hamming_3d(n1,n2,n3)
% Return a symmetrical 3D hamming window of dimensions n1 x n2 x n3 
% The centre is at the origin [n1/2 n2/2 n3/2] 
% The filtered is scaled such that max(out(:)) = 1 
% 
% Jack J. Miller, University of Oxford 2017 
if ~(isscalar(n1) && isscalar(n2) && isscalar(n3))
    error('Input arguments should be scalar'); 
end
if ((n2 == 1) && (n3 ==1))
    out = hamming(n1); 
elseif n3 == 1 
    out = bsxfun(@times,hamming(n1),hamming(n2).'); 
else
    out = bsxfun(@times,bsxfun(@times,hamming(n1),hamming(n2).'),permute(hamming(n3),[3 2 1])); 
end
out = out./(max(out(:))); 
end
function h = hamming(x) % 1D window of length x
% C.f. https://en.wikipedia.org/wiki/Window_function#Hamming_window
% This overrides Matlab's inbuilt window function. 

alpha=.54;%.54;  %0.54 hamming ! 0.5 hanning! 
beta=1-alpha; 
h = alpha-beta*cos(2*pi*[0:(x-1)]./(x-1));    % hamming
% h = alpha*(1-cos(2*pi*[0:(x-1)]./(x-1)));      %hanning
h = h.'; 
% h = kaiser(x,2.5);
end