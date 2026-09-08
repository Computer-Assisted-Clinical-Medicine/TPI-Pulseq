
    function [imc]=my_coilcombine2(IMs,cdim)
    % SoS combination
    imc=sqrt(sum(abs(IMs.^2),cdim));
    imph = sum(angle(IMs).*abs(IMs),cdim)./imc;
    imc = imc.*exp(1i*imph);
    end