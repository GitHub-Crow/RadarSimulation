function [vectorMag, vectorMag_dB, vectorMag_dB_Norm] = processVector(vectorArg)
vectorMag = abs(vectorArg);
vectorMag_dB = 20*log10(vectorMag + eps);
vectorMag_dB_Norm = 20*log10((vectorMag + eps)/max(vectorMag));
end

