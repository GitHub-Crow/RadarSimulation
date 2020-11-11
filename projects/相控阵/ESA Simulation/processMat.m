function [mag, magdB, magdB_Norm] = processMat(mat)
mag = abs(mat);
magMax = max(max(mag));
magdB = 10*log10(mag);
magdB_Norm = magdB - 10*log10(magMax);
end