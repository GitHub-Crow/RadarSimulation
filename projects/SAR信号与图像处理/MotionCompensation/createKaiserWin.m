function win = createKaiserWin(WL, beta, len, center)
im = zeros(len, 1);
im(center) = 1;
win = conv(im, ones(WL, 1), 'same'); 
WL = WL + 2*(WL - sum(win));
factor = kaiser(WL, beta);
win = conv(im, factor, 'same');
end