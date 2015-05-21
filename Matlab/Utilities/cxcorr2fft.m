function R=cxcorr2fft(im1,im2)

R=ifft2(fft2(im1/norm(im1(:))).*conj(fft2(im2/norm(im2(:)))));

if(~isreal(R))
        R=abs(R);
end