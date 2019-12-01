function y= water()

host=imread('Lena.bmp');
if length(size(host))>2
    host=rgb2gray(host);
end
host=imresize(host,[512, 512]);
%subplot(1,3,1)
%imshow(host);
%title('Original Image');
[LL1,LH1,HL1,HH1]=dwt2(host,'haar');
[LL2,LH2,HL2,HH2]=dwt2(LL1,'haar');
[LL3,LH3,HL3,HH3]=dwt2(LL2,'haar');

water_mark=imread('Peppers.bmp');
if length(size(water_mark))>2
    host=rgb2gray(water_mark);
end
water_mark=imbinarize(water_mark);
water_mark=imresize(water_mark,[32, 32]);

beta = zeros(1024,1);

for i=0:31
    beta(i*32+1:i*32+32,1)= water_mark(i+1,1:32);
end

q= 78;
num=1;
for i=1:64
    for j=1:8:64
        C= LH3(i,j:j+7);
        z = floor(mean(C)/(q*1.0)+0.5);
        z1 = floor(mean(C)/(q*1.0));
        z = mod(z,2);
        y=0;
        if z== beta(num,1)
            y= z*q;
        elseif z==z1
            y= (z+1)*q;
        else
            y= (z-1)*q;
        end
        LH3(i,j:j+7)=C-(mean(C)-y);
        num=num+1;
    end
end

for i=1:64
    for j=1:8:64
        C= HL3(i,j:j+7);
        z = floor(mean(C)/(q*1.0)+0.5);
        z1 = floor(mean(C)/(q*1.0));
        z = mod(z,2);
        y=0;
        if z== beta(num,1)
            y= z*q;
        elseif z==z1
            y= (z+1)*q;
        else
            y= (z-1)*q;
        end
        HL3(i,j:j+7)=C-(mean(C)-y);
        num=num+1;
    end
end


LL2= idwt2(LL3,LH3,HL3,HH3,'haar');
LL1= idwt2(LL2,LH2,HL2,HH2,'haar');
watermarked = idwt2(LL1,LH1,HL1,HH1,'haar');
subplot(1,3,1)
imshow(uint8(watermarked));
title('Watermarked image');

imwrite(uint8(watermarked),'Watermarked.png');
y='Watermarked succesfully';
PSNR = psnr(host, uint8(watermarked))
end