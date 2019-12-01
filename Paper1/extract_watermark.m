function y=ext_water()

water_marked = imread('Watermarked.png');

if length(size(water_marked))>2
    water_marked=rgb2gray(water_marked);
end

water_marked = imresize(water_marked,[512,512]);

[LL1,LH1,HL1,HH1]=dwt2(water_marked,'haar');
[LL2,LH2,HL2,HH2]=dwt2(LL1,'haar');
[LL3,LH3,HL3,HH3]=dwt2(LL2,'haar');

beta = zeros(1024,1);

num=1;
q=78;
for i=1:64
    for j=1:8:64
        C= LH3(i,j:j+7);
        z = floor(mean(C)/(q*1.0)+0.5);
        beta(num,1) = mod(z,2);
        num=num+1;
    end
end

for i=1:64
    for j=1:8:64
        C= HL3(i,j:j+7);
        z = floor(mean(C)/(q*1.0)+0.5);
        beta(num,1) = mod(z,2);
        num=num+1;
    end
end

watermark = zeros(32,32);

for i=0:31
    watermark(i+1,1:32) = beta(i*32+1:i*32+32,1);
end

water_mark=imread('Peppers.bmp');
if length(size(water_mark))>2
    host=rgb2gray(water_mark);
end
water_mark=imbinarize(water_mark);
water_mark=imresize(water_mark,[32, 32]);

[number,ratio]=biterr(water_mark,watermark);
disp(number);
disp(ratio);

y= 'watermark extracted';
end


