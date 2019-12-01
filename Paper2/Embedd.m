function cost=Embedd(img,wt,pop)

first=pop(1,1);
second=pop(1,2);

k=1; dr=0; dc=0;
% dr is to address 1:8 row every time for new block in x
% dc is to address 1:8 column every time for new block in x
% k is to change the no. of cell

%%%%%%%%%%%%%%%%% To divide image in to 4096---8X8 blocks %%%%%%%%%%%%%%%%%%
x={};
for ii=1:8:512 % To address row -- 8X8 blocks of image
    for jj=1:8:512 % To address columns -- 8X8 blocks of image
        for i=ii:(ii+7) % To address rows of blocks
            dr=dr+1;
            for j=jj:(jj+7) % To address columns of block
                dc=dc+1;
                z(dr,dc)=img(i,j);
            end
            dc=0;
        end
        x{k}=dct2(z); k=k+1;
        z=[]; dr=0;
    end
end

% To find the ideal (i,j) which have min argument
ijs = zeros(4032,2);
for k=1:4032
    kx=(x{k});
    kx1=(x{k+64});
    minarg=1000;
    for i=1:8
        for j=1:8
            val = abs(kx(i,j)-kx1(i,j));
            if(val<minarg)
                minarg=val;
                ijs(k,1)=i;
                ijs(k,2)=j;
            end
        end
    end
end


welem=1; wmrk={};
for i=1:63
    for j=1:64
        wmrk{welem}=wt(i,j);
        welem=welem+1;
    end
end

welem=welem-1;

for k=1:4032
    kx=x{k};
    kx1=x{k+64};
    i=ijs(k,1); j=ijs(k,2);
    if ((wmrk{k}==1 && kx(i,j)<kx1(i,j)) || (wmrk{k}==0 && kx(i,j)>kx1(i,j)))
        temp=kx(i,j);
        kx(i,j)=kx1(i,j);
        kx1(i,j)=temp;
    end
    
    if (wmrk{k}==1 && abs(kx(i,j)-kx1(i,j))<first)
        kx(i,j)=kx(i,j)+second*first;
        kx1(i,j)=kx1(i,j)-(1-second)*first;
    elseif (wmrk{k}==0 && abs(kx(i,j)-kx1(i,j))>first)
        kx(i,j)=kx(i,j)-second*first;
        kx1(i,j)=kx1(i,j)+(1-second)*first;
    end
    
    x{k}=kx;
    x{k+64}=kx1;   
end

for k=1:4096
    val= idct2(x{k});
    x{k}= val;
end

water_marked=zeros(512,512);

for k=1:4096
    row=(ceil(k/64)-1)*8+1;
    col=(mod(k-1,64))*8+1;
    water_marked(row:row+7,col:col+7)=x{k};
end

subplot(1,3,1)
imshow(uint8(water_marked));
title('Watermarked image');

imwrite(uint8(water_marked),'Watermarked.png');
y='Watermarked succesfully';
cost = -psnr(img,uint8(water_marked)); 