function [cost1,cost2]= CalCost2(Cover_Image,Secret_Image,pop)


blk=8;
cof_Cover=16;  %No of Coefficient
cof_Secret=8; % No of selected coefficients in secret image
[row col]=size(Cover_Image);
order=[1 9 2 3 10 17 25 18 11 4 5 12 19 26 33 41 34 27 20 13 6 7 14 21 28 35 42 49 57 50 43 36 29 22 15 8 16 23 30 37 44 51 58 59 52 45 38 31 24 32 39 46 53 60 61 54 47 40 48 55 62 63 56 64];
fun1=@dct2;
fun2=@idct2;

% Cover Image Processing 
J_Cover = blkproc(Cover_Image,[blk blk],fun1); % For Cover Image
x_Cover=im2col(J_Cover,[blk blk],'distinct');
x_Cover=x_Cover';  
x_Cover_zig=x_Cover(:,order);
x_Cover_Reduced=x_Cover_zig(:,1:cof_Cover); %Taking cof coefficients
no_block=size(x_Cover_Reduced,1);

% Secret Image Processing
J_Secret = blkproc(Secret_Image,[blk blk],fun1); % for secret image
x_Secret=im2col(J_Secret,[blk blk],'distinct');
x_Secret=x_Secret';  
x_Secret_zig=x_Secret(:,order);
x_Secret_Reduced=x_Secret_zig(:,1:cof_Secret);

% Embedding of secret message----------
%alpha=0.018177;
alpha= pop(1,1);
%reduced_x(:,cof+1)= reduced_x(:,cof)+ 0.001*reduced_xsecret(:,1);

for i=1:no_block
    x_Cover_Reduced(i,cof_Cover+1)= x_Cover_Reduced(i,cof_Cover)+ alpha*x_Secret_Reduced(i,1);
 %  x_Cover_Reduced(i,cof_Cover+1)= alpha*x_Secret_Reduced(i,1);
    max_c=max(x_Cover_zig(i,cof_Cover+2:cof_Cover+cof_Secret+2));
%     temp_min = x_Cover_zig(i,cof_Cover+2:cof_Cover+cof_Secret+2);
%     ind = temp_min > 0;
%     min_c = min(temp_min(ind));
    min_c=min(x_Cover_zig(i,cof_Cover+2:cof_Cover+cof_Secret+2));
    x_Cover_Reduced(i,cof_Cover+2)=max_c;
    x_Cover_Reduced(i,cof_Cover+cof_Secret+2)=min_c;
    t1=cof_Cover+3;
    t2=cof_Cover+cof_Secret+1;
   x_Cover_Reduced(i,t1:t2)=(x_Secret_Reduced(i,2:cof_Secret)-min_c)/(max_c-min_c);
   %x_Cover_Reduced(i,t1:t2)=0.01*x_Secret_Reduced(i,2:cof_Secret);
    %x_Cover_Reduced(i,cof_Cover+cof_Secret+2)=min_c;        
end

%% Stego Image Generation 
J_Stego=zeros(no_block,blk*blk);
J_Stego(:,order(1:cof_Cover+cof_Secret+2))=x_Cover_Reduced(:,:);

x_Stego=col2im(J_Stego',[blk blk],[row col],'distinct');
Stego_Image_DCTcompressed=uint8(blkproc(x_Stego,[blk blk],fun2));
imwrite(Stego_Image_DCTcompressed,'Stego_Image_DCTcompressed_Lena.jpeg');
PSNR_Steg_Image=psnr(Cover_Image,Stego_Image_DCTcompressed,255);
%figure,imshow(Stego_Image_DCTcompressed),title('Stego Image');
%figure,imhist(Cover_Image);
%figure,imhist(Stego_Image_DCTcompressed);
NAE_Cover = NAE(Cover_Image,Stego_Image_DCTcompressed);
%%%%% Recovery of secret Image 
J_Stego = blkproc(Stego_Image_DCTcompressed,[blk blk],fun1);
x_Stego=im2col(J_Stego,[blk blk],'distinct');  % taking blk x blk matrix as column.....
x_Stego=x_Stego';
x_Stego_zig=x_Stego(:,order);
Secret_recover_DC=zeros(no_block,blk*blk);

for i=1:no_block
    Secret_recover_DC(i,1) = (1/alpha)*(x_Stego_zig(i,cof_Cover+1) - x_Stego_zig(i,cof_Cover));
    %Secret_recover_DC(i,1) = (1/alpha)*x_Stego_zig(i,cof_Cover+1);
    max_c=max(x_Stego_zig(i,cof_Cover+2:cof_Cover+cof_Secret+2));
%     temp_min = x_Stego_zig(i,cof_Cover+2:cof_Cover+cof_Secret+2);
%     ind = temp_min > 0;
%     min_c = min(temp_min(ind));
     min_c=min(x_Stego_zig(i,cof_Cover+2:cof_Cover+cof_Secret+2));
    t1=cof_Cover+3;
    t2=cof_Cover+cof_Secret+1;
    Secret_recover_DC(i,order(2:cof_Secret)) = x_Stego_zig(i,t1:t2) * (max_c-min_c) + min_c;
    %Secret_recover_DC(i,order(2:cof_Secret)) = 100 * x_Stego_zig(i,t1:t2);
end

x_Secret_Recover=col2im(Secret_recover_DC',[blk blk],[row col],'distinct');
Secret_Recover_Image_DCTcompressed=uint8(blkproc(x_Secret_Recover,[blk blk],fun2));
imwrite(Secret_Recover_Image_DCTcompressed,'Secret_Recover_Image_DCTcompressed_Pepper.jpeg');
PSNR_Secret_Image=psnr(Secret_Image,Secret_Recover_Image_DCTcompressed,255);
%figure,imshow(Secret_Recover_Image_DCTcompressed),title('Secret Recovered Image');
%figure,imhist(Secret_Image);
%figure,imhist(Secret_Recover_Image_DCTcompressed);
NAE_Secret = NAE(Secret_Image,Secret_Recover_Image_DCTcompressed);
TAF_Secret = TAF(Secret_Image,Secret_Recover_Image_DCTcompressed);

cost1=-PSNR_Steg_Image;
cost2=NAE_Secret;