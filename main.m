%============================================================
%               demo2 - denoise an image图像去噪
% this is a run_file the demonstrate how to denoise an image, 文件演示用字典如何去噪
% using dictionaries. The methods implemented here are the same
% one as described in "Image Denoising Via Sparse and Redundant
% representations over Learned Dictionaries", (appeared in the 
% IEEE Trans. on Image Processing, Vol. 15, no. 12, December 2006).
% ============================================================
 
clear
bb=8; % block size块大小
RR=4; % redundancy factor  冗余因素
K=RR*bb^2; % number of atoms in the dictionary字典原子的数目.k=256,运算顺序先幂运算再乘法。
sigma = 50; %噪声均方误差
%pathForImages ='';图像路径=''？；
%imageName = 'barbara.png';图像名字；
%   [IMin0,pp]=imread('cameraman.tif');IMin0表示原始图像
 [IMin0,pp]=imread('lena.jpg');%imread读取图像数据。
IMin0=im2double(IMin0);%像素值被归一化
%图像数据由矩阵表示。不同格式数据不同。
%彩色图片（rgb）格式包含三种颜色信息，每个颜色包含颜色深度。
%灰度图像（gray）格式只包含颜色深度。
%彩色图像的矩阵维度多于灰度图像。
if (length(size(IMin0))>2)
    IMin0 = rgb2gray(IMin0);%提取图像维度，并判断，若是彩色图片，则转变为灰度图像
end
if (max(IMin0(:))<2)
    IMin0 = IMin0*255;%像素值转换
    %判断图像深度的表示方式，常规下是整数。如果是小数表示（最大值不超过2）就转换为整数形式。
end
%IMin0表示原始图像
%IMin表示加噪图像
IMin=IMin0+sigma*randn(size(IMin0));%%%%%%此处有随机函数%%%%%50*随机矩阵
PSNRIn = 20*log10(255/sqrt(mean((IMin(:)-IMin0(:)).^2)));
%输入图像的峰值信噪比=20lg(255/加噪图像与原始图像差值平方的均值再开方)
tic%开始一个秒表计时器
[IoutAdaptive,output] = denoiseImageKSVD(IMin, sigma,K);
 %调用函数式
PSNROut = 20*log10(255/sqrt(mean((IoutAdaptive(:)-IMin0(:)).^2)));
%输出图像的峰值信噪比=20lg（255/输出自适应图像与原始图像差值平方的均值再开方）
%ooooooooooooooooooooooooooooooooooooooooooooooooo%
figure;%图形
subplot(1,3,1); imshow(IMin0,[]); title('输入图像');
subplot(1,3,2); imshow(IMin,[]); title(strcat(['噪声图像,',num2str(PSNRIn),'dB']));
subplot(1,3,3); imshow(IoutAdaptive,[]); title(strcat(['去噪后图像, ',num2str(PSNROut),'dB']));
figure;
%在一行三列的画布上，第一列原始图像，第二列加噪图像，第三列输出自适应图像。
%其中tittle(strcat表示先转换为字符串，再和字符连接。
%oooooooooooooooooooooooooooooooooooooooooooooooo%
I = displayDictionaryElementsAsImage(output.D, floor(sqrt(K)),floor(size(output.D,2)/floor(sqrt(K))),bb,bb);
title('从噪声图像中训练的块字典');
%/////////////////////////////////////////////////////////上式
toc%tic到toc的时间间隔，单位：s。