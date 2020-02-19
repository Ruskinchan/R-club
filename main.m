%============================================================
%               demo2 - denoise an imageͼ��ȥ��
% this is a run_file the demonstrate how to denoise an image, �ļ���ʾ���ֵ����ȥ��
% using dictionaries. The methods implemented here are the same
% one as described in "Image Denoising Via Sparse and Redundant
% representations over Learned Dictionaries", (appeared in the 
% IEEE Trans. on Image Processing, Vol. 15, no. 12, December 2006).
% ============================================================
 
clear
bb=8; % block size���С
RR=4; % redundancy factor  ��������
K=RR*bb^2; % number of atoms in the dictionary�ֵ�ԭ�ӵ���Ŀ.k=256,����˳�����������ٳ˷���
sigma = 50; %�����������
%pathForImages ='';ͼ��·��=''����
%imageName = 'barbara.png';ͼ�����֣�
%   [IMin0,pp]=imread('cameraman.tif');IMin0��ʾԭʼͼ��
 [IMin0,pp]=imread('lena.jpg');%imread��ȡͼ�����ݡ�
IMin0=im2double(IMin0);%����ֵ����һ��
%ͼ�������ɾ����ʾ����ͬ��ʽ���ݲ�ͬ��
%��ɫͼƬ��rgb����ʽ����������ɫ��Ϣ��ÿ����ɫ������ɫ��ȡ�
%�Ҷ�ͼ��gray����ʽֻ������ɫ��ȡ�
%��ɫͼ��ľ���ά�ȶ��ڻҶ�ͼ��
if (length(size(IMin0))>2)
    IMin0 = rgb2gray(IMin0);%��ȡͼ��ά�ȣ����жϣ����ǲ�ɫͼƬ����ת��Ϊ�Ҷ�ͼ��
end
if (max(IMin0(:))<2)
    IMin0 = IMin0*255;%����ֵת��
    %�ж�ͼ����ȵı�ʾ��ʽ���������������������С����ʾ�����ֵ������2����ת��Ϊ������ʽ��
end
%IMin0��ʾԭʼͼ��
%IMin��ʾ����ͼ��
IMin=IMin0+sigma*randn(size(IMin0));%%%%%%�˴����������%%%%%50*�������
PSNRIn = 20*log10(255/sqrt(mean((IMin(:)-IMin0(:)).^2)));
%����ͼ��ķ�ֵ�����=20lg(255/����ͼ����ԭʼͼ���ֵƽ���ľ�ֵ�ٿ���)
tic%��ʼһ������ʱ��
[IoutAdaptive,output] = denoiseImageKSVD(IMin, sigma,K);
 %���ú���ʽ
PSNROut = 20*log10(255/sqrt(mean((IoutAdaptive(:)-IMin0(:)).^2)));
%���ͼ��ķ�ֵ�����=20lg��255/�������Ӧͼ����ԭʼͼ���ֵƽ���ľ�ֵ�ٿ�����
%ooooooooooooooooooooooooooooooooooooooooooooooooo%
figure;%ͼ��
subplot(1,3,1); imshow(IMin0,[]); title('����ͼ��');
subplot(1,3,2); imshow(IMin,[]); title(strcat(['����ͼ��,',num2str(PSNRIn),'dB']));
subplot(1,3,3); imshow(IoutAdaptive,[]); title(strcat(['ȥ���ͼ��, ',num2str(PSNROut),'dB']));
figure;
%��һ�����еĻ����ϣ���һ��ԭʼͼ�񣬵ڶ��м���ͼ�񣬵������������Ӧͼ��
%����tittle(strcat��ʾ��ת��Ϊ�ַ������ٺ��ַ����ӡ�
%oooooooooooooooooooooooooooooooooooooooooooooooo%
I = displayDictionaryElementsAsImage(output.D, floor(sqrt(K)),floor(size(output.D,2)/floor(sqrt(K))),bb,bb);
title('������ͼ����ѵ���Ŀ��ֵ�');
%/////////////////////////////////////////////////////////��ʽ
toc%tic��toc��ʱ��������λ��s��