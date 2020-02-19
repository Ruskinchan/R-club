function [IOut,output] = denoiseImageKSVD(Image,sigma,K,varargin)%变长度输入宗量（实际输入的参量）  
%==========================================================================  
%   P E R F O R M   D E N O I S I N G   U S I N G   A  D I C T  I O N A R Y  
%                  T R A I N E D   O N   N O I S Y   I M A G E
%                  使用字典对噪声图像进行降噪
%==========================================================================  
% function IOut = denoiseImageKSVD(Image,sigma,K,varargin)  
% denoise an image by sparsely representing each block with the  
% already overcomplete trained Dictionary, and averaging the represented parts.  
% Detailed description can be found in "Image Denoising Via Sparse and Redundant  
% representations over Learned Dictionaries", (appeared in the   
% IEEE Trans. on Image Processing, Vol. 15, no. 12, December 2006).  
% This function may take some time to process. Possible factor that effect  
% the processing time are:  
%  1. number of KSVD iterations - the default number of iterations is 10.  
%  However, fewer iterations may, in most cases, result an acceleration in  
%  the process, without effecting  the result too much. Therefore, when  
%  required, this parameter may be re-set.  
%  2. maxBlocksToConsider - The maximal number of blocks to train on. If this   
%  number is larger the number of blocks in the image, random blocks  
%  from the image will be selected for training.   
% ===================================================================  
% INPUT ARGUMENTS : Image - the noisy image (gray-level scale)  
%                   sigma - the s.d. of the noise (assume to be white Gaussian).  
%                   K - the number of atoms in the trained dictionary. 训练字典的原子数目 
%    Optional arguments:                
%                  'blockSize' - the size of the blocks the algorithm  
%                       works. All blocks are squares, therefore the given  
%                       parameter should be one number (width or height).  
%                       Default value: 8.  块大小
%                       'errorFactor' - a factor that multiplies sigma in order  
%                       to set the allowed representation error. In the  
%                       experiments presented in the paper, it was set to 1.15  
%                       (which is also the default  value here).误差因子  
%                  'maxBlocksToConsider' - maximal number of blocks that  
%                       can be processed. This number is dependent on the memory  
%                       capabilities of the machine, and performances?  
%                       considerations. If the number of available blocks in the  
%                       image is larger than 'maxBlocksToConsider', the sliding  
%                       distance between the blocks increases. The default value  
%                       is: 250000.  
%                  'slidingFactor' - the sliding distance between processed  
%                       blocks. Default value is 1. However, if the image is  
%                       large, this number increases automatically (because of  
%                       memory requirements). Larger values result faster  
%                       performances (because of fewer processed blocks). 滑动因子 
%                  'numKSVDIters' - the number of KSVD iterations processed  
%                       blocks from the noisy image. If the number of  
%                       blocks in the image is larger than this number,  
%                       random blocks from all available blocks will be  
%                       selected. The default value for this parameter is:  
%                       10 if sigma > 5, and 5 otherwise.  
%                  'maxNumBlocksToTrainOn' - the maximal number of blocks  
%                       to train on. The default value for this parameter is  
%                       65000. However, it might not be enough for very large  
%                       images  
%                  'displayFlag' - if this flag is switched on,  
%                       announcement after finishing each iteration will appear,  
%                       as also a measure concerning the progress of the  
%                       algorithm (the average number of required coefficients  
%                       for representation). The default value is 1 (on). 显示标志，每次迭代结束之后显示 
%                  'waitBarOn' - can be set to either 1 or 0. If  
%                       waitBarOn==1 a waitbar, presenting the progress of the  
%                       algorithm will be displayed.  等待栏（算法进程将显示）
% OUTPUT ARGUMENTS : Iout - a 2-dimensional array in the same size of the  
%                       input image, that contains the cleaned image.与输入图像尺寸相同的矩阵，包含清理后的图像。  
%                    output.D - the trained dictionary.训练字典  
% =========================================================================  
   
% first, train a dictionary on the noisy image首先在噪声图像上进行字典训练  
   
reduceDC = 1;  
[NN1,NN2] = size(Image);%图像大小，  
waitBarOn = 1;  %等待栏
if (sigma > 5)%%%sigma=50   numIterOfKsvd = 10;  sigma是高斯噪声的标准差
    numIterOfKsvd = 10;  %默认KSVD迭代次数为10次
else  
    numIterOfKsvd = 5;  
end  
C = 1.15;  %误差因子
maxBlocksToConsider = 260000;  %考虑最大块数目？
slidingDis = 1;  %滑动因子，可随图像的增大而增大。滑动分解，重叠交叉
bb = 8;%分解块的大小  
maxNumBlocksToTrainOn = 65000;%需要训练的原子数目  
displayFlag = 1;  
hh=length(varargin)%输入参数的长度
%TF=strcmp(a,b)是字符串比较函数。相同返回值为1，否则为0.
% for argI = 1:2:length(varargin)         %实际输入的参量
%     if (strcmp(varargin{argI}, 'slidingFactor'))  
%         slidingDis = varargin{argI+1};  
%     end  
%     if (strcmp(varargin{argI}, 'errorFactor'))  
%         C = varargin{argI+1};  
%     end  
%     if (strcmp(varargin{argI}, 'maxBlocksToConsider'))  
%         maxBlocksToConsider = varargin{argI+1};  
%     end  
%     if (strcmp(varargin{argI}, 'numKSVDIters'))  
%         numIterOfKsvd = varargin{argI+1};  
%     end  
%     if (strcmp(varargin{argI}, 'blockSize'))  
%         bb = varargin{argI+1};  
%     end  
%     if (strcmp(varargin{argI}, 'maxNumBlocksToTrainOn'))  
%         maxNumBlocksToTrainOn = varargin{argI+1};  
%     end  
%     if (strcmp(varargin{argI}, 'displayFlag'))  
%         displayFlag = varargin{argI+1};  
%     end  
%     if (strcmp(varargin{argI}, 'waitBarOn'))  
%         waitBarOn = varargin{argI+1};  
%     end  
% end  
   
if (sigma <= 5)  
    numIterOfKsvd = 5;  
end  
   
% first, train a dictionary on blocks from the noisy image  
   
if(prod([NN1,NN2]-bb+1)> maxNumBlocksToTrainOn) %prob函数返回元素乘积，im2col函数中矩阵分解，子矩阵大小为bb*bb。prob返回分解后矩阵的列数。 
    randPermutation =  randperm(prod([NN1,NN2]-bb+1)); %随机置换序列。 
    selectedBlocks = randPermutation(1:maxNumBlocksToTrainOn);  
   
    blkMatrix = zeros(bb^2,maxNumBlocksToTrainOn); %全0矩阵64*65000 
    for i = 1:maxNumBlocksToTrainOn  %65000次循环
        [row,col] = ind2sub(size(Image)-bb+1,selectedBlocks(i)); %线性索引转换成相应的下标 
        currBlock = Image(row:row+bb-1,col:col+bb-1);  
        blkMatrix(:,i) = currBlock(:);  
    end  
else  
    blkMatrix = im2col(Image,[bb,bb],'sliding');%%%%%%%8*8=64   所以blkMatrix矩阵大小为：64*[（NN1-bb+1）*(NN2-bb+1)]  
end  
   
param.K = K;%%%K=256  4*8*8=256  
param.numIteration = numIterOfKsvd ;%sigma=50   所以numIterOfKsvd = 10;  
   
param.errorFlag = 1; % decompose signals until a certain error is reached. do not use fix number of coefficients.  
param.errorGoal = sigma*C;  
param.preserveDCAtom = 0;  
   
Pn=ceil(sqrt(K));%%Pn=16  %向上取整（恰好256整数开方）
DCT=zeros(bb,Pn);%%bb=8  8*16的全0矩阵
for k=0:1:Pn-1,  
    V=cos([0:1:bb-1]'*k*pi/Pn);  
    if k>0, V=V-mean(V); end;  
    DCT(:,k+1)=V/norm(V);  %V是向量，norm函数返回V的二范数。元素平方和开方。
end;  
DCT=kron(DCT,DCT);%%%%%跟DCT中的代码一样的   kronecker乘积  64*256的矩阵  
   
param.initialDictionary = DCT(:,1:param.K );%%%% 取了256列。也就是全部都取了  
param.InitializationMethod =  'GivenMatrix';  
   
if (reduceDC)%%reduceDC=1  
    vecOfMeans = mean(blkMatrix);  %求块矩阵每一列均值
    blkMatrix = blkMatrix-ones(size(blkMatrix,1),1)*vecOfMeans;%%%减去平均数  blkMatrix矩阵大小为：64*[（NN1-bb+1）*(NN2-bb+1)]  
end  
   
if (waitBarOn)%waitBarOn=1  
    counterForWaitBar = param.numIteration+1;%param.numIteration = numIterOfKsvd ;  =10  
    h = waitbar(0,'Denoising In Process ...');  
    param.waitBarHandle = h;  
    param.counterForWaitBar = counterForWaitBar;  
end  
   
   
param.displayProgress = displayFlag;%displayFlag = 1;  
[Dictionary,output] = KSVD(blkMatrix,param);%%%%%%%最核心的函数%%%%%%%%%%%
output.D = Dictionary;  
   
if (displayFlag)%displayFlag = 1;  
    disp('finished Trainning dictionary');  
end  
   
   
   
% denoise the image using the resulted dictionary  
errT = sigma*C;  
IMout=zeros(NN1,NN2);  
Weight=zeros(NN1,NN2);  
%blocks = im2col(Image,[NN1,NN2],[bb,bb],'sliding');  
while (prod(floor((size(Image)-bb)/slidingDis)+1)>maxBlocksToConsider)  
    slidingDis = slidingDis+1;  
end  
[blocks,idx] = my_im2col(Image,[bb,bb],slidingDis);  
   
if (waitBarOn)  
    newCounterForWaitBar = (param.numIteration+1)*size(blocks,2);  
end  
   
   
% go with jumps of 30000  
for jj = 1:30000:size(blocks,2)  
    if (waitBarOn)  
        waitbar(((param.numIteration*size(blocks,2))+jj)/newCounterForWaitBar);  
    end  
    jumpSize = min(jj+30000-1,size(blocks,2));  
    if (reduceDC)  
        vecOfMeans = mean(blocks(:,jj:jumpSize));  
        blocks(:,jj:jumpSize) = blocks(:,jj:jumpSize) - repmat(vecOfMeans,size(blocks,1),1);  
    end  
      
    %Coefs = mexOMPerrIterative(blocks(:,jj:jumpSize),Dictionary,errT);  
    Coefs = OMPerr(Dictionary,blocks(:,jj:jumpSize),errT);  
    if (reduceDC)  
        blocks(:,jj:jumpSize)= Dictionary*Coefs + ones(size(blocks,1),1) * vecOfMeans;  
    else  
        blocks(:,jj:jumpSize)= Dictionary*Coefs ;  
    end  
end  
   
count = 1;  
Weight = zeros(NN1,NN2);  
IMout = zeros(NN1,NN2);  
[rows,cols] = ind2sub(size(Image)-bb+1,idx);  
for i  = 1:length(cols)  
    col = cols(i); row = rows(i);          
    block =reshape(blocks(:,count),[bb,bb]);  
    IMout(row:row+bb-1,col:col+bb-1)=IMout(row:row+bb-1,col:col+bb-1)+block;  
    Weight(row:row+bb-1,col:col+bb-1)=Weight(row:row+bb-1,col:col+bb-1)+ones(bb);  
    count = count+1;  
end;  
   
if (waitBarOn)  
    close(h);  
end  
IOut = (Image+0.034*sigma*IMout)./(1+0.034*sigma*Weight);  
