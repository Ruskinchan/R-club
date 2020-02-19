    function [blocks,idx] = my_im2col(I,blkSize,slidingDis)
       
    if (slidingDis==1)  
        blocks = im2col(I,blkSize,'sliding');%行为blksize元素的总个数，列为(m-bb+1) x (n-bb+1)=62001     
        % http://fuda641.blog.163.com/blog/static/20751421620135483846711/  
        idx =1:size(blocks,2); %size函数取新矩阵的列数，idx取值范围从1到62001 
        return  
    end  
       
    idxMat = zeros(size(I)-blkSize+1);  %生成0矩阵
    idxMat([1:slidingDis:end-1,end],[1:slidingDis:end-1,end]) = 1;
    % take blocks in distances of 'slidingDis',从滑行距离中取矩阵，但始终取每行每列的第一个和最后一个
    % but always take the first and last one (in each row and column).  
    idx = find(idxMat);  %idx即矩阵中非零元素的位置
    [rows,cols] = ind2sub(size(idxMat),idx);  
    blocks = zeros(prod(blkSize),length(idx));%prob子矩阵列元素的乘积，返回组成行向量。  
    for i = 1:length(idx)  
        currBlock = I(rows(i):rows(i)+blkSize(1)-1,cols(i):cols(i)+blkSize(2)-1);  
        blocks(:,i) = currBlock(:);  
    end  