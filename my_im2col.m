    function [blocks,idx] = my_im2col(I,blkSize,slidingDis)
       
    if (slidingDis==1)  
        blocks = im2col(I,blkSize,'sliding');%��ΪblksizeԪ�ص��ܸ�������Ϊ(m-bb+1) x (n-bb+1)=62001     
        % http://fuda641.blog.163.com/blog/static/20751421620135483846711/  
        idx =1:size(blocks,2); %size����ȡ�¾����������idxȡֵ��Χ��1��62001 
        return  
    end  
       
    idxMat = zeros(size(I)-blkSize+1);  %����0����
    idxMat([1:slidingDis:end-1,end],[1:slidingDis:end-1,end]) = 1;
    % take blocks in distances of 'slidingDis',�ӻ��о�����ȡ���󣬵�ʼ��ȡÿ��ÿ�еĵ�һ�������һ��
    % but always take the first and last one (in each row and column).  
    idx = find(idxMat);  %idx�������з���Ԫ�ص�λ��
    [rows,cols] = ind2sub(size(idxMat),idx);  
    blocks = zeros(prod(blkSize),length(idx));%prob�Ӿ�����Ԫ�صĳ˻������������������  
    for i = 1:length(idx)  
        currBlock = I(rows(i):rows(i)+blkSize(1)-1,cols(i):cols(i)+blkSize(2)-1);  
        blocks(:,i) = currBlock(:);  
    end  