function [Dictionary,output] = KSVD(...  
    Data,...   
    param)  
  
if (~isfield(param,'displayProgress'))  
    param.displayProgress = 0;  
end  
totalerr(1) = 99999;  
if (isfield(param,'errorFlag')==0)
    param.errorFlag = 0;  
end  
   
if (isfield(param,'TrueDictionary')) 
    displayErrorWithTrueDictionary = 1;  
    ErrorBetweenDictionaries = zeros(param.numIteration+1,1);  
    ratio = zeros(param.numIteration+1,1);  
else  
    displayErrorWithTrueDictionary = 0;  
    ratio = 0;%
end  
if (param.preserveDCAtom>0)  %param.preserveDCAtom = 0;  
    FixedDictionaryElement(1:size(Data,1),1) = 1/sqrt(size(Data,1));  
else  
    FixedDictionaryElement = [];  
end  
  
if (size(Data,2) < param.K)%K=256    size(Data,2)=249*249    
    disp('Size of data is smaller than the dictionary size. Trivial solution...');  
    Dictionary = Data(:,1:size(Data,2));  
    return;  
elseif (strcmp(param.InitializationMethod,'DataElements'))  
    Dictionary(:,1:param.K-param.preserveDCAtom) = Data(:,1:param.K-param.preserveDCAtom);  
elseif (strcmp(param.InitializationMethod,'GivenMatrix')) 
    Dictionary(:,1:param.K-param.preserveDCAtom) = param.initialDictionary(:,1:param.K-param.preserveDCAtom);CT(:,1:param.K );
end  
if (param.preserveDCAtom) 
    tmpMat = FixedDictionaryElement \ Dictionary;  
    Dictionary = Dictionary - FixedDictionaryElement*tmpMat;  
end  

Dictionary = Dictionary*diag(1./sqrt(sum(Dictionary.*Dictionary)));%64*256      
Dictionary = Dictionary.*repmat(sign(Dictionary(1,:)),size(Dictionary,1),1);
totalErr = zeros(1,param.numIteration);
  
for iterNum = 1:param.numIteration 
    if (param.errorFlag==0)   
        CoefMatrix = OMP([FixedDictionaryElement,Dictionary],Data, param.L); 
    else    
        CoefMatrix = OMPerr([FixedDictionaryElement,Dictionary],Data, param.errorGoal); 
        param.L = 1;  
    end  
      
    replacedVectorCounter = 0;  
    rPerm = randperm(size(Dictionary,2)) 
    for j = rPerm   
        [betterDictionaryElement,CoefMatrix,addedNewVector] = I_findBetterDictionaryElement(Data,...     
            [FixedDictionaryElement,Dictionary],j+size(FixedDictionaryElement,2),...  
            CoefMatrix,param.L);  
        Dictionary(:,j) = betterDictionaryElement;
        if (param.preserveDCAtom) 
            tmpCoef = FixedDictionaryElement\betterDictionaryElement;  
            Dictionary(:,j) = betterDictionaryElement - FixedDictionaryElement*tmpCoef;  
            Dictionary(:,j) = Dictionary(:,j)./sqrt(Dictionary(:,j)'*Dictionary(:,j));  
        end  
        replacedVectorCounter = replacedVectorCounter+addedNewVector; 
    end  
   
      
    if (iterNum>1 & param.displayProgress)
        if (param.errorFlag==0)  
            output.totalerr(iterNum-1) = sqrt(sum(sum((Data-[FixedDictionaryElement,Dictionary]*CoefMatrix).^2))/prod(size(Data)));  
            disp(['Iteration   ',num2str(iterNum),'   Total error is: ',num2str(output.totalerr(iterNum-1))]);  
        else %执行此句  
            output.numCoef(iterNum-1) = length(find(CoefMatrix))/size(Data,2); 
            disp(['Iteration   ',num2str(iterNum),'   Average number of coefficients: ',num2str(output.numCoef(iterNum-1))]);  
        end  
    end  
    if (displayErrorWithTrueDictionary ) 
        [ratio(iterNum+1),ErrorBetweenDictionaries(iterNum+1)] = I_findDistanseBetweenDictionaries(param.TrueDictionary,Dictionary)
        disp(strcat(['Iteration  ', num2str(iterNum),' ratio of restored elements: ',num2str(ratio(iterNum+1))]));  
        output.ratio = ratio;  
    end  
      
   Dictionary = I_clearDictionary(Dictionary,CoefMatrix(size(FixedDictionaryElement,2)+1:end,:),Data);
         
%     h = waitbar(0,'Denoising In Process ...');  
%     param.waitBarHandle = h;  
    if (isfield(param,'waitBarHandle'))  
        waitbar(iterNum/param.counterForWaitBar);  
    end  
end  
   
output.CoefMatrix = CoefMatrix;  
Dictionary = [FixedDictionaryElement,Dictionary];
function [betterDictionaryElement,CoefMatrix,NewVectorAdded] = I_findBetterDictionaryElement(Data,Dictionary,j,CoefMatrix,numCoefUsed)  
if (length(who('numCoefUsed'))==0)  
    numCoefUsed = 1;  
  
end  
relevantDataIndices = find(CoefMatrix(j,:));
if (length(relevantDataIndices)<1) 
    ErrorMat = Data-Dictionary*CoefMatrix;  
    ErrorNormVec = sum(ErrorMat.^2);  
    [d,i] = max(ErrorNormVec);  
    betterDictionaryElement = Data(:,i);%ErrorMat(:,i); %  
    betterDictionaryElement = betterDictionaryElement./sqrt(betterDictionaryElement'*betterDictionaryElement);  
    betterDictionaryElement = betterDictionaryElement.*sign(betterDictionaryElement(1));  
    CoefMatrix(j,:) = 0;  
    NewVectorAdded = 1
    return;  
end  
   
NewVectorAdded = 0;  
tmpCoefMatrix = CoefMatrix(:,relevantDataIndices);  
tmpCoefMatrix(j,:) = 0; 
errors =(Data(:,relevantDataIndices) - Dictionary*tmpCoefMatrix);

[betterDictionaryElement,singularValue,betaVector] = svds(errors,1);
%a=[1 2 3 4;5 6 7 8;9 10 11 12;2 4 6 7.99999]; [u,s,v]=svds(a)   u*s*v'    [u,s,v]=svds(a,1):取出的第一主成分   
%对于svds函数：a为M*N的矩阵，那么u:M*M   S:M*N(简写成M*M)   V=N*M    V'=M*N  
%对于svd函数：a为M*N的矩阵， 那么u:M*M   S:M*N             V=N*N    V'=N*N  
%将字典原子D的解定义为U中的第一列，将系数向量CoefMatrix的解定义为V的第一列与S(1，1)的乘积  
CoefMatrix(j,relevantDataIndices) = singularValue*betaVector';

function [ratio,totalDistances] = I_findDistanseBetweenDictionaries(original,new)  
catchCounter = 0;  
totalDistances = 0;  
for i = 1:size(new,2)  
    new(:,i) = sign(new(1,i))*new(:,i);  
end  
for i = 1:size(original,2)  
    d = sign(original(1,i))*original(:,i);  
    distances =sum ( (new-repmat(d,1,size(new,2))).^2);  
    [minValue,index] = min(distances);  
    errorOfElement = 1-abs(new(:,index)'*d);  
    totalDistances = totalDistances+errorOfElement;  
    catchCounter = catchCounter+(errorOfElement<0.01);  
end  
ratio = 100*catchCounter/size(original,2);  

function Dictionary = I_clearDictionary(Dictionary,CoefMatrix,Data)  
T2 = 0.99;  
T1 = 3;  
K=size(Dictionary,2);   
Er=sum((Data-Dictionary*CoefMatrix).^2,1); *relevantDataIndices  
G=Dictionary'*Dictionary; 
G = G-diag(diag(G));
for jj=1:1:K,  
    if max(G(jj,:))>T2 | length(find(abs(CoefMatrix(jj,:))>1e-7))<=T1 ,  
        [val,pos]=max(Er);  
        clearDictionary=1 
        Er(pos(1))=0;%
        Dictionary(:,jj)=Data(:,pos(1))/norm(Data(:,pos(1)));
        G=Dictionary'*Dictionary;  
        G = G-diag(diag(G));  
    end;  
end;
