
clc;
clear;
close all;

CostFunction=@(img,wt,x) CalCost(img,wt,x);        % Cost Function

nVar=2;             % Number of Decision Variables  originally nVar=5
VarSize=[1 nVar];   % Decision Variables Matrix Size

%VarMin= -10;         % Decision Variables Lower Bound  % Originally VarMin = 0;
%VarMax= 10;         % Decision Variables Upper Bound  % VarMax = 1;

MaxIt=10;              % Maximum Number of Iterations
nPop=30;               % Population Size (Colony Size)
nOnlooker=nPop;         % Number of Onlooker Bees
L=6; % Abandonment Limit Parameter (Trial Limit)
a=1;                    % Acceleration Coefficient Upper Bound
%fun = @(block_struct) dct2(block_struct);
im=imread('Lena.bmp'); % im=imread(im);
if length(size(im))>2
    im=rgb2gray(im);
end
im = imresize(im,[512 512]); % Resize image
%dct_img=blockproc(im,[8,8],fun);% DCT of image using 8X8 block

wt = imread('tiffany.bmp');
if length(size(wt))>2
    wt=rgb2gray(wt);
end
wt = imbinarize(wt);
imwrite(uint8(wt),'WatermarkOrg.png');
watermark = imresize(wt,[63 64]);% Resize and Change in binary 




%% Initialization
% Empty Bee Structureempty_bee.Position=[];
empty_bee.Cost=[];
% Initialize Population Array
pop=repmat(empty_bee,nPop,1);

% Initialize Best Solution Ever Found
BestSol.Cost=[inf 0];
flag=0;
% Create Initial Population
for i=1:nPop
    %pop(i).Position=unifrnd(VarMin,VarMax,VarSize);%pop(i).Position=randi([VarMin,VarMax],VarSize); for 0/1
    pop(i).Position(1,1)= 0+rand*80;
    pop(i).Position(1,2)= rand;
    [pop(i).Cost(1),pop(i).Cost(2)]=CostFunction(im,watermark,pop(i).Position);
    %disp(pop(i).Cost);
    if (flag==0 && pop(i).Cost(1)>-36 && pop(i).Cost(1)<=BestSol.Cost(1))
        BestSol=pop(i);
    elseif (flag==0 && pop(i).Cost(1)<=-36 && pop(i).Cost(1)>=-49)
        BestSol=pop(i);
        flag=1;
    elseif (flag==1 && pop(i).Cost(2)<BestSol.Cost(2) && pop(i).Cost(1)>=-49)
        BestSol=pop(i);
    end
end

% Abandonment Counter
C=zeros(nPop,1);

% Array to Hold Best Cost Values
BestCost=zeros(MaxIt,2);


%% ABC Main Loop

for it=1:MaxIt
    
    % Recruited Bees
    for i=1:nPop
        
        % Choose k randomly, not equal to i
        K=[1:i-1 i+1:nPop];
        b= numel(K);
        k=K(randi([1,b]));
        
        % Define Acceleration Coeff.
        phi=a*unifrnd(-1,+1,VarSize);  % 
        
        % New Bee Position
        newbee.Position=pop(i).Position+phi.*(pop(i).Position-pop(k).Position);
        
        % Evaluation
        [newbee.Cost(1),newbee.Cost(2)]=CostFunction(im,watermark,newbee.Position);
        
        % Comparision
        if (newbee.Cost(1)>-36 && pop(i).Cost(1)>-36 && newbee.Cost(1)<=pop(i).Cost(1))
            pop(i)=newbee;
        elseif (newbee.Cost(1)<=-36 && pop(i).Cost(1)>-36 && newbee.Cost(1)>=-49)
            pop(i)=newbee;
        elseif (newbee.Cost(1)>=-49 && pop(i).Cost(1)<-49)
            pop(i)=newbee;
        elseif (newbee.Cost(1)<=-36 && pop(i).Cost(1)<=-36 && newbee.Cost(2)<pop(i).Cost(2) && newbee.Cost(1)>=-49 && pop(i).Cost(1)>=-49)
            pop(i)=newbee;
        else
            C(i)=C(i)+1;
        end
        
    end
    % Calculate Fitness Values and Selection Probabilities
    F=zeros(nPop,1);
    if flag==0
        sum1=0;
        for i=1:nPop
            sum1=sum1+pop(i).Cost(1);
        end   
        MeanCost = sum1/nPop;
        for i=1:nPop
            F(i) = exp(-pop(i).Cost(1)/MeanCost); % Convert Cost to Fitness
        end
    else
        sum1=0;
        for i=1:nPop
            sum1=sum1+pop(i).Cost(2);
        end   
        MeanCost = sum1/nPop;
        for i=1:nPop
            F(i) = exp(-pop(i).Cost(2)/MeanCost); % Convert Cost to Fitness
        end
    end
    P=F/sum(F);
    
    % Onlooker Bees
    for m=1:nOnlooker
        
        % Select Source Site
        i=RouletteWheelSelection(P);
        
        % Choose k randomly, not equal to i
        K=[1:i-1 i+1:nPop];
        b1=numel(K);
        %fprintf('%d, ',b1);
        k=K(randi([1,b1]));
        
        % Define Acceleration Coeff.
        phi=a*unifrnd(-1,+1,VarSize);
        
        % New Bee Position
        newbee.Position=pop(i).Position+phi.*(pop(i).Position-pop(k).Position);
        
        % Evaluation
        [newbee.Cost(1),newbee.Cost(2)]=CostFunction(im,watermark,newbee.Position);
        
        % Comparision
        if (newbee.Cost(1)>-36 && pop(i).Cost(1)>-36 && newbee.Cost(1)<=pop(i).Cost(1))
            pop(i)=newbee;
        elseif (newbee.Cost(1)<=-36 && pop(i).Cost(1)>-36 && newbee.Cost(1)>=-49)
            pop(i)=newbee;
        elseif (newbee.Cost(1)>=-49 && pop(i).Cost(1)<-49)
            pop(i)=newbee;
        elseif (newbee.Cost(1)<=-36 && pop(i).Cost(1)<=-36 && newbee.Cost(2)<pop(i).Cost(2) && newbee.Cost(1)>=-49 && pop(i).Cost(1)>=-49)
            pop(i)=newbee;
        else
            C(i)=C(i)+1;
        end
    end
    
    % Scout Bees
    for i=1:nPop
        if C(i)>=L
            pop(i).Position(1,1)= rand*80;
            pop(i).Position(1,2)= rand; %pop(i).Position=randi([VarMin,VarMax],VarSize);for 0/1
            [pop(i).Cost(1),pop(i).Cost(2)]=CostFunction(im,watermark,pop(i).Position);
            C(i)=0;
        end
    end
    
    % Update Best Solution Ever Found
    for i=1:nPop
        if (flag==0 && pop(i).Cost(1)>-36 && pop(i).Cost(1)<=BestSol.Cost(1))
            BestSol=pop(i);
        elseif (flag==0 && pop(i).Cost(1)<=-36 && pop(i).Cost(1)>=-49)
            BestSol=pop(i);
            flag=1;
        elseif (flag==1 && pop(i).Cost(2)<BestSol.Cost(2) && pop(i).Cost(1)>=-49)
            BestSol=pop(i);
        end
    end
    
    % Store Best Cost Ever Found
    BestCost(it,1)=BestSol.Cost(1);
    BestCost(it,2)=BestSol.Cost(2);
    
    % Display Iteration Information
    disp(['Iteration ' num2str(it) ': Best Cost = [' num2str(BestCost(it,1)) ',' num2str(BestCost(it,2)) '], k= ' num2str(BestSol.Position(1,1)) ',p=' num2str(BestSol.Position(1,2))]);
    
end


%% Results

figure;
xaxis=linspace(0,10,10);
%plot(BestCost,'LineWidth',2);
%semilogy(BestCost(:,1),'LineWidth',2);
%semilogy(BestCost(:,1),'LineWidth',2);
xlabel('Iteration');
ylabel('Best Cost');
plot(xaxis,BestCost(:,1),xaxis,BestCost(:,2))
grid on;


%%%%%%%%%%%%%%% Embedding %%%%%%
PSNR= Embedd(im,watermark,BestSol.Position)
