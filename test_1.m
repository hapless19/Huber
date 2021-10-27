
clc;clear all;close all;

%----------------------------------
addpath('\data\');
%----------------------------------
load init.txt;
initstate = init;%initial tracker

%----------------------------Set path
img_dir = dir('\data\*.png');
%---------------------------
img = imread(img_dir(1).name);
img = double(img(:,:,1));
alpha = 1;
% %--------------------------------------------------------
x = initstate(1);% x axis at the Top left corner
y = initstate(2);
w = initstate(3);% width of the rectangle
h = initstate(4);% height of the rectangle
pos1(1,:) = initstate;

p = pos1(1,:);  %david

opt = struct('numsample', 600, 'condenssig',1, 'posNum',10, 'negNum',30, 'graError',0.1, 'maxGraIteNum',15, ...
    'tmplsize', [16, 16], 'batchsize',1,'r',2, 'tmplNum', 2, 'nrml_threshold',1, 'gridSpacing',1, 'c',10, ...
    'patchSize',8, 'maxImSize',16, 'gamma', 0.15, 'ispool', 0, ...
    'lambda', 0.0001, 'pyramid', [1, 2], 'conf', 0.8, ...
    'bk_cons_mat', [1,1,0,0; 0,0,1,1; -0.8, 0.8,0,0; 0,0,-0.8,0.8]);

opt.affsig = [4, 4, 0.01, 0,0,0]; 
opt.posSig = [4,4,.01,.02,.002,.001];

frame =  double(img)/256;

param0 = [p(1), p(2), p(3)/opt.tmplsize(2), p(5), p(4)/p(3), 0];
param0 = affparam2mat(param0);
param.est = param0';

Dir = img_dir;
fNum = size(Dir,1)-3; 

paramRes = zeros(6, fNum);    
paramRes(:,1) = param.est;
paramRes_predict = zeros(4, fNum); 
paramRes_predict(:,1) = [0 0 0 0]';

[posSample, posPosition] = SelectPos(frame, opt.tmplsize, param.est, opt.posNum, opt.posSig);

for i = 1:size(posSample, 3)
     t1 = posSample(:,:,i);
     template(:,i) = t1(:);
end

tmpl.D = template;
 
imshow(uint8(img));
    rectangle('Position',initstate,'LineWidth',4,'EdgeColor','r');
    hold on;
    text(5, 18, strcat('#',num2str(1)), 'Color','y', 'FontWeight','bold', 'FontSize',20);
    set(gca,'position',[0 0 1 1]); 
    pause(0.00001); 
    hold off;
    
% %--------------------------------------------------------
 frame_index = 2;
    img = imread(img_dir(frame_index).name);
    imgSr = img;% imgSr is used for showing tracking results.
    img = double(img(:,:,1));
    frame = double(img)/256;
    
    f_frame = sprintf('%d',frame_index);
    fprintf('tracking frame: the %sth frame\n', f_frame); 
    
    n = opt.numsample;
    sz = opt.tmplsize;

    param.param = repmat(affparam2geom(param.est(:)), [1,n]);
    %prediction
    param.param = param.param + randn(6,n).*repmat(opt.affsig(:),[1,n]); 
    % show the predicted region by particle filtering
    boundary = zeros(4,n);
    tmp = zeros(4,1);
    tmp_conv = [1,0,-0.5,0; 1,0,0.5,0; 0,1,0,-0.5;0,1,0,0.5];
    for i=1:n
        tmp(1:2) = param.param(1:2,i);
        tmp(3) = param.param(3,i)*opt.tmplsize(2);
        tmp(4) = param.param(5,i)*tmp(3);    
        boundary(:,i) = tmp_conv*tmp;    
    end
    left = min(boundary(1,:));
    right = max(boundary(2,:));
    top = min(boundary(3,:));
    down = max(boundary(4,:));
    bkBox = [left, right, top, down];
    bkBox = round(bkBox);
    param.bkBox = bkBox;

    frm = frame;
    wimgs = warpimg(frm, affparam2mat(param.param), sz);
    
    X = wimgs;
    patchsize=opt.patchSize;
    gridspacing=opt.gridSpacing;
    num_samples = size(X,3);
    
    for i=1:num_samples 
        I = X(:,:,i);
        All_fea(:,i) = I(:);      
    end

    pos.rho=1.0;
    pos.tolerance=1e-4;
    pos.lambda1=0.1;
    pos.lambda2=0.1;
    pos.lambda3=0.1;
    pos.max_iter=100;

    [ZP, ZE,fun_value_difference,fun_value]=ADMM_l11_Huber(tmpl.D, All_fea,pos);
    W1 = ZP;
 
    eta_max	= -inf;
    for i=1:size(All_fea,2)
        W2 = W1(:,i);
        D_s = (All_fea(:,i) - tmpl.D*W2).^2;
        eta_1(i) = exp(-alpha*sqrt(sum(D_s)));
        if(eta_1(i)>eta_max)
            id_max	= i;
            c_max	= W2;      
            eta_max = eta_1(i);
        end
    end
    
    maxidx = id_max;
    param.est = affparam2mat(param.param(:,maxidx)); %
    param.wimg = wimgs(:,:,maxidx);
    param.obser = W1(:,maxidx);
    
    paramRes(:,frame_index) = param.est;
    paramRes_predict(:,frame_index) = param.bkBox;
    
    
    posPosition2 = paramRes(:,frame_index);
    M_p = [posPosition2(1,1) posPosition2(3,1) posPosition2(4,1); posPosition2(2,1) posPosition2(5,1) posPosition2(6,1)];
    w = opt.tmplsize(1);
    h = opt.tmplsize(2);
    corners = [ 1,-w/2,-h/2; 1,w/2,-h/2; 1,w/2,h/2; 1,-w/2,h/2; 1,-w/2,-h/2 ]';
    corners1 = M_p * corners;
    track_res(1, frame_index) = floor(corners1(1,1));
    track_res(2, frame_index) = floor(corners1(2,1));
    track_res(3, frame_index) = floor(corners1(1,2) - corners1(1,1));
    track_res(4, frame_index) = floor(corners1(2,3) - corners1(2,1));
 
    initstate = [track_res(1, frame_index) track_res(2, frame_index) track_res(3, frame_index) track_res(4, frame_index)];
    pos1(frame_index,:) = initstate;

    imshow(uint8(imgSr));
    rectangle('Position',initstate,'LineWidth',4,'EdgeColor','r');
    hold on;
    text(5, 18, strcat('#',num2str(frame_index)), 'Color','y', 'FontWeight','bold', 'FontSize',20);
    set(gca,'position',[0 0 1 1]); 
    pause(0.00001); 
    hold off;
