cover_dir = '/data-x/g12/zhangjiansong/stable_diffusion_G/gen1000_7.5000'
mean_dir 
var_dir
payload
stego_dir
lambda
coef
if ~exist(stego_dir,'dir')
	mkdir(stego_dir);
end
% load cover image
parfor i=1:10000
    coverPath = fullfile(cover_dir, [num2str(i),'.pgm']);
    stegoPath = fullfile(stego_dir, [num2str(i),'.pgm']);
    meanPath = fullfile(mean_dir, [num2str(i),'.npy']);
    varPath = fullfile(var_dir, [num2str(i),'.npy']);
    cover = double(imread(coverPath));
    mean = readNPY(meanPath);
    var = readNPY(varPath);
%% Run default embedding
    stego = HILLEmbed(cover, payload,var,mean,lambda,coef);
    
    stego(stego>255)=255;
    stego(stego<0)=0;
    stego = uint8(stego);
    imwrite(stego,stegoPath);
end
fprintf(' - DONE');


