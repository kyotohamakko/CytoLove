function [ei, orient,AODF_F, ei2, orient2] = sdeconv(im,varargin)
shape = size(im);
orig_shape = shape;

show_kernel=false;


numits=50;
reg=1;
aniso=0.8;
wid=1;
lam=0.1/4;
L=6;


cpp=false;



gpu=false;



for k = 1:2:length(varargin),
        eval(sprintf('%s=varargin{k+1};',varargin{k}));
end;

L = L+1;

if cpp
    gpu=false;
end;


if show_kernel
    ker = single(fspecial('gaussian',[33 33],wid));
    ker = forwardOp(ker,L,lam,cpp);
    
    A_rn = getAODF(ker);
    disp(size(A_rn));
    figure; imagesc(squeeze(A_rn(1, :,:)));title("Rotate");
    
    
    %ker = sum(sum(ker,3),4);
    %figure(2);
%     subplot(1,2,1); imagesc(real(ker));
%     subplot(1,2,2); imagesc(imag(ker));
    %figure; imagesc(real(ker));
    %figure; imagesc(imag(ker));
    drawnow;
    return
end



im=single(im);

if ~exist('mask','var')
    if mod(shape(1),2) == 0,
        im(end+1,:) = 0;
    end;
    if mod(shape(2),2) == 0,
        im(:,end+1) = 0;
    end;
    shape = size(im);
else
    if mod(shape(1),2) == 0,
        im(end+1,:) = 0;
        mask(end+1,:) = 0;
    end;
    if mod(shape(2),2) == 0,
        im(:,end+1) = 0;
        mask(:,end+1) = 0;
    end;
    shape = size(im);

    mask = single(mask);
    ims = im;
    ims(not(mask(:))) = mean(ims(mask(:)>0));
    ims = real(ifft2(fft2(ims).*fft2(fftshift(fspecial('gaussian',shape,5)))));
    ims(mask(:)>0) = im(mask(:)>0);
    im = ims;
end


gauss = single(circshift(fspecial('gaussian',shape,wid),[1 1]));

if gpu
   gauss=gpuArray(gauss) ;
   im=gpuArray(im) ;
end


ftgauss = fft2(fftshift(gauss));
ftgauss2 = ftgauss.^2;





imsm = (real(ifft2(fft2(im).*ftgauss)));
% imsm = ifft2(fft2(im).*ftgauss);
% figure;imagesc(imsm);colormap gray;title("imsm");
Ay = forwardOp(imsm,L,lam,cpp);

x = single(zeros([shape,L,2]));

%if gpu
%   Ay=gpuArray(Ay) ;
%   x=gpuArray(x) ;
%   ftgauss2=gpuArray(ftgauss2) ;
%end

x = conjgrad(@(z)(theOp(z,L,lam,reg,aniso,ftgauss2,cpp)),Ay,x,numits);

if gpu
    x=gather(x);
end

[ei, orient] = getRidgeIndicator(x);
ei=ei(1:orig_shape(1),1:orig_shape(2));
orient=orient(1:orig_shape(1),1:orig_shape(2));
%if gpu
 %  ei=gather(ei);
 %  orient=gather(orient);
%end

AODF_F = getAODF(x);
AODF_F=AODF_F(:, 1:orig_shape(1),1:orig_shape(2));

if nargout>2
    [ei2, orient2] = getRidgeIndicator(-x);
    ei2=ei2(1:orig_shape(1),1:orig_shape(2));
    orient2=orient2(1:orig_shape(1),1:orig_shape(2));
 %   if gpu
  %     ei2=gather(ei2);
   %    orient2=gather(orient2);
    %end 
end

%[ei2,orient2] = getRidgeIndicator(-x) ;

%figure(1);
%subplot(2,2,1);
%imagesc(ei);
%subplot(2,2,2);
%imagesc(-ei2);
%subplot(2,2,3);
%imagesc(ei-ei2);
%subplot(2,2,4);
%imagesc(im);



%figure(3);
%q = backwardOp(enhance(x),L,lam,cpp);
%q = real(q(:,:,1));
%imagesc(q.*(q>0));


function y = theOp(x,L,lam,reg,aniso,ftgauss2,cpp)
    Ax = backwardOp(x,L,lam,cpp);
    Ax = ifft2(fft2(squeeze(Ax(:,:,1)+Ax(:,:,2))).*ftgauss2);
    AtAx = forwardOp(Ax,L,lam,cpp);

if cpp
    y = -reg*FCdiffusionProp(x,aniso) + AtAx;
else
    y= -reg*cat(4,sd_diffprop(squeeze(x(:,:,:,1)),aniso,1),sd_diffprop(squeeze(x(:,:,:,2)),aniso,-1))+ AtAx;
end



function q = enhance(f)

L = size(f,3);
d = cat(3,squeeze(f(:,:,:,1)),zeros(size(f,1),size(f,2),20),flipdim(squeeze(f(:,:,2:end,2)),3));
d(:,:,1) = d(:,:,1) + f(:,:,1,2);
d = permute(d,[3 1 2]);
fd = fft(d,L);

fd = fd.^1;

fd = ifft(fd);
d = ipermute(fd,[3 1 2]);
q(:,:,:,1) = d(:,:,1:L);
q(:,:,2:L,2) = flipdim(d(:,:,end-1:-1:(end-(L-1))),3);


function [ei, orient] = getRidgeIndicator(f)

d = cat(3,squeeze(f(:,:,:,1)),zeros(size(f,1),size(f,2),2000),flipdim(squeeze(f(:,:,2:end,2)),3));
d(:,:,1) = d(:,:,1) + f(:,:,1,2);

d = permute(d,[3 1 2]);
L = size(d,1);
fd = ((fft(d,L)));

id = 1:L;
theta = angle(exp(2*pi*1i*id/L))*180/pi;
bk_r_mid = size(theta(theta>0), 2);
[md_ id_] = max(real(fd(1:bk_r_mid,:,:)));
orient = squeeze(exp(2*pi*1i*id_/L));

[md id] = max(real(fd));

% orient = squeeze(exp(2*pi*1i*id/L));
md = squeeze(md);
ei = md;

function [AODF_F] = getAODF(f)

d = cat(3,squeeze(f(:,:,:,1)),zeros(size(f,1),size(f,2),200),flipdim(squeeze(f(:,:,2:end,2)),3));
d(:,:,1) = d(:,:,1) + f(:,:,1,2);

d = permute(d,[3 1 2]);
L = size(d,1);
fd = ((fft(d,L)));

AODF_F = real(fd);

% [md id] = max(real(fd));
% 
% orient = squeeze(exp(4*pi*1i*id/L));
% md = squeeze(md);
% ei = md;


function showridge(f,mask)
ei = getRidgeIndicator(f,mask) ;

sfigure(1); clf;
imagesc(ei);



function showquiv(f)
d = cat(3,squeeze(f(:,:,:,1)),zeros(size(f,1),size(f,2),20),flipdim(squeeze(f(:,:,2:end,2)),3));
d = permute(d,[3 1 2]);
%d(1,:,:) = 0;
figure(1);
clf;
L = size(d,1);

subplot(2,1,1);
fd = abs(imag(fft(d,L)));
[md id] = max(fd); 
md = squeeze(md);
dir = squeeze(exp(i*2*pi*id/L));
imagesc(squeeze(md)); hold on;

subplot(2,1,2);
fd = (real(fft(d,L)));
[md id] = max(fd); 
md = squeeze(md);
dir = squeeze(exp(i*pi*id/L));
imagesc(squeeze(md)); hold on;

%imagesc(squeeze(abs(d(end-1,:,:))));
%imagesc(squeeze(sum(abs(d).^2)))
%quiver(imag(dir).*md,real(dir).*md,0.5,'.r'); 
%quiver(-imag(dir).*md,-real(dir).*md,0.5,'.r'); hold off;

function [Ay] = forwardOp(img,A,freq,cpp)

order = 1;
typ = 0;

% p=freq/4;
K = 0:(A-1);
alpha1 = single(freq.^K ./  factorial(K));
alpha2 = single(((-freq).^K ./  factorial(K))');

%alpha1(2:2:end) = 0;
%alpha2(2:2:end) = 0;

if cpp
    Ay(:,:,:,1) = ComplexFilterHolo(img, alpha1, order, typ);
    Ay(:,:,:,2) = ComplexFilterHolo(img, alpha2, order, typ) ;
else
    Ay(:,:,:,1) =sd_filter(img, alpha1, 'filter');
    Ay(:,:,:,2) =sd_filter(img, alpha2, 'filter');
end




function [Ax] = backwardOp(x,A,freq,cpp)

order = 1;
typ = 0;

% p=freq/4;

K = 0:(A-1);
alpha1 = single(((-freq).^K ./  factorial(K))');
alpha2 = single(((freq).^K ./  factorial(K)));

%alpha1(2:2:end) = 0;
%alpha2(2:2:end) = 0;


if cpp
    Ax(:,:,1) = ComplexVoterHolo(x(:,:,:,1)+i*0.000000000001, alpha1, order, typ);
    Ax(:,:,2) = ComplexVoterHolo(x(:,:,:,2)+i*0.000000000001, alpha2, order, typ) ;
else
    Ax(:,:,1) = sd_filter(x(:,:,:,1)+1i*0.000000000001, alpha1, 'voter');
    Ax(:,:,2) = sd_filter(x(:,:,:,2)+1i*0.000000000001, alpha2, 'voter');
end





function [x] = conjgrad(Aop,b,x,its)
r=b-Aop(x);
p=r;
 rsold=r(:)'*r(:);
 
 for i=1:its
        Ap=Aop(p);
        alpha=(r(:)'*r(:))/(p(:)'*Ap(:));
        x=x+alpha*p;
        r=r-alpha*Ap;
        rsnew=r(:)'*r(:);
        
        if sqrt(rsnew)<1e-10
              break;
        end
        p=r+rsnew/rsold*p;
        rsold=rsnew;
end


   
   