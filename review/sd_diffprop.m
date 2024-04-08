function y=sd_diffprop(x,coeff,sg)

filter_mode='circular';

degree=size(x,3);

fac=coeff(1)/2;

y=zeros(size(x),'like',x);


for i=0:(degree-1)
    y(:,:,i+1)=y(:,:,i+1)+lap(squeeze(x(:,:,i+1)),filter_mode);
    
    if (i<(degree-2))
        y(:,:,i+1)=y(:,:,i+1)+fac*deriv_zz(x(:,:,i+1+2),-sg,filter_mode);
    end
    if i>1
        y(:,:,i+1)=y(:,:,i+1)+fac*deriv_zz(x(:,:,i+1-2),sg,filter_mode);
    else
        y(:,:,i+1)=y(:,:,i+1)+fac*deriv_zz(x(:,:,2-i+1),sg,filter_mode);
    end
end


function x=lap(x,filter_mode)

lap_kernel=...
           [0,1,0;
           1,-4,1;
           0,1,0];

x=imfilter(x,lap_kernel,filter_mode);



function x=deriv_zz(x,con,filter_mode)


zz_kernel=...
           [1i,-1,-1i;
           1,0,1;
           -1i,-1,1i];        


w_kernel=...
           [0.5,1,0.5;
           1,1,1;
           0.5,1,0.5];       
       

zz_kernel=zz_kernel.*w_kernel;
       
if ~(con>0)
    zz_kernel=conj(zz_kernel);
end
       
x=  imfilter(x,zz_kernel,filter_mode);       
       


