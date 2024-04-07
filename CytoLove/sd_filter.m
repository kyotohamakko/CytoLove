function y=sd_filter(x,coeff,type)

order=1;
%filter_mode='replicate';
filter_mode='circular';


degree=numel(coeff);

shape=size(x);
if numel(shape)==3
    shape=shape(1:2);
end

G{1}=zeros(shape,'like',x);
G{2}=zeros(shape,'like',x);

sign=1;        
if size(coeff,1)>1
    sign=-1;
end



in=1;
out=2;

switch type
    case 'filter'
        J=0:degree-1;
        test=@(j)(j>0);
        y=zeros([shape,degree],'like',x);
        %fprintf('allocating [%d %d %d] (%d GB)\n',[shape,degree],ceil(prod([shape,degree])*8/1024/1024));
        G{1}=x;
    case 'voter'
        J=degree-1:-1:0;    
        test=@(j)(j<(degree-1));
end


%gpu=strcmp(class(x),'gpuArray');



for j=J
    if test(j)
        G{out}=sd_deriv(G{in},sign,order,filter_mode);
        [in,out]=swap(in,out);
    end
    switch type
       case 'filter'
            y(:,:,j+1)=coeff(j+1)*G{in};
       case 'voter' 
            if j==0
                y=G{in}+coeff(j+1)*squeeze(x(:,:,j+1));
            else
                G{in}=G{in}+coeff(j+1)*squeeze(x(:,:,j+1));
            end;
    end
end



function [y,x]=swap(x,y)


function x=sd_deriv(x,con,order,filter_mode)


switch (order)
    case 2
    complex_derivatives=...
           [0,0,(-1i)*0.2,0,0;
           0,0,1i,0,0;
           (-1)*0.2,1 ,0,-1, (1)*0.2 ;
           0,0,-1i,0,0;
           0,0,(1i)*0.2,0,0;]*1.5;
       
    otherwise
    complex_derivatives=...
           [0,-1i,0;
           (-1),0,1;
           0,1i,0];

        
end;
       
if ~(con>0)
    complex_derivatives=conj(complex_derivatives);
end
       


x=imfilter(x,complex_derivatives,filter_mode);
        
