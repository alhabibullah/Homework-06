%% Jacobi Method
clear all
close all
clc

Lx=2;
Ly=1;
n=10;
tolerance=1e-20;

x=linspace(0,Lx,n);
y=linspace(0,Ly,n);

dx=Lx/(n-1);
dy=Ly/(n-1);

%Initialization

u=zeros(n,n);
u_exact=zeros(n,n);
u(:,1)=0;
u(:,n)=y';
u(1,:)=u(2,:);
u(n,:)=u(n-1,:);
u_old=u;
err=1;
n_iter=1;
k=2*(dx^2+dy^2)/(dx^2*dy^2);

% %loop setup
% 
while err>tolerance
    for i=2:(n-1)
        for j=2:(n-1)
            u(i,j)=(1/k)*((u_old(i-1,j)+u_old(i+1,j))/dx^2+...
                (u_old(i,j-1)+u_old(i,j+1)/dy^2));
        end
    end
    u(1,:)=u(2,:);
    u(n,:)=u(n-1,:);
    err=max(max(abs(u-u_old)));
    n_iter=n_iter+1;
    u_old=u;

end
figure(1)
contourf(x,y,u,'ShowText','on');
colormap(jet);
xlabel('X axis')
ylabel('Y axis')
title_text=sprintf(['Jacobi method||' 'Total iterations=%d'],n_iter);
title(title_text);

%% Gauss Seidel Method
clear vars;

Lx=2;
Ly=1;
n=10;
tolerance=1e-20;

x=linspace(0,Lx,n);
y=linspace(0,Ly,n);

dx=Lx/(n-1);
dy=Ly/(n-1);

%Initialization

u=zeros(n,n);
u_exact=zeros(n,n);
u(:,1)=0;
u(:,n)=y';
u(1,:)=u(2,:);
u(n,:)=u(n-1,:);
u_old=u;
err=1;
n_iter=1;
k=2*(dx^2+dy^2)/(dx^2*dy^2);

% %loop setup
% 
while err>tolerance
    for i=2:(n-1)
        for j=2:(n-1)
            %Gauss Sidel setup using new prvious boundary value
            u(i,j)=(1/k)*((u(i-1,j)+u_old(i+1,j))/dx^2+...
                (u(i,j-1)+u_old(i,j+1)/dy^2));
        end
    end
    u(1,:)=u(2,:);
    u(n,:)=u(n-1,:);
    err=max(max(abs(u-u_old)));
    n_iter=n_iter+1;
    u_old=u;

end
figure(2)
contourf(x,y,u,'ShowText','on');
colormap(jet);
xlabel('X axis')
ylabel('Y axis')
title_text=sprintf(['Gauss Seidel method||' 'Total iterations=%d'],n_iter);
title(title_text);


%% SOR method
clear vars;
Lx=2;
Ly=1;
n=10;
tolerance=1e-20;

x=linspace(0,Lx,n);
y=linspace(0,Ly,n);

dx=Lx/(n-1);
dy=Ly/(n-1);

%Initialization

u=zeros(n,n);
u_exact=zeros(n,n);
u(:,1)=0;
u(:,n)=y';
u(1,:)=u(2,:);
u(n,:)=u(n-1,:);
u_old=u;
err=1;
n_iter=1;
k=2*(dx^2+dy^2)/(dx^2*dy^2);
u_1=zeros(n,n);
rel_factor=.9; % relaxation factor

% %loop setup
% 
while err>tolerance
    for i=2:(n-1)
        for j=2:(n-1)
            %SOR
            u_1(i,j)=(1/k)*((u(i-1,j)+u_old(i+1,j))/dx^2+...
                (u(i,j-1)+u_old(i,j+1)/dy^2));
            u(i,j)=u_old(i,j)*(1-rel_factor)+rel_factor*u_1(i,j);
        end
    end
    u(1,:)=u(2,:);
    u(n,:)=u(n-1,:);
    err=max(max(abs(u-u_old)));
    n_iter=n_iter+1;
    u_old=u;

end
figure(3)
contourf(x,y,u,'ShowText','on');
colormap(jet);
xlabel('X axis')
ylabel('Y axis')
title_text=sprintf(['SOR method||' 'Total iterations=%d'],n_iter);
title(title_text);





