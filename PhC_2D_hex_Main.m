%% INITIALIZE OF HEXAGONAL AND HONEYCOMB LATTICE

% Developed by: Ibrahim Issah. 
% Msc. University of Eastern Finland
% Ref. Igor Photonic crystals
clear all
close all
clc



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DASHBOARD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c=2.99792458e8;

%% DASHBOARD 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%  Plotting parameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

band=1;               %% Band structures 
Field=1;              %% Plot field at M-Gamma-K
higher_order_gamma=1; %% Plot the field at the varied gamma_points
Afield = 0;           %% Additional Band to plot using the higher_order_gamma
Epsilon=0;            %% Crystal structure
scale_factor =20;     %% arrows length for quiver plot
autoscale    =0.5;     %% quiver


AAbs=0;               %% Plot abs(E)
RReal=1;              %% Plot real(E)
IImag=0;              %% Plot imag(E)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hex         = 0;      %% Actual hex
comb        = 0;      %% Actual honeycomb
triangle    = 0;      %% triangular lattice
square      = 1;      %% square lattice
rectangular = 0;      %% rectangular
rhombic     = 0;      %% rhombic lattice
hexagonal   = 0;      %% hexagonal lattice
ellipse     = 0;      %% ellipse lattice

TM=0;
TE=1;

Nx=128;         % number of points on the x grid % has to be a power of 2 (32,64,128,256,512,...)
Ny=128;         % number of points on the y grid % has to be a power of 2 (32,64,128,256,512,...)
NGx=10;        % number of harmonics % has to be 2 times -1 smaller than x
NGy=15;        % number of harmonics % has to be 2 times -1 smaller than y

Nkx=40;        % number of points on the k space for the dispersion
Nky=Nkx;       % number of points on the k space for the dispersion

nmodes=3;     % number of solutions needed 
nmodes_sub = 1; % number of solutions to neglect in the expansion
Np=1;          % number of period to plot for the Field

ep1 = 1.46^2;
ep2 = 0.17^2; 
n1 =sqrt(ep1);       %% optical index material 1
n2 = sqrt(ep2);      %% optical index material 2


NormUnits=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Building of the index Geometry %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if NormUnits==1
  L=1;
elseif NormUnits==0
  L=1e-6;  
end


Lx=L;
Ly=L*sqrt(3)/2;
x=linspace(-Lx/2, Lx/2, Nx);
y=linspace(-Ly/2, Ly/2, Ny);
a1=Lx*[1 0];
a2=Lx*[  1/2  sqrt(3)/2 ];

count=1;
for jj=0:Nx-1
  	for j=0:Ny-1
      AAA(count,:) = jj*a1/(Nx-1) + j*a2/(Ny-1) ;
      count=count+1;
   end
end

Xhex=reshape(AAA(:,1),Ny,Nx);
Yhex=reshape(AAA(:,2),Ny,Nx);

dx=Xhex(1,2)-Xhex(1,1);
dy=Yhex(2,1)-Yhex(1,1);


if (hex==1) && (comb==1) || (hex==0) && (comb==0)
    disp('Error: Select hexagonal lattice or honey-comb lattice')
    %break
end

if hex==1
    a=0.1;%0.495;
    idx1 =  ( (Xhex-Lx*3/4).^2 + (Yhex-Ly/2).^2 ) < (a*L)^2;
    idx2 =  ( (Xhex-Lx*3/4+Lx/2).^2 + (Yhex-Ly/2+Ly).^2 ) < (a*L)^2;
    idx3 =  ( (Xhex-Lx*3/4-Lx/2).^2 + (Yhex-Ly/2-Ly).^2 ) < (a*L)^2;
    idx4 =  ( (Xhex-Lx*3/4-Lx).^2 + (Yhex-Ly/2).^2 ) < (a*L)^2;
    idx5 =  ( (Xhex-Lx*3/4+Lx).^2 + (Yhex-Ly/2).^2 ) < (a*L)^2;
    
    idx=idx1+idx2+idx3+idx4+idx5;

    eps = idx*n2^2 + (1-idx)*n1^2 ;
end

if comb==1
    a=0.1;
    idx1a =  ( (Xhex-Lx*3/4).^2 + (Yhex-Ly/5).^2 ) < (a*L)^2;
    idx1b =  ( (Xhex-Lx*3/4+Lx/2).^2 + (Yhex-Ly/5-Lx/sqrt(3)+Ly).^2 ) < (a*L)^2;
    idx2a =  ( (Xhex-Lx*3/4).^2 + (Yhex-Ly/5-Lx/sqrt(3)).^2 ) < (a*L)^2;
    idx2b =  ( (Xhex-Lx*3/4-Lx/2).^2 + (Yhex-Ly/5-Ly).^2 ) < (a*L)^2;
    
    idx=idx1a+idx1b+idx2a+idx2b;

    eps = idx*n2^2 + (1-idx)*n1^2 ;
end 

if ellipse == 1 
    a =1; 
    Rx=4;       %% radius in the the x-direction of the ellipse [m]
    Ry=8;       %% radius in the the y-direction of the ellipse [m]
    x0=0;y0=0;   
    [XX, YY] = meshgrid(y, x);
    %% center positions of the ellipse [m]
    idx= ( abs(XX-x0)<(a*Lx) ) .* ( abs(YY-y0)<(a*Ly) ) ;
    idxa=((XX-x0)/Rx).^2 + ((YY-y0)/Ry).^2 < .001*a;
    idxaa= idxa.*idx;
    eps = idxaa*n2^2 + (1-idxaa)*n1^2 ;
end 

if hexagonal == 1 
    a = .45; 
    [XX,YY]=meshgrid(x,y);
    idx1=  (abs(XX)<a*sqrt(3)/2);
    idx2=(tan(pi/6)*XX+a>YY) .* (tan(pi/6)*XX-a<YY) .* (-tan(pi/6)*XX-a<YY) .* (-tan(pi/6)*XX+a>YY);
    idx=idx1.*idx2;
    eps = idx*n2^2 + (1-idx)*n1^2 ;
end 

if rhombic ==1 
    a =.35;
    [XX,YY]=meshgrid(x,y);
    idx1= (abs(XX)<a*sqrt(3)/2);
    idx2=(tan(pi/6)*XX+a>YY) .* (tan(pi/6)*XX-a<YY);
    idx=idx1.*idx2;
    idxa = (XX.^2 + YY.^2)>.01;
    idxaa = idx.*idxa; 
    eps = idxaa*n2^2 + (1-idxaa)*n1^2 ;
end

if square == 1 
    [XX,YY] = meshgrid(x,y);
    a=.1;
    idx = (XX.^2 + YY.^2) < (a*L)^2;
    eps = idx*n2^2 + (1-idx)*n1^2 ;
end 

if rectangular == 1
    x0=0;y0=0;
    a= .2; 
    [XX,YY] = meshgrid(x,y);
    idx= (abs(XX-x0)<(a*Lx)) .* ( abs(YY-y0)<(a*Ly) ) ;
    idxa = (XX.^2+ YY.^2) >(.1*a*L)^2;
    idxaa= idxa.*idx;
    eps = idxaa*n2^2 + (1-idxaa)*n1^2 ;
end 

if triangle == 1
a = 1; 
w   = 0.8 * a*L;
eps  = n1^2* ones(Nx,Ny);
h   = sqrt(w^2 - (w/2)^2);    %% Height of triangle
ny  = round(h/dy);
ny1 = 1 + floor((Ny - ny)/2);
ny2 = ny1 + ny - 1;
for ny = ny1 : ny2
    f = (ny - ny1 + 1)/(ny2 - ny1 + 1);
    nx = round(f*w/dx);
    nx1 = 1 + floor((Nx - nx)/2);
    nx2 = nx1 + nx - 1;
    eps(ny, nx1:nx2) = n2^2;
end 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Reciprocal lattice vectors %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NGx = 2*floor(NGx/2);           %% round to lower even number
NGy = 2*floor(NGy/2);           %% round to lower even number

b1=2*pi/Lx*[1  -sqrt(3)/3];
b2=2*pi/Lx*[0 2*sqrt(3)/3];

count=1;
GGG=[];

for jj=-NGx:NGx
for j=-NGy:NGy
    GGG(count,:)=jj*b1+j*b2;
    count=count+1;
end
end

Gxhex=reshape(GGG(:,1),2*NGy+1,2*NGx+1);
Gyhex=reshape(GGG(:,2),2*NGy+1,2*NGx+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Hexagonal Fourier Transform %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Gamma=1./eps;
f=Gamma;

for jj=1:length(Gxhex(1,:))
for j=1:length(Gyhex(:,1))
        whex = exp( -1i*(   Gxhex(1,jj) *(Xhex-Xhex(1))*(Nx-1)/Nx + ( Gyhex(j,jj) )*(Yhex-Yhex(1))*(Ny-1)/Ny));
        Ghex(j,jj) = sum(sum(f.*whex));
end
end

Gammak = Ghex*dx*dy/Lx/Ly ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Building of the reciproque lattice vector %%%% again %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

count=1;
GGG=[];
for jj=-NGx/2:NGx/2
    for j=-NGy/2:NGy/2
        GGG(count,:)=jj*b1+j*b2;
        count=count+1;
    end
end

Gxhex=reshape(GGG(:,1),NGy+1,NGx+1);
Gyhex=reshape(GGG(:,2),NGy+1,NGx+1);

NGx=length(Gxhex(1,:));
NGy=length(Gyhex(:,1));
NG=NGx*NGy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Building of k-space vector %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kx=linspace( 0 , pi/L , Nkx)*2/3;
ky=linspace( 0 , pi/Ly , Nky);

k=[
sort(kx,'descend')'    sort(ky,'descend')'
ky'*0                 ky'   
%kx'                   ky(end)+kx'*0

];
k(Nkx, :)=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% NOTHING TO CHANGE ANYMORE!!! %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (TE==1 && TM==1) || (TE==0 && TM==0)
  disp('Error: Select "TM" or "TE"');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Building first part of Hamitonian that is not depending on k %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HHH=zeros(NGy,NGx,NGy,NGx);

for ix=1:NGx
for jx=1:NGx
    for iy=1:NGy
    for jy=1:NGy
        HHH(iy,ix,jy,jx) = Gammak(iy-jy+NGy,ix-jx+NGx );
    end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(k(:,1))
  
  [psi,f0]=PhC2D_hex_PWE_f(Xhex,Yhex,Gxhex,Gyhex,k(i,:),HHH,nmodes,TE,TM);
  
  E(:,:,:,i)=psi;
  
  if NormUnits==1
    FF(:,i) = f0 * Lx / (2*pi);
  elseif NormUnits==0
    FF(:,i) = f0 * c / (2*pi) *1e-12;     % Convert in THz
    lambda(:,i)=2*pi./f0*1e6;             % Convert in wavelength (um)
  end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if AAbs==1
  EE=abs(E);
end
if RReal==1
  EE=real(E);
end
if IImag==1
  EE=imag(E);
end

if NormUnits==0
  Xhex=Xhex*1e6;
  Yhex=Yhex*1e6;
  Lx=Lx*1e6;
  Ly=Ly*1e6;
  k=k*1e-6;
end

if Field==1
    
    if TM==1 && TE==0
      figure('position',[100 50 1000 1000],'name','Ez');
    elseif TE==1 && TM==0
      figure('position',[100 50 1000 1000],'name','Exy');
    end
    colormap(jet)
    
    for ii=0:nmodes-1
        for i=1:Np
            for j=1:Np
                subplot(nmodes,3,1+3*ii)
                hold on
                pcolor(x+(i-1/2)*Lx , y+(j-1/2)*Ly , EE(:,:,ii+1,1)); 
                contour( x+(i-1/2)*Lx , y+(j-1/2)*Ly,abs(eps),1,'linewidth',2,'linecolor','w')
            end       
        end
        shading interp
        
        if RReal==1 || IImag==1
          caxis([-1 1])
        elseif AAbs==1
          caxis([0 1])
        end
        if NormUnits==1
          title(strcat('M: w=' , num2str(FF(1+ii,1), '%.2f') ))
          xlabel('x (a.u)')
          ylabel('y (a.u)')
        elseif NormUnits==0 
          title(strcat('M: \lambda=' , num2str(lambda(1+ii,1), '%.2f') , 'um' ))
          xlabel('x (um)')
          ylabel('y (um)')
        end
        xlim([0 1*Np*Lx])
        ylim([0  Np*Ly])
        
    end
    for ii=0:nmodes-1
        for i=1:Np
            for j=1:Np
                subplot(nmodes,3,2+3*ii)
                hold on
                pcolor(x+(i-1/2)*Lx , y+(j-1/2)*Ly , EE(:,:,ii+1,1*Nkx))
                %imagesc((0:size(Xhex, 1)-1)./size(Xhex, 1), (0:size(Yhex, 1)-1)./size(Yhex,1), EE(:, :,ii+1, 3))
                %pcolor( (0:size(Xhex, 1)-1)./size(Xhex, 1) , (0:size(Yhex, 1)-1)./size(Yhex,1) , EE(:,:,ii+1,1*Nkx) )       
                contour( x+(i-1/2)*Lx , y+(j-1/2)*Ly,abs(eps),1,'linewidth',2,'linecolor','w')
            end
        end
        shading interp
        %colorbar
        if RReal==1 || IImag==1
          caxis([-1 1])
        elseif AAbs==1
          caxis([0 1])
        end
        if NormUnits==1
          title(strcat('\Gamma: w=' , num2str(FF(1+ii,ceil(length(k)/2)), '%.2f') ))
          xlabel('x (a.u)')
          ylabel('y (a.u)')
        elseif NormUnits==0 
          title(strcat('\Gamma: \lambda=' , num2str(lambda(1+ii,ceil(length(k)/2)), '%.2f') , 'um' ))
          xlabel('x (um)')
          ylabel('y (um)')
        end
       xlim([0 1*Np*Lx])
       ylim([0     Np*Ly])
    end
    for ii=0:nmodes-1
        for i=1:Np
            for j=1:Np
                subplot(nmodes,3,3+3*ii)
                hold on
                pcolor(x+(i-1/2)*Lx , y+(j-1/2)*Ly , EE(:,:,ii+1,2*Nkx-1));
                contour( x+(i-1/2)*Lx , y+(j-1/2)*Ly,abs(eps),1,'linewidth',2,'linecolor','w');
             
            end
        end
        shading interp
        colorbar
        if RReal==1 || IImag==1
          caxis([-1 1])
        elseif AAbs==1
          caxis([0 1])
        end
        if NormUnits==1
          title(strcat('K: w=' , num2str(FF(1+ii,ceil(length(k)*2/2)), '%.2f') ));
          xlabel('x (a.u)');
          ylabel('y (a.u)');
        elseif NormUnits==0 
          title(strcat('K: \lambda=' , num2str(lambda(1+ii,ciel(length(k)*2/2)), '%.2f') , 'um' ));
          xlabel('x (um)')
          ylabel('y (um)')
        end
        xlim([0 1*Np*Lx]);
        ylim([0 Np*Ly]);
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Epsilon==1
  
  figure('position',[1100 50 500 400])
  subplot(111)
  hold on
     
  for i=1:Np
    for j=1:Np
        pcolor( Xhex+(i-1+(j-1)/2)*Lx,Yhex+(j-1)*Ly,real(eps));
    end
  end
  shading flat
  colormap(jet)
  c=colorbar;
  title(c,'Epsilon')
  if NormUnits==1
    xlabel('x (a.u)')
    ylabel('y (a.u)')
  elseif NormUnits==0 
    xlabel('x (um)')
    ylabel('y (um)')
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if band==1
    
    figure1= figure('position',[43 62 1194 643]);
    
    subplot(1,3,1)
    hold on; grid on;
    
    plot(0:length(k)-1,real(FF(1:nmodes,:))','o-','linewidth', 2);
    
    yscale=get(gca,'ylim');
    xlim([0 length(k)-1]);
    
    plot( [0/3*length(k)    0/3*length(k)] , yscale , 'k');
    plot( [1/2*length(k)    1/2*length(k)] , yscale , 'k');
    plot( [3/3*length(k)    3/3*length(k)] , yscale , 'k');
    
    text(0/3*length(k) , -0.05*yscale(2) , 'M');
    text(1/2*length(k) , -0.05*yscale(2) , '\Gamma');
    text(2/2*length(k) , -0.05*yscale(2) , 'K');
    %text(3/3*length(k) , -0.05*yscale(2) , ' M')
    %xlabel('k')
    set(gca,'xticklabel',[])
        
    if NormUnits==1
      ylabel('w (2\pi/Ltot)')
    elseif NormUnits==0 
      ylabel('f (THz)')
    end
    title(strcat('R/a=',num2str(a),'; n1=',num2str(n1,'%.2f'),'; n2=',num2str(n2,'%.2f')  ))
    subplot(1,3,2)  
    kkk=0:length(k)-1;
    my_index = kkk > 25 & kkk < 55; 
    plot(kkk(my_index), FF(2:nmodes-nmodes_sub, my_index),'o-', 'linewidth', 2); 
    title('Magnifide band of the lattice'); 
    colorbar
    grid on

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    subplot(1,3,3)
    hold on;grid on;
    LW=2;
    
    for i=1:nmodes
      plot3( k(:,1)*Lx/pi , k(:,2)*Ly/pi , real(FF(i,:))','o-');
    end
    set(gca,'xticklabel',[]);
    set(gca,'yticklabel',[]);
    xlabel('kx (pi/Lx)');
    ylabel('ky (pi/Ly)');
    if NormUnits==1
      zlabel('w (2\pi/Ltot)');
    elseif NormUnits==0 
      zlabel('f (THz)');
    end
    view(-7, 7);
    
    plot3( [-1 1]*2/3 , +[1 1] , [0 0] ,'b', 'linewidth',LW );
    plot3( [-1 1]*2/3 , -[1 1] , [0 0] ,'b', 'linewidth',LW );
    plot3( +[1 2]*2/3 , +[1 0] , [0 0] ,'b', 'linewidth',LW );
    plot3( +[1 2]*2/3 , -[1 0] , [0 0] ,'b', 'linewidth',LW );
    plot3( -[2 1]*2/3 , -[0 1] , [0 0] ,'b', 'linewidth',LW );
    plot3( -[2 1]*2/3 , +[0 1] , [0 0] ,'b', 'linewidth',LW );
    
    plot3( [-1 1]*2/3 , +[1 1] , [1 1]*max(real(FF(nmodes,:))) ,'b', 'linewidth',LW );
    plot3( [-1 1]*2/3 , -[1 1] , [1 1]*max(real(FF(nmodes,:))) ,'b', 'linewidth',LW );
    plot3( +[1 2]*2/3 , +[1 0] , [1 1]*max(real(FF(nmodes,:))) ,'b', 'linewidth',LW );
    plot3( +[1 2]*2/3 , -[1 0] , [1 1]*max(real(FF(nmodes,:))) ,'b', 'linewidth',LW );
    plot3( -[2 1]*2/3 , -[0 1] , [1 1]*max(real(FF(nmodes,:))) ,'b', 'linewidth',LW );
    plot3( -[2 1]*2/3 , +[0 1] , [1 1]*max(real(FF(nmodes,:))) ,'b', 'linewidth',LW );
    
    plot3( +[1 1]*2/3 , +[1 1] , [0 1]*max(real(FF(nmodes,:))) ,'b', 'linewidth',LW );
    plot3( +[1 1]*2/3 , -[1 1] , [0 1]*max(real(FF(nmodes,:))) ,'b', 'linewidth',LW );
    plot3( -[1 1]*2/3 , +[1 1] , [0 1]*max(real(FF(nmodes,:))) ,'b', 'linewidth',LW );
    plot3( -[1 1]*2/3 , -[1 1] , [0 1]*max(real(FF(nmodes,:))) ,'b', 'linewidth',LW );
    plot3( +2*[1 1]*2/3 ,0*[1 1] , [0 1]*max(real(FF(nmodes,:))) ,'b', 'linewidth',LW);
    plot3( -2*[1 1]*2/3 ,0*[1 1] , [0 1]*max(real(FF(nmodes,:))) ,'b', 'linewidth',LW);
    
    annotation(figure1,'textbox',...
    [0.216736040609137 0.440124416796268 0.0387631133671743 0.108864696734059],...
    'LineWidth',2,...
    'LineStyle','--',...
    'FitBoxToText','off',...
    'EdgeColor',[0 0 1]);

% Create arrow
annotation(figure1,'arrow',[0.254653130287648 0.409475465313029],...
    [0.55054432348367 0.91601866251944],'Color',[1 0 0],'LineWidth',2,...
    'LineStyle','--');

% Create arrow
annotation(figure1,'arrow',[0.25703125 0.41015625],...
    [0.441747572815534 0.116504854368932],'Color',[1 0 0],'LineWidth',2,...
    'LineStyle','--');

% Create textbox
annotation(figure1,'textbox',...
    [0.517187500000001 0.475351351351352 0.053125 0.0422297297297298],...
    'String',{'6 modes'},...
    'FontSize',12,...
    'EdgeColor',[1 1 1]);
end
colormap('jet'); 
if higher_order_gamma == 1
    figure4= figure('position',[43 62 1194 643]);
    for ii=1:nmodes-1
            for i=1:Np
                for j=1:Np
                    subplot(3,3, ii)
                    hold on
                    pcolor(x+(i-1/2)*Lx , y+(j-1/2)*Ly ,(EE(:,:,ii+1+Afield,1*Nkx)));  
                    [U, V] = curl(EE(:, :, ii+1, 1*Nkx), (EE(:, :, ii+1, 1*Nkx)));
                    h1 = quiver(x(1:scale_factor:end)+(i-1/2)*Lx, y(1:scale_factor:end)+(j-1/2)*Ly, U((1:scale_factor:end), (1:scale_factor:end)),V((1:scale_factor:end), (1:scale_factor:end)));
                    set(h1,'AutoScale','on', 'AutoScaleFactor',autoscale)
                    contour( x+(i-1/2)*Lx , y+(j-1/2)*Ly,abs(eps),1,'linewidth',2,'linecolor','w')
                end
            end
            shading interp
            %colorbar
            if RReal==1 || IImag==1
              caxis([-1 1])
            elseif AAbs==1
              caxis([0 1])
            end
            if NormUnits==1
              title(strcat('\Gamma: w=' , num2str(FF(1+ii,ceil(length(k)/2)), '%.2f') ))
              xlabel('x (a.u)')
              ylabel('y (a.u)')
            elseif NormUnits==0 
              title(strcat('\Gamma: \lambda=' , num2str(lambda(1+ii,ceil(length(k)/2)), '%.2f') , 'um' ))
              xlabel('x (um)')
              ylabel('y (um)')
            end
           xlim([0 1*Np*Lx])
           ylim([0  Np*Ly])
           axis equal tight
           axis square
           colorbar
           colormap(jet);
    end
end

%% Section for addition of Higher order modes

% Higher_modes = ones(size(EE(:, :, 1, Nkx))); 
% for ii = 1:nmodes-1
%      Higher_modes =Higher_modes.*EE(:, :, ii+1, Nkx);
% end
% figure(23)
% pcolor(Higher_modes./max(Higher_modes))
% shading interp
% colorbar 
% axis equal tight
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%