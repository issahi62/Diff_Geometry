%% INITIALIZE OF SQUARE AND TRIANGULAR LATTICE

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
Field=0;              %% Plot field at M-Gamma-K
higher_order_gamma=1; %% Plot the field at the varied gamma_points
Epsilon=1;            %% Crystal structure

triangle     =1;       %% triangular lattice
square      = 0;      %% square lattice
rectangular = 0;      %% rectangular
rhombic     = 0;      %% rhombic lattice
hexagonal   = 0;      %% hexagonal lattice
ellipse     = 0;      %% ellipse lattice

AAbs=0;               %% Plot abs(E)
RReal=1;              %% Plot real(E)
IImag=0;              %% Plot imag(E)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TM=1;
TE=0;

Nx=128;         % number of points on the x grid % has to be a power of 2 (32,64,128,256,512,...)
Ny=128;         % number of points on the y grid % has to be a power of 2 (32,64,128,256,512,...)
NGx=10;        % number of harmonics % has to be 2 times -1 smaller than x
NGy=11;        % number of harmonics % has to be 2 times -1 smaller than y

Nkx=30;        % number of points on the k space for the dispersion
Nky=Nkx;       % number of points on the k space for the dispersion

nmodes=25;      % number of solutions asked
nmodes_sub = 3; % number of solutions to neglect in the expansion
Np=1;           % number of period to plot for the Field

n1 =1;         %% optical index material 1
n2 = 1.52; %% optical index material 2

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
Ly=L;
x=linspace(-Lx/2, Lx/2, Nx);
y=linspace(-Ly/2, Ly/2, Ny);
dx=x(2)-x(1);
dy=y(2)-y(1);

if (ellipse==1 && hexagonal==1 && square ==1 && rectangular ==1 && rhombic == 1 && triangle ==1) || (ellipse==0 && hexagonal==0 && square ==0 && rectangular ==0 && rhombic == 0 && triangle ==0 )
  disp('Error: Select "one lattice"')
  %break
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
a=0.1;
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
%%%%%%%%%%%%%%%%%%%%%% Building Epsilon in Fourier space %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NGx = 2*floor(NGx/2);           %% round to lower even number
NGy = 2*floor(NGy/2);           %% round to lower even number

Gamma=1./eps;
Gammak = fftshift(fft2(Gamma))*dx*dy/Lx/Ly;
Gammak = Gammak(Ny/2-NGy+1:Ny/2+NGy+1 , Nx/2-NGx+1:Nx/2+NGx+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Reciprocal lattice vectors %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Gx = (-NGx/2:NGx/2)'*2*pi/Lx;
Gy = (-NGy/2:NGy/2)'*2*pi/Ly;

NGx=length(Gx);
NGy=length(Gy);
NG=NGx*NGy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Building of k-space vector %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kx=linspace( 0 , pi/Lx , Nkx);
ky=linspace( 0 , pi/Ly , Nky);

k=[
sort(kx,'descend')'    sort(ky,'descend')'
kx'                    kx'*0
%ky'*0+kx(end)          ky'   

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
  disp('Error: Select "TM" or "TE"')
  %break
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
  
  [psi,f0]=PhC2D_sq_PWE_f(x,y,Gx,Gy,k(i,:),HHH,nmodes,TE,TM);
  
  E(:,:,:,i)=psi;
  
  if NormUnits==1
    FF(:,i) = f0 * Lx / (2*pi);%/sqrt(5);
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
  x=x*1e6;
  y=y*1e6;
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
                pcolor( x+(i-1/2)*Lx , y+(j-1/2)*Ly , EE(:,:,ii+1,1) )
                %contour( (x+(i-1/2)*Lx), y+(j-1/2)*Ly ,abs(eps),1,'linewidth',2,'linecolor','w')
            end       
        end
        shading flat
        
        %colorbar
        if RReal==1 || IImag==1
          caxis([-1 1])
        elseif AAbs==1
          caxis([0 1])
        end
        if NormUnits==1
          title(strcat('M: w=' , num2str(FF(1+ii,1), '%.2f') ))
          xlabel('x (norm. units)')
          ylabel('y (norm. units)')
        elseif NormUnits==0 
          title(strcat('M: \lambda=' , num2str(lambda(1+ii,1), '%.2f') , 'um' ))
          xlabel('x (um)')
          ylabel('y (um)')
        end
        xlim([0 Np*Lx])
        ylim([0 Np*Ly])
    end
    for ii=0:nmodes-1
        for i=1:Np
            for j=1:Np
                subplot(nmodes,3,2+3*ii)
                hold on
                pcolor( x+(i-1/2)*Lx , y+(j-1/2)*Ly , EE(:,:,ii+1,1*Nkx) )
                %contour( (x+(i-1/2)*Lx), y+(j-1/2)*Ly ,abs(eps),1,'linewidth',2,'linecolor','w')
            end
        end
        shading flat
        %colorbar
        if RReal==1 || IImag==1
          caxis([-1 1])
        elseif AAbs==1
          caxis([0 1])
        end
        if NormUnits==1
          title(strcat('\Gamma: w=' , num2str(FF(1+ii,ceil(length(k)/2)), '%.2f') ))
          xlabel('x (norm. units)')
          ylabel('y (norm. units)')
        elseif NormUnits==0 
          title(strcat('\Gamma: \lambda=' , num2str(lambda(1+ii,ceil(length(k)/2)), '%.2f') , 'um' ))
          xlabel('x (um)')
          ylabel('y (um)')
        end
        
        
        xlim([0 Np*Lx])
        ylim([0 Np*Ly])
       
    end
    for ii=0:nmodes-1
        for i=1:Np
            for j=1:Np
                subplot(nmodes,3,3+3*ii)
                hold on
                pcolor( x+(i-1/2)*Lx , y+(j-1/2)*Ly , EE(:,:,ii+1,2*Nkx-1) )
                %contour( (x+(i-1/2)*Lx), y+(j-1/2)*Ly ,abs(eps),1,'linewidth',2,'linecolor','w')
            end
        end
        shading flat
        %colorbar
        if RReal==1 || IImag==1
          caxis([-1 1])
        elseif AAbs==1
          caxis([0 1])
        end
        if NormUnits==1
          title(strcat('X: w=' , num2str(FF(1+ii,length(k)*2/2), '%.2f') ))
          xlabel('x (norm. units)')
          ylabel('y (norm. units)')
        elseif NormUnits==0 
          title(strcat('X: \lambda=' , num2str(lambda(1+ii,length(k)*2/2), '%.2f') , 'um' ))
          xlabel('x (um)')
          ylabel('y (um)')
        end
        xlim([0 Np*Lx])
        ylim([0 Np*Ly])
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Epsilon==1
  
  figure('WindowStyle', 'normal', 'position',[1100 50 500 400])
  subplot(111)
  hold on
     
  for i=1:Np
    for j=1:Np
        pcolor( x+(i-1/2)*Lx , y+(j-1/2)*Ly , real(eps) )
    end
  end
  shading flat
  colormap(jet)
  c=colorbar;
  title(c,'Epsilon')
  if NormUnits==1
    xlabel('x (norm. units)')
    ylabel('y (norm. units)')
  elseif NormUnits==0 
    xlabel('x (um)')
    ylabel('y (um)')
  end
  %axis equal

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if band==1
    
    figure1= figure('WindowStyle', 'normal','position',[50 62 913 666]);
    subplot(1,3,1)
    hold on; grid on;
    
    plot(0:length(k)-1,real(FF(1:nmodes,:))','o-', 'linewidth', 2)
   
    yscale=get(gca,'ylim');
    xlim([0 length(k)-1])
    
    plot( [0/3*length(k)    0/3*length(k)] , yscale , 'k')
    plot( [1/2*length(k)    1/2*length(k)] , yscale , 'k')
    plot( [3/3*length(k)    3/3*length(k)] , yscale , 'k')
    
    text(0/3*length(k) , -0.05*yscale(2) , ' M')
    text(1/2*length(k) , -0.05*yscale(2) , ' \Gamma'     )
    text(2/2*length(k) , -0.05*yscale(2) , ' X'     )
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
    grid on; 
    kkk=0:length(k)-1;
    my_index = kkk > 20 & kkk < 39; 
    plot(kkk(my_index), FF(2:nmodes-nmodes_sub, my_index),'o-', 'linewidth', 2); 
    title('Magnifide band of the lattice'); 
    colorbar
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    subplot(1,3,3)
    hold on;grid on;
    LW=2.5;
    
    for i=1:nmodes
      plot3( k(:,1)*Lx/pi , k(:,2)*Ly/pi , real(FF(i,:))','o-')
    end
    set(gca,'xticklabel',[])
    set(gca,'yticklabel',[])
    xlabel('kx (pi/Lx)')
    ylabel('ky (pi/Ly)')
    if NormUnits==1
      zlabel('w (2\pi/Ltot)')
    elseif NormUnits==0 
      zlabel('f (THz)')
    end
    xlim([-1 1]*max(k(:,1))*Lx/pi)
    ylim([-1 1]*max(k(:,2))*Ly/pi)
    view(-65,5)
    
    plot3( [-1 1] , +[1 1] , [0 0] ,'b', 'linewidth',LW )
    plot3( [-1 1] , -[1 1] , [0 0] ,'b', 'linewidth',LW )
    plot3( +[1 1] , [-1 1] , [0 0] ,'b', 'linewidth',LW )
    plot3( -[1 1] , [-1 1] , [0 0] ,'b', 'linewidth',LW )
    
    plot3( [-1 1] , +[1 1] , [1 1]*max(real(FF(nmodes,:))) ,'b', 'linewidth',LW )
    plot3( [-1 1] , -[1 1] , [1 1]*max(real(FF(nmodes,:))) ,'b', 'linewidth',LW )
    plot3( +[1 1] , [-1 1] , [1 1]*max(real(FF(nmodes,:))) ,'b', 'linewidth',LW )
    plot3( -[1 1] , [-1 1] , [1 1]*max(real(FF(nmodes,:))) ,'b', 'linewidth',LW )
    
    plot3( +[1 1] , +[1 1] , [0 1]*max(real(FF(nmodes,:))) ,'b', 'linewidth',LW )
    plot3( -[1 1] , -[1 1] , [0 1]*max(real(FF(nmodes,:))) ,'b', 'linewidth',LW )
    plot3( +[1 1] , -[1 1] , [0 1]*max(real(FF(nmodes,:))) ,'b', 'linewidth',LW )
    plot3( -[1 1] , +[1 1] , [0 1]*max(real(FF(nmodes,:))) ,'b', 'linewidth',LW )
    
    annotation(figure1,'rectangle',...
    [0.216772179627601 0.578078078078078 0.0417163198247535 0.0870870870870869],...
    'Color',[1 0 0],...
    'LineWidth',2,...
    'LineStyle','--');

% Create arrow
annotation(figure1,'arrow',[0.255202628696605 0.411829134720701],...
    [0.666666666666667 0.924924924924925],'Color',[0 0 1],'LineWidth',2,...
    'LineStyle','--');

% Create arrow
annotation(figure1,'arrow',[0.25859375 0.4125],...
    [0.578796561604584 0.107449856733524],'Color',[0 0 1],'LineWidth',2,...
    'LineStyle','--');
% Create textbox
annotation(figure1,'textbox',...
    [0.438526526560789 0.495495495495495 0.0592409638554218 0.027027027027027],...
    'String',{'4 modes'},...
    'LineWidth',2,...
    'LineStyle','--',...
    'HorizontalAlignment','center',...
    'FontSize',10,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);
end

if higher_order_gamma == 1
    figure4= figure('WindowStyle', 'normal','position',[43 62 1194 643]);
    for ii=1:9
        for i=1:Np
            for j=1:Np
                subplot(3,3, ii)
                hold on
                pcolor( x+(i-1/2)*Lx , y+(j-1/2)*Ly , EE(:,:,ii+16,1*Nkx) )
                %contour( (x+(i-1/2)*Lx), y+(j-1/2)*Ly ,abs(eps),1,'linewidth',2,'linecolor','w')
            end
        end
        shading flat
        %colorbar
        if RReal==1 || IImag==1
          caxis([-1 1])
        elseif AAbs==1
          caxis([0 1])
        end
        if NormUnits==1
          title(strcat('\Gamma: w=' , num2str(FF(1+ii,ceil(length(k)/2)), '%.2f') ))
          xlabel('x (norm. units)')
          ylabel('y (norm. units)')
        elseif NormUnits==0 
          title(strcat('\Gamma: \lambda=' , num2str(lambda(1+ii,ceil(length(k)/2)), '%.2f') , 'um' ))
          xlabel('x (um)')
          ylabel('y (um)')
        end
        xlim([0 Np*Lx])
        ylim([0 Np*Ly])
        axis equal tight
        colorbar
        colormap(jet);
    end
end

%% Section for addition of Higher order modes
% Higher_modes = ones(size(EE(:, :, 1, Nkx))); 
% for ii = 1:nmodes-1
%      Higher_modes =Higher_modes.*EE(:, :, ii+1, Nkx);
% end
% figure(23) $x^2+e^{\pi i}$
% pcolor(Higher_modes./max(Higher_modes))
% shading interp
% colorbar 
% axis equal tight

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%