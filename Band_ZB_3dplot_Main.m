%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% last update 4January2022, lne %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Here, you have to choose your material among the following %%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Material='GaAs'; SpinOrbit=1;
Material='Si'; SpinOrbit=1;
%Material='Ge'; SpinOrbit=0;    % there are no spin-orbit parameter for Ge!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Library
ExtractParameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nk=30;              %% number of k-points in each direction

if strcmp(M{1},'GaAs')==1
    kmax=0.1;
end
if strcmp(M{1},'Si')==1
    kmax=0.8;
end
if strcmp(M{1},'Ge')==1
    kmax=0.6;
end

% band are computed in the first positive part from 0 to kmax and then concatenate 
kx = linspace(0,kmax,Nk);
ky = linspace(0,kmax,Nk);
kz = linspace(0,kmax,Nk);

[Ky,Kx,Kz] = meshgrid(ky,kx,kz);
k = [Kx(:) Ky(:) Kz(:)] * 2*pi/a;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
for ii=1:length(k(:,1))
    if SpinOrbit==0
        Ek = TightBinding_f(M{4},a,k(ii,:),SpinOrbit);
        Ekk(ii) = Ek(5); % only the first conduction band Energy is stored
    elseif SpinOrbit==1
        Ek = TightBinding_f(M{5},a,k(ii,:),SpinOrbit);
        Ekk(ii) = Ek(9); % only the first conduction band Energy is stored
    end
end
toc

Ekkk = reshape(Ekk,[length(kx) length(ky) length(kz)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here, I concatenate Ekkk into a bigger matrix E thanks to the symmetry
% Doing so, the code is 8 times faster

kx = linspace(-kmax,kmax,2*length(kx));
ky = linspace(-kmax,kmax,2*length(ky));
kz = linspace(-kmax,kmax,2*length(kz));

[Ky,Kx,Kz] = meshgrid(ky,kx,kz);

E=[
Ekkk(end:-1:1,end:-1:1,:) Ekkk(end:-1:1,:,:)
Ekkk(:,end:-1:1,:)        Ekkk
];

E = cat( 3 , E(:,:,end:-1:1) , E );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% real lattice construction
L=2;
shift=L/2;

c = [
0 0 0 ; L 0 0 ; L L 0 ; 0 L 0 ; 0 0 0
0 0 L ; L 0 L ; L L L ; 0 L L ; 0 0 L
L 0 L ; L 0 0 ; L 0 L ; L L L ; L L 0 ; L L L ; 0 L L ; 0 L 0 
] - shift;

cc = [ L/2 L/2 L/2  ] - shift;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reciproque lattice construction

bb1 = [ 0 L/4 L/2 ; 0 L/2 3*L/4 ; 0 3*L/4 L/2 ; 0 L/2 1*L/4 ; 0 L/4 L/2 ] - shift;
bb2 = [ L L/4 L/2 ; L L/2 3*L/4 ; L 3*L/4 L/2 ; L L/2 1*L/4 ; L L/4 L/2 ] - shift;
bb3 = [ L/4 0 L/2 ; L/2 0 3*L/4 ; 3*L/4 0 L/2 ; L/2 0 1*L/4 ; L/4 0 L/2 ] - shift;
bb4 = [ L/4 L L/2 ; L/2 L 3*L/4 ; 3*L/4 L L/2 ; L/2 L 1*L/4 ; L/4 L L/2 ] - shift;
bb5 = [ L/4 L/2 0 ; L/2 3*L/4 0 ; 3*L/4 L/2 0 ; L/2 1*L/4 0 ; L/4 L/2 0 ] - shift;
bb6 = [ L/4 L/2 L ; L/2 3*L/4 L ; 3*L/4 L/2 L ; L/2 1*L/4 L ; L/4 L/2 L ] - shift;

bbb1  = [ 0 L/4 L/2 ; L/4 0 L/2 ] - shift;
bbb2  = [ 0 L/2 L/4 ; L/4 L/2 0 ] - shift;
bbb3  = [ L/2 L/4 0 ; L/2 0 L/4 ] - shift;

bbb4  = [ 0 L/2 3*L/4 ; L/4 L/2 L ] - shift;
bbb5  = [ L/2 L/4 L ; L/2 0 3*L/4 ] - shift;

bbb6  = [ 3*L/4 L/2 L ; L L/2 3*L/4 ] - shift;
bbb7  = [ 3*L/4 0 L/2 ; L L/4 L/2 ] - shift;
bbb8  = [ 3*L/4 L/2 0 ; L L/2 L/4 ] - shift;
bbb9  = [ L 3*L/4 L/2 ; 3*L/4 L L/2 ] - shift;
bbb10 = [ L/2 3*L/4 L ; L/2 L 3*L/4 ] - shift;
bbb11 = [ L/2 3*L/4 0 ; L/2 L L/4 ] - shift;
bbb12 = [ L/4 L L/2 ; 0 3*L/4 L/2 ] - shift;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('position',[100 100 1000 700],'color','w')
FS=15;
LW=2;
MS=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,1,1,'fontsize',FS)
hold on;grid on;box on;

if strcmp(M{1},'Si')==1
   v=1.3; % value a bit higher than the material bandgap in eV
end
if strcmp(M{1},'GaAs')==1
   v=1.7; % value a bit higher than the material bandgap in eV
end
if strcmp(M{1},'Ge')==1
   v=0.79; % value a bit higher than the material bandgap in eV
end

p = patch(isosurface(kx,+ky,+kz,E,v));
isonormals(kx,+ky,+kz,E, p)

set(p, 'facecolor','red','edgecolor','none')

camlight;
lighting gouraud;
light ('Position', [-1 -1 0]);

% xlabel('kx')
% ylabel('ky')
% zlabel('kz')

axis equal
xlim([-1 1])
ylim([-1 1])
zlim([-1 1])
view(-30,10)
title(strcat(M{1},'; isosurface=',num2str(v),'eV'));
set(gca,'xticklabel',[],'yticklabel',[],'zticklabel',[])
%set(gca,'Box','off','Visible','off')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plot3( c(:,1) , c(:,2) , c(:,3) , 'k.-','markersize',MS,'linewidth',LW)
%plot3( cc(:,1) , cc(:,2) , cc(:,3) , 'k.','markersize',MS,'linewidth',LW)
%plot3( b(:,1) , b(:,2) , b(:,3) , 'b.','markersize',20,'linewidth',2)

plot3( bb1(:,1) , bb1(:,2) , bb1(:,3) , 'b','linewidth',LW)
plot3( bb2(:,1) , bb2(:,2) , bb2(:,3) , 'b','linewidth',LW)
plot3( bb3(:,1) , bb3(:,2) , bb3(:,3) , 'b','linewidth',LW)
plot3( bb4(:,1) , bb4(:,2) , bb4(:,3) , 'b','linewidth',LW)
plot3( bb5(:,1) , bb5(:,2) , bb5(:,3) , 'b','linewidth',LW)
plot3( bb6(:,1) , bb6(:,2) , bb6(:,3) , 'b','linewidth',LW)

plot3( bbb1(:,1) , bbb1(:,2) , bbb1(:,3) , 'b','linewidth',LW)
plot3( bbb2(:,1) , bbb2(:,2) , bbb2(:,3) , 'b','linewidth',LW)
plot3( bbb3(:,1) , bbb3(:,2) , bbb3(:,3) , 'b','linewidth',LW)
plot3( bbb4(:,1) , bbb4(:,2) , bbb4(:,3) , 'b','linewidth',LW)
plot3( bbb5(:,1) , bbb5(:,2) , bbb5(:,3) , 'b','linewidth',LW)
plot3( bbb6(:,1) , bbb6(:,2) , bbb6(:,3) , 'b','linewidth',LW)
plot3( bbb7(:,1) , bbb7(:,2) , bbb7(:,3) , 'b','linewidth',LW)
plot3( bbb8(:,1) , bbb8(:,2) , bbb8(:,3) , 'b','linewidth',LW)
plot3( bbb9(:,1) , bbb9(:,2) , bbb9(:,3) , 'b','linewidth',LW)
plot3( bbb10(:,1) , bbb10(:,2) , bbb10(:,3) , 'b','linewidth',LW)
plot3( bbb11(:,1) , bbb11(:,2) , bbb11(:,3) , 'b','linewidth',LW)
plot3( bbb12(:,1) , bbb12(:,2) , bbb12(:,3) , 'b','linewidth',LW)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
