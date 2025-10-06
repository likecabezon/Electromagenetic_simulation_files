close all
clear
clc

physical_constants;
unit = 1e-3; % all length in mm

feed.pos = 100;
feed.R = 50;

element.len = 930;
element.width = 30;
element.thicknes = 0.5;
element.cells = 4;

transformer.Spacing = 30;
transformer.extra_len = -60;

% size of the simulation box
SimBox = [2400 1500 2400];


f0 = 98e6; % center frequency
fc = 40e6; % 20 dB corner frequency
FDTD = InitFDTD('NrTS', 1000000 );
FDTD = SetGaussExcite( FDTD, f0, fc );
BC = {'MUR' 'MUR' 'MUR' 'MUR' 'MUR' 'MUR'}; % boundary conditions
FDTD = SetBoundaryCond( FDTD, BC );

max_res = c0 / (f0+fc) / unit / 20; % cell size: lambda/20
CSX = InitCSX();

%create fixed lines for the simulation box, substrate and port
mesh.x = [((-element.len/2) + 2*element.width + transformer.Spacing) ((-element.len/2) + element.width + transformer.Spacing) -SimBox(1)/2 SimBox(1)/2 -element.len/2 ((-element.len/2) + element.width) ((element.len/2) - element.width) element.len/2 -feed.pos (-feed.pos - element.width)];
mesh.x = SmoothMeshLines( mesh.x, max_res, 1.4); % create a smooth mesh between specified fixed mesh lines

mesh.y = [-SimBox(2)/2 SimBox(2)/2 linspace(-element.thicknes/2 , element.thicknes/2, element.cells)];
mesh.y = SmoothMeshLines( mesh.y, max_res, 1.4 );

%create fixed lines for the simulation box and given number of lines inside the substrate
mesh.z = [transformer.extra_len ((-element.len/2) + 2*element.width + transformer.Spacing) ((-element.len/2) + element.width + transformer.Spacing) -SimBox(3)/2 -element.len/2 ((-element.len/2) + element.width) ((element.len/2) - element.width) element.len/2 SimBox(3)/2 ];
mesh.z = SmoothMeshLines( mesh.z, max_res, 1.4 );

CSX = DefineRectGrid( CSX, unit, mesh );


%% create elements
CSX = AddMetal( CSX, 'arm' ); % create a perfect electric conductor (PEC)
start = [-element.len/2 -element.thicknes/2 -element.len/2];
stop  = [((-element.len/2) + element.width)  element.thicknes/2 element.len/2];
CSX = AddBox(CSX,'arm',10,start,stop); % add a box-primitive to the metal property 'patch'

start = [-element.len/2 -element.thicknes/2 ((element.len/2) - element.width)];
stop  = [element.len/2  element.thicknes/2 element.len/2];
CSX = AddBox(CSX,'arm',10,start,stop);

start = [((element.len/2) - element.width) -element.thicknes/2 -element.len/2];
stop  = [element.len/2  element.thicknes/2 element.len/2];
CSX = AddBox(CSX,'arm',10,start,stop);

start = [-element.len/2 -element.thicknes/2 -element.len/2];
stop  = [element.len/2  element.thicknes/2 ((-element.len/2) + element.width)];
CSX = AddBox(CSX,'arm',10,start,stop);


% matching network transformer

start = [-element.width/2 -element.thicknes/2 ((-element.len/2) + element.width)];
stop  = [element.width/2  element.thicknes/2 ((-element.len/2) + element.width + transformer.Spacing)];
CSX = AddBox(CSX,'arm',10,start,stop);

start = [((-element.len/2) + element.width + transformer.Spacing) -element.thicknes/2 ((-element.len/2) + element.width + transformer.Spacing)];
stop  = [element.width/2  element.thicknes/2 ((-element.len/2) + 2*element.width + transformer.Spacing)];
CSX = AddBox(CSX,'arm',10,start,stop);

start = [((-element.len/2) + element.width + transformer.Spacing) -element.thicknes/2 ((-element.len/2) + element.width + transformer.Spacing)];
stop  = [((-element.len/2) + 2*element.width + transformer.Spacing)  element.thicknes/2 transformer.extra_len];
CSX = AddBox(CSX,'arm',10,start,stop);

%% add lumped ports
start = [-feed.pos  -element.thicknes/2 ((-element.len/2) + element.width)];
stop  = [((-feed.pos) - element.width)  element.thicknes/2 ((-element.len/2) + element.width + transformer.Spacing)];
[CSX port] = AddLumpedPort(CSX, 5 ,1 ,feed.R, start, stop, [0 0 1], true);

SimBox = SimBox - max_res * 4; %reduced SimBox size for nf2ff box
[CSX nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', -SimBox/2, SimBox/2);

%% prepare simulation folder
Sim_Path = 'tmp';
Sim_CSX = 'metal_sheet_square_loop.xml';

[status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
[status, message, messageid] = mkdir( Sim_Path ); % create empty simulation folder

%% write openEMS compatible xml-file
WriteOpenEMS( [Sim_Path '/' Sim_CSX], FDTD, CSX );

%% show the structure
CSXGeomPlot( [Sim_Path '/' Sim_CSX] );

%% run openEMS
RunOpenEMS( Sim_Path, Sim_CSX );

%% postprocessing & do the plots
freq = linspace( f0-fc, f0+fc, 501 );
port = calcPort(port, Sim_Path, freq);

Zin = port.uf.tot ./ port.if.tot;
s11 = port.uf.ref ./ port.uf.inc;
P_in = 0.5 * port.uf.inc .* conj( port.if.inc ); % antenna feed power

% plot best archibable reflection coefficient
figure
best_ref_coef = (Zin - real(Zin))./(Zin + real(Zin));
plot(freq/1e6, 20*log10(abs(best_ref_coef)));
hold on
grid on
title( 'best_ref_coef S_{11}' );
xlabel( 'frequency f / MHz' );
ylabel( 'reflection coefficient |S_{11}|' );

% plot feed point impedance
figure
plot( freq/1e6, real(Zin), 'k-', 'Linewidth', 2 );
hold on
grid on
plot( freq/1e6, imag(Zin), 'r--', 'Linewidth', 2 );
title( 'feed point impedance' );
xlabel( 'frequency f / MHz' );
ylabel( 'impedance Z_{in} / Ohm' );
legend( 'real', 'imag' );

% plot reflection coefficient S11
figure
plot( freq/1e6, 20*log10(abs(s11)), 'k-', 'Linewidth', 2 );
grid on
title( 'reflection coefficient S_{11}' );
xlabel( 'frequency f / MHz' );
ylabel( 'reflection coefficient |S_{11}|' );

drawnow

%% NFFF contour plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find resonance frequency from s11
f_res_ind = find(s11==min(s11));
f_res = freq(f_res_ind);

% calculate the far field at phi=0 degrees and at phi=90 degrees
disp( 'calculating far field at phi=[0 90] deg...' );

nf2ff = CalcNF2FF(nf2ff, Sim_Path, f_res, [-180:2:180]*pi/180, [0 90]*pi/180);

% display power and directivity
disp( ['radiated power: Prad = ' num2str(nf2ff.Prad) ' Watt']);
disp( ['directivity: Dmax = ' num2str(nf2ff.Dmax) ' (' num2str(10*log10(nf2ff.Dmax)) ' dBi)'] );
disp( ['efficiency: nu_rad = ' num2str(100*nf2ff.Prad./real(P_in(f_res_ind))) ' %']);

% normalized directivity as polar plot
figure
polarFF(nf2ff,'xaxis','theta','param',[1 2],'normalize',1)

% log-scale directivity plot
figure
plotFFdB(nf2ff,'xaxis','theta','param',[1 2])
% conventional plot approach
% plot( nf2ff.theta*180/pi, 20*log10(nf2ff.E_norm{1}/max(nf2ff.E_norm{1}(:)))+10*log10(nf2ff.Dmax));

drawnow

%%
disp( 'calculating 3D far field pattern and dumping to vtk (use Paraview to visualize)...' );
thetaRange = (0:2:180);
phiRange = (0:2:360) - 180;
nf2ff = CalcNF2FF(nf2ff, Sim_Path, f_res, thetaRange*pi/180, phiRange*pi/180,'Verbose',1,'Outfile','3D_Pattern.h5');

figure
plotFF3D(nf2ff,'logscale',-20);


E_far_normalized = nf2ff.E_norm{1} / max(nf2ff.E_norm{1}(:)) * nf2ff.Dmax;
DumpFF2VTK([Sim_Path '/3D_Pattern.vtk'],E_far_normalized,thetaRange,phiRange,'scale',1e-3);
