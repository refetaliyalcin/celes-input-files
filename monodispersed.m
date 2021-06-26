clc
clear all
close all
lamda=500;%nanometer wavelength
addpath(genpath('./src'))
ne_repeat=10;
n_sio2=1.5;
k_sio2=0;
% f_v=[0.01]
f_v = [0.01,0.02,0.05,0.1,0.15,0.2];
length_fv=length(f_v);

cyl_radius=3000;
cyl_height=2000;
cyl_vol=pi*cyl_height*cyl_radius^2;
rect_side=cyl_radius*2*2;
rect_height=cyl_height*2;
Vol_rect=rect_side*rect_side*rect_height;
radius = 150;
Vol_sphere=(4/3)*pi*radius^3;
% rng('shuffle')
% rng(1,'twister');
T=zeros(length_fv,ne_repeat);
R=zeros(length_fv,ne_repeat);
for i=4:length_fv
for tekrar=1:ne_repeat
rng('shuffle')
   

    number_of_spheres_in_rect = round(Vol_rect*f_v(i)/Vol_sphere);
    vol_spheres=Vol_sphere*number_of_spheres_in_rect;
    rect_fv_amac=f_v(i)
    sayac=1;
   
    x=rand()*rect_side-rect_side/2;
    y=rand()*rect_side-rect_side/2;
    z=rand()*rect_height-rect_height/2;
    data=[x,y,z,radius];

    i_sphere=2;

    tried=0;
    while i_sphere<number_of_spheres_in_rect
        x=rand()*rect_side-rect_side/2;
        y=rand()*rect_side-rect_side/2;
        z=rand()*rect_height-rect_height/2;
        cakisma=0;
        for j=1:size(data,1)
            if (data(j,1)-x)^2+(data(j,2)-y)^2+(data(j,3)-z)^2<(radius+data(j,4))^2
               cakisma=1;
               break
            end
        end
        if (cakisma==0)
           data=[data;x,y,z,radius];
           i_sphere=i_sphere+1;
           tried=0;
        end
        tried=tried+1;
        if tried>10^6
            disp('Failed to add, skip')
           i_sphere=i_sphere+1;
            tried=0;
        end
    end
    rect_fv_sonuc=size(data,1)*Vol_sphere/Vol_rect
    inside_cylinder=data(:,1).^2+data(:,2).^2<(cyl_radius).^2;
    data=data(inside_cylinder,:);
    inside_cylinder=abs(data(:,3))<(0.5*cyl_height);
    data=data(inside_cylinder,:);
    number_of_spheres_in_cyl=size(data,1);
    cyl_fv_sonuc=number_of_spheres_in_cyl*Vol_sphere/cyl_vol
   
    data=[data,n_sio2*ones(size(data,1),1), k_sio2*ones(size(data,1),1)];


    % initialize the CELES class instances

    % initialize particle class instance
    %   - positionArray:        Nx3 array (float) in [x,y,z] format
    %   - refractiveIndexArray: Nx1 array (complex) of complex refractive indices
    %   - radiusArray:          Nx1 array (float) of sphere radii
    particles = celes_particles('positionArray',        data(:,1:3), ...
                                'refractiveIndexArray', data(:,5)+1i*data(:,6), ...
                                'radiusArray',          data(:,4) ...
                                );

    % initialize initial field class instance
    %   - polarAngle:           scalar (float) polar angle of incoming beam/wave,
    %                           in radians. for Gaussian beams, only 0 or pi are
    %                           currently possible
    %   - azimuthalAngle:       scalar (float) azimuthal angle of incoming
    %                           beam/wave, in radians
    %   - polarization:         string (char) polarization of incoming beam/wave
    %                           ('TE' or 'TM')
    %   - beamWidth:            scalar (float) width of beam waist. use 0 or inf
    %                           for plane wave
    %   - focalPoint:           1x3 array (float) focal point
    initialField = celes_initialField('polarAngle',     0, ...
                                      'azimuthalAngle', 0, ...
                                      'polarization',   'TE', ...
                                      'beamWidth',      500, ...
                                      'focalPoint',     [0,0,0] ...
                                      );

    % initialize input class instance
    %   - wavelength:           scalar (float) vacuum wavelength, same unit as
    %                           particle positions and radii
    %   - mediumRefractiveIndex: scalar (complex) refractive index of environment
    %   - particles:            valid instance of celes_particles class
    %   - initialField:         valid instance of celes_initialField class
    input = celes_input('wavelength',                   lamda, ...
                        'mediumRefractiveIndex',        1, ...
                        'particles',                    particles, ...
                        'initialField',                 initialField ...
                        );

    % initialize preconditioner class instance
    %   - type:                 string (char) type of preconditioner (currently
    %                           only 'blockdiagonal' and 'none' available)
    %   - partitionEdgeSizes    1x3 array (float) edge size of partitioning cuboids
    %                           (applies to 'blockdiagonal' type only)
    precnd = celes_preconditioner('type',               'blockdiagonal', ...
                                  'partitionEdgeSizes', [500,500,500] ...
                                  );

    % initialize solver class instance
    %   - type:                 string (char) solver type (currently 'BiCGStab' or
    %                           'GMRES' are implemented)
    %   - tolerance:            scalar (float) target relative accuracy of solution
    %   - maxIter:              scalar (int) maximum number of iterations allowed
    %   - restart:              scalar (int) restart parameter (applies only to
    %                           GMRES solver)
    %   - preconditioner:       valid instance of celes_preconditioner class
    solver = celes_solver('type',                       'GMRES', ...
                          'tolerance',                  5e-4, ...
                          'maxIter',                    1000, ...
                          'restart',                    200, ...
                          'preconditioner',             precnd);

    % initialize numerics class instance
    %   - lmax:                 scalar (int) maximal expansion order of scattered
    %                           fields around particle center
    %   - polarAnglesArray:     1xN array (float) sampling of polar angles in the
    %                           plane wave patterns, in radians
    %   - azimuthalAnglesArray: sampling of azimuthal angles in the plane wave
    %                           patterns, in radians
    %   - gpuFlag:              scalar (bool) set to false if you experience GPU
    %                           memory problems at evaluation time (translation
    %                           operator always runs on GPU, though)
    %   - particleDistanceResolution: scalar (float) resolution of lookup table for
    %                           spherical Hankel function (same unit as wavelength)
    %   - solver:               valid instance of celes_solver class
    numerics = celes_numerics('lmax',                   3, ...
                              'polarAnglesArray',       0:pi/5e3:pi, ...
                              'azimuthalAnglesArray',   0:pi/1e2:2*pi, ...
                              'gpuFlag',                true, ...
                              'particleDistanceResolution', 1, ...
                              'solver',                 solver);

    % define a grid of points where the field will be evaulated
    [x,z] = meshgrid(-cyl_radius*1.5:25:cyl_radius*1.5, -cyl_height*1.5:25:cyl_height*1.5); y = zeros(size(x));
    % initialize output class instance
    %   - fieldPoints:          Nx3 array (float) points where to evaluate the
    %                           electric near field
    %   - fieldPointsArrayDims: 1x2 array (int) dimensions of the array, needed to
    %                           recompose the computed field as a n-by-m image
    output = celes_output('fieldPoints',                [x(:),y(:),z(:)], ...
                          'fieldPointsArrayDims',       size(x));

    % initialize simulation class instance
    %   - input:                valid instance of celes_input class
    %   - numerics:             valid instance of celes_input class
    %   - output:               valid instance of celes_output class
    simul = celes_simulation('input',                   input, ...
                             'numerics',                numerics, ...
                             'output',                  output);

%     plot results
%     display particles
%     figure('Name','Particle positions','NumberTitle','off');
%     plot_spheres(gca,simul.input.particles.positionArray, ...
%                      simul.input.particles.radiusArray, ...
%                      simul.input.particles.refractiveIndexArray)
% 
% 
%     error('stop')
    % run simulation
    simul.run;

    % evaluate transmitted and reflected power
    simul.evaluatePower;
    T(i,tekrar)=simul.output.totalFieldForwardPower/simul.output.initialFieldPower*100
    R(i,tekrar)=simul.output.totalFieldBackwardPower/simul.output.initialFieldPower*100
    sonuclar=[T(i,tekrar),R(i,tekrar),real_fv_cyl,cyl_radius];
     writematrix(sonuclar,'silicaplswork.txt','WriteMode','append')
    %dlmwrite(sonuclar,'WriteMode','append')
    % evaluate field at output.fieldPoints
    % simul.evaluateFields;


    %
    %  plot near field
    % figure('Name','Near-field cross-cut','NumberTitle','off');
    % plot_field(gca,simul,'abs E','Total field')
    % caxis([0,2])

    %  export animated gif
    % figure('Name','Animated near-field cross-cut','NumberTitle','off');
    % plot_field(gca,simul,'real Ey','Total field','Ey_total.gif')
end
end
sum(T,2)/ne_repeat
std(T)
