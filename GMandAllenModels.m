
% basic RD Gierer – Meinhardt (hereafter GM) pattern formation code
% on a torus (plane with periodic boundary conditions at each edge
% Ref:  Miyazawa, Seita, Michitoshi Okamoto, and Shigeru Kondo. "Blending
% of animal colour patterns by hybridization." Nature communications 1
% (2010): 66.
% u = activator concentration
% v = inhibitor concentration
% diffusion code adapted from:  Diffusion in 1D and 2D
% version 1.0.0.0 (3.44 KB) by Suraj Shankar
% https://www.mathworks.com/matlabcentral/fileexchange/38088-diffusion-in-1d-and-2d
% Simulating the 2-D Diffusion equation by the Finite Difference
...Method
    % Numerical scheme used is a first order upwind in time and a second
...order central difference in space (Implicit and Explicit)
% parameters
% R*A = basal activator growth rate
% R*B = activator decay rate
% R = overall factor multipling growth and decay terms in both activator
% and inhibitor
% in this code, there is no separate inhibitor decay rate parameter
% Du = activator diffusion constant; must be << Dv for pattern formation
% Dv = inhibitor diffusion constant
% C = carrying capacity -- when = 0, reduces to the basic GM model; when >
% 0, limits the activator growth rate to R/C
% min_pattern_change : stop when the pattern changes by < this between the
% Ntimeview plots (I found values 0.00001 to 0.000001 worked pretty well,
% but use trial and error to check this for your cases
% -------------------------------------------------------------------------
% Example values: (default values for parameters not given)
% white spots on a black background: A = 0, B = 0.85, C = 0 and An=Bn=Cn=1
% labyrinthine: A = 0.05, B = 1.5, C = 0.5and An=Bn=Cn=1
% black spots on a white background: A = 0, B = 01, C = 0.36 and An=Bn=Cn=1
% -------------------------------------------------------------------------
% initialization
clear; close all;
% change directory to store results in
cd(uigetdir);
%         
% example values:
% set size (number of elements) of spatial grid on which we do calculations
% 7 = 128 <--smaller value used in fish hybrid paper
% 8 = 256; 9 = 512; 10 = 1024 <-- nice figure
% use powers of two for later processing
% min_pattern_change : stop when pattern changes by < this
% I found that this value 0.0001 worked well for a labyrinthine pattern
% example values from Miyazawa et al. 2010:
% Fig. 1:  R = 1.0, A = 0.08, B = 1.5, Du = 0.5, Dv = 10.0 ( Fig. 1 ) and
% Fig. 3:  R = 0.2, A = 0.08, B = 1.5, D u = 0.5 Du = 1.0 Dv = 20.0
% -------------------------------------------------------------------------
% input pattern formation specifications



list = {'GM Model','Allen Cat Model','Allen Snake Model'}; 
[indx,tf] = listdlg('ListString',list);


if ismember(1,indx)
    
    prompt = {'min. activator growth rate (A):',... %answer 1
        'max. activator growth rate (A):',...       %answer 2
        'number A values:',...                      %answer 3
        'min. inhibitor growth rate (B):',...       %answer 4
        'max. inhibitor growth rate (B):',...       %answer 5
        'number B values:',...                      %answer 6
        'min. carrying capacity (C):',...           %answer 7
        'max. carrying capacity (C):',...           %answer 8
        'number C values:',...                      %answer 9
        };
    title = 'Pattern formation specifications';
    dims = [1 35];              % input field specifications
    definput = {'0','0.2','1','0.85','2.0','3','0.359','1','1'};   % default values
    answer = inputdlg(prompt,title,dims,definput);
    
    prompt2 = {'R factor:',...                       %answer 10
        'activator diffusion constant:',...         %answer 11
        'inhibitor diffusion constant:',...         %answer 12
        'minimum percent change in pattern:',...    %answer 13
        'timesteps per pattern simulation:',...     %answer 14
        'number of grids per edge (power of 2)',... %answer 15
        'pattern edge size',...                     %answer 16
        'random number generator seed (0=random):',... % answer 17
        'display patterns (Y or N):',...            %answer 18
        'number of loops:',...                      %answer 19
        };
    title = 'Pattern formation specifications';
    dims = [1 35];              % input field specifications
    definput2 = {'0.2','1.0','20','0','30000','2^7','2','0','Y','1'};   % default values
    answer2 = inputdlg(prompt2,title,dims,definput2);
    
    [Amin,statusflag] = str2num(answer{1});
    [Amax,statusflag] = str2num(answer{2});
    [An,statusflag] = str2num(answer{3});
    [Bmin,statusflag] = str2num(answer{4});
    [Bmax,statusflag] = str2num(answer{5});
    [Bn,statusflag] = str2num(answer{6});
    [Cmin,statusflag] = str2num(answer{7});
    [Cmax,statusflag] = str2num(answer{8});
    [Cn,statusflag] = str2num(answer{9});
    [R,statusflag] = str2num(answer2{1});
    [Du,statusflag] = str2num(answer2{2});
    [Dv,statusflag] = str2num(answer2{3});
    [min_pattern_change,statusflag] = str2num(answer2{4});
    [Ntimesteps,statusflag] = str2num(answer2{5}); % total number of tiemsteps simulated:
    % Miyazawa et al 2010 used 20K to 80K
    [ngrid,statusflag] = str2num(answer2{6});
    [spatial_scale,statusflag] = str2num(answer2{7});
    [rand_seed,statusflag] = str2num(answer2{8});
    draw_pattern = char(answer2(9));
    [num_loops,statusflag] = str2num(answer2{10});
    % ----------------------------------
    % -------------------------------------------------------------------------
    % time parameters
    %num_loops = 2; % number of times to repeat each simulation
    Ntimeview = round(Ntimesteps*0.01);  % display pattern every 1% runtime
    dt=0.01;                             % timestep
    
    binary_pattern = 0;                % display binary final image grid (0 = grayscale)
    % -------------------------------------------------------------------------
    show_percent_change = 1;           % display percentage pattern change
    show_percent_change = 0;
    
    % Initialize array of A, B, C values to eventually save with PSM values for
    % each image
    %Initialize storCount to index into ABCValues
    storCount = 0;
    ABCValues = zeros(An*Bn*Cn*num_loops,3);
    
    %%
    
    % create arrays to loop over
    if An > 1
        DeltaA = (Amax - Amin)/(An-1);
        Aarray = Amin:DeltaA:Amax;
    else
        Aarray = Amin;
    end
    if Bn > 1
        DeltaB = (Bmax - Bmin)/(Bn-1);
        Barray = Bmin:DeltaB:Bmax;
    else
        Barray = Bmin;
    end
    if Cn > 1
        DeltaC = (Cmax - Cmin)/(Cn-1);
        Carray = Cmin:DeltaC:Cmax;
    else
        Carray = Cmin;
    end
    
    % set up spatial grid
    nx=ngrid;                        %Number of steps in space(x)
    ny=ngrid;
    %Number of steps in space(y)
    dx=nx/(nx-1);                       %Width of space step(x)
    dy=ny/(ny-1);                       %Width of space step(y)
    x=0:dx:nx;                          %x grid points
    y=0:dy:ny;                          %y grid points
    % scale x and y to range [0,spatial_scale]
    x = spatial_scale*x/max(x);y = spatial_scale*y/max(y);
    u=zeros(nx,ny);                  %Preallocating u -- saves runtime
    un=zeros(nx,ny);                 %Preallocating un
    v=zeros(nx,ny);                  %Preallocating u -- saves runtime
    vn=zeros(nx,ny);                 %Preallocating un
    % -------------------------------------------------------------------------
    filename = cell(An,Bn,Cn);      % create a cell array to hold filenames
    % create graphics window
    for ii = 1:An
        for jj = 1:Bn
            for kk = 1:Cn
                iter_num = 0;
                while iter_num < num_loops
                    iter_num = iter_num + 1;
                    fprintf('analyzing the %d of %d A, %d of %d B,%d of %d C values\n',ii,An,jj,Bn,kk,Cn);
                    A = Aarray(ii); B = Barray(jj); C = Carray(kk);
                    if draw_pattern=='Y'
                        fig = figure('Name',strcat('activat. conc. (white=high): A=',num2str(A),';B=',num2str(B),';C=',num2str(C))','Position',[300,300,512,512]);
                        % make sure figure is square for saving
                        fig.PaperSize = [8.5 8.5];
                        fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 8.5 8.5];
                    end
                    % -------------------------------------------------------------
                    % start of time loop
                    %Initial Conditions for u and v concentrations
                    if rand_seed == 0           % random, based on time
                        rng('shuffle');
                    else                        % user defined for reproducibility
                        rng(rand_seed);
                    end
                    % "seed" (initialize)the random number generator
                    % so the simulations are repeatable; for later testing try out with
                    % different seeds (40 used in original paper)
                    %rng(10) is interesting notes Rebeckah
                    
                    %random initial conditions for u and v (t = 0 values)
                    % create each as an n^2 vector and reshape into n x n matrix
                    u=rand(ngrid*ngrid,1); u = reshape(u,ngrid,ngrid);
                    v=rand(ngrid*ngrid,1); v = reshape(v,ngrid,ngrid);
                    
                    % Calculating the spatial field index
                    % 1 and nx, ny = ngrid will be set by periodic boundary conditions
                    i=2:nx-1;
                    j=2:ny-1;
                    ij = 1:ngrid;
                    for t=0:Ntimesteps
                        uold = u;
                        vold = v;
                        un = u;               % un is a copy of the values at the current timestep
                        vn = v;
                        % u will become the updated value at the next timestep
                        % just growth and decay terms
                        un = un+R*(un .* un ./ (vn .* (1.0 + C .* un .* un)) + A - B .* un)*dt;
                        vn = vn+R*(un .* un - vn)*dt;
                        % diffusion for all but edges
                        u(i,j)=un(i,j)+...
                            (Du*dt*(un(i+1,j)-2*un(i,j)+un(i-1,j))/(dx*dx))+...
                            (Du*dt*(un(i,j+1)-2*un(i,j)+un(i,j-1))/(dy*dy));
                        % periodic boundary conditions
                        % activator
                        % diffusion at right = left edge (periodic bc's)
                        u(i,1)=un(i,1)+...
                            (Du*dt*(un(i+1,1)-2*un(i,1)+un(i-1,1))/(dx*dx))+...
                            (Du*dt*(un(i,2)-2*un(i,1)+un(i,ngrid))/(dy*dy));
                        % diffusion at top = bottom edge (periodic bc's)
                        u(1,j)=un(1,j)+...
                            (Du*dt*(un(2,j)-2*un(1,j)+un(ngrid,j))/(dx*dx))+...
                            (Du*dt*(un(1,j+1)-2*un(1,j)+un(1,j-1))/(dy*dy));
                        % deal with corners (1,1),(1,ngrid),(ngrid,1),(ngrid,ngrid)
                        u(1,1)=un(1,1)+...
                            (Du*dt*(un(2,1)-2*un(1,1)+un(ngrid,1))/(dx*dx))+...
                            (Du*dt*(un(1,2)-2*un(1,1)+un(1,ngrid))/(dy*dy));
                        u(1,ngrid)=u(1,1);u(ngrid,1)=u(1,1);u(ngrid,ngrid)=u(1,1);
                        % copy updated values to other edge
                        u(ngrid,:) = u(1,:);
                        u(:,ngrid) = u(:,1);
                        % inhibitor
                        % diffusion for all but edges
                        v(i,j)=vn(i,j)+...
                            (Dv*dt*(vn(i+1,j)-2*vn(i,j)+vn(i-1,j))/(dx*dx))+...
                            (Dv*dt*(vn(i,j+1)-2*vn(i,j)+vn(i,j-1))/(dy*dy));
                        % periodic boundary conditions
                        % diffusion at left edge
                        v(i,1)=vn(i,1)+...
                            (Dv*dt*(vn(i+1,1)-2*vn(i,1)+vn(i-1,1))/(dx*dx))+...
                            (Dv*dt*(vn(i,2)-2*vn(i,1)+vn(i,ngrid))/(dy*dy));
                        % diffusion at top
                        v(1,j)=vn(1,j)+...
                            (Dv*dt*(vn(2,j)-2*vn(1,j)+vn(ngrid,j))/(dx*dx))+...
                            (Dv*dt*(vn(1,j+1)-2*vn(1,j)+vn(1,j-1))/(dy*dy));
                        % deal with corners (1,1),(1,ngrid),(ngrid,1),(ngrid,ngrid)
                        v(1,1)=vn(1,1)+...
                            (Dv*dt*(vn(2,1)-2*vn(1,1)+vn(ngrid,1))/(dx*dx))+...
                            (Dv*dt*(vn(1,2)-2*vn(1,1)+vn(1,ngrid))/(dy*dy));
                        v(1,ngrid)=v(1,1);v(ngrid,1)=v(1,1);v(ngrid,ngrid)=v(1,1);
                        % copy updated values to other edge
                        v(ngrid,:) = v(1,:);
                        v(:,ngrid) = v(:,1);
                        % draw activator and inhibitor patterns
                        if draw_pattern=='Y'
                            if (mod(t,Ntimeview)==0)  % only draw every Ntimeview timesteps
                                % shows the image with u scaled to the full
                                % grayscale range--note it may look lighter when
                                % saved for this reason.
                                patternfig = image([0,ngrid],[0,ngrid],u,'CDataMapping','scaled');
                                colormap('gray');
                                %             colorbar;
                                %                         text(0,-20,strcat('percent time steps = ',num2str(t/Ntimesteps)));
                                drawnow;
                                refreshdata(patternfig);
                            end
                            if (mod(t,Ntimeview)==0)  % only draw every Ntimeview timesteps
                                if t > 2*Ntimeview
                                    diff = abs(u - uold);
                                    if show_percent_change
                                        fprintf('percent activator pattern change = %f; goal = %f\n',mean2(diff(:)),min_pattern_change * mean2(u));
                                    end
                                    if mean2(diff(:)) < min_pattern_change * mean2(u)
                                        break;
                                    end
                                end
                            end
                        end
                    end
                    % end of time loop
                    % -------------------------------------------------------------
                    % save pattern figure as binary
                    % saved as .tif, .png, and .jpg to compare pattern measures
                    % before and after conversion. All but one of these
                    % should be commented out unless deugging
                    img = u;
                    filename{ii,jj,kk,iter_num} = strcat('A_',num2str(A),'_B_',num2str(B),'_C_',num2str(C),'_Du_',num2str(Du),'_Dv_',num2str(Dv),'_iter_',num2str(iter_num),'.tif');
                    imwrite(img,filename{ii,jj,kk,iter_num});  % save version w/o white border
                    filename{ii,jj,kk,iter_num} = strcat('A_',num2str(A),'_B_',num2str(B),'_C_',num2str(C),'_Du_',num2str(Du),'_Dv_',num2str(Dv),'_iter_',num2str(iter_num),'.jpg');
                    imwrite(img,filename{ii,jj,kk,iter_num});  % save version w/o white border
                    filename{ii,jj,kk,iter_num} = strcat('A_',num2str(A),'_B_',num2str(B),'_C_',num2str(C),'_Du_',num2str(Du),'_Dv_',num2str(Dv),'_iter_',num2str(iter_num),'.png');
                    imwrite(img,filename{ii,jj,kk,iter_num});  % save version w/o white border
                    % save each image in an array to be used in making a montage
                    % (square array of images) -- note these are reversed in order
                    % so they display in the correct order in the montage)
                    image_montage_array{ii,jj,kk,iter_num} = img;
                    close;
                    % Increment storCount, store A, B, C values into ABCValues at
                    % index storCount
                    storCount = storCount + 1;
                    ABCValues(storCount,1) = A;
                    ABCValues(storCount,2) = B;
                    ABCValues(storCount,3) = C;
                    
                end
            end
        end
    end
    % =========================================================================
    % if desired, compute the distance in pattern measurements space between a
    % sample pattern and all the patterns generated
    answer = 'Yes';
    switch answer
        case 'Yes'
            numFiles = An*Bn*Cn*num_loops;                    % number of images
            k = 1;                                  % create single image array
            for ii = 1:An
                for jj = 1:Bn
                    for kk = 1:Cn
                        for  indind = 1:num_loops
                            image_cell_array{k} = image_montage_array{ii,jj,kk,indind};
                            k = k + 1;
                        end
                    end
                end
            end
            % -------------------------------------------------------------------------
            % compute the threshold for each image
            for k = 1:(numFiles)
                gray_threshold(k) = graythresh(image_cell_array{k});
                gray_threshold_comp(k) = graythresh(imcomplement(image_cell_array{k}));
            end
            binary_threshold = median(gray_threshold);  % this should help avoid enhancing contrast where there is none
            binary_threshold_comp = median(gray_threshold_comp);
            % -------------------------------------------------------------------------
            % compute all pattern measurement space coordinates
            %First turns image into greyscale, then binary
            %Then finds area of regions and the mean area. Uses mean area to filter out
            %small spots.
            %Filters out insignificant areas
            %Finds eccentricity of remaining regions and other measures
            %Saves in vector
            for k = 1:(numFiles)
                % make color images into grayscale
                [rows, columns, numberOfColorChannels] = size(image_cell_array{k});
                if numberOfColorChannels > 1                % convert to grayscale
                    img = rgb2gray(image_cell_array{k});
                else
                    img = image_cell_array{k};
                end
                imgaussfilt(img);           % blur to reduce noise and pixelation
                % compute fraction pigmentation
                fraction_pigmented(k,1) = mean2(img)/256;% 0-0.5: mostly black; > 0.5 mostly white
                % do the following line only if the images are natural photographs (not
                % simulated images) *%*
                img = adapthisteq(img);            % use adaptive histogram equalization to account for unequal illumination
                BW=im2bw(img,binary_threshold);    % turn grayscale image into binary (0's and 1's = black and white)
                % now 0 = black = pigmented and 1 = white = unpigmented, because we are
                % interested in the formation of pigmented regions so they are the
                % natural choice for measuring properties of the pattern
                
                BW = bwareafilt(BW,[3 Inf]);        % for natural photographs only %*% remove single pixel objects %*% black or white?  Do for both
                BW = imcomplement(BW); % now 1 = pigmented and 0 = unpigmented
                BW = bwareafilt(BW,[3 Inf]);        % remove single pixels %*%
                
                % now compute the fraction pigmented as the percent that's black
                num_pixels = numel(BW(:));          % total number of pixels in the images
                fraction_pigmented(k,1) = sum(BW(:))/num_pixels;
                %----------------------------------------------------------------------
                % find all white objects and their attributes
                s = regionprops(BW,'Eccentricity','ConvexArea','Centroid',...
                    'Area','Perimeter','EquivDiam','MajorAxis');
                BW_img{k} = BW;
                %----------------------------------------------------------------------
                % area-weighted measure of how compact the objects are on average
                % also known as isoperimetric quotient or Polsby-Popper compactness
                Area = cat(1, s.Area);
                num_objects = numel(Area);
                if num_objects > 0
                    obj_num_meas(k,1) = 1/sqrt(num_objects);
                else
                    obj_num_meas(k,1) = 1; % one uniform area
                end
                p_Area = Area/sum(Area);          % fractional area weighting
                % correct for any NaN values due to zero total area -- replace with
                % 0, for no objects
                if isnan(p_Area)
                    p_Area = 0;
                end
                %----------------------------------------------------------------------
                % area-weighted mean Polsby-Popper compactness:  range [0,1]
                Perimeter = cat(1, s.Perimeter);
                compactness(k,1) = 4*pi*sum(p_Area.*Area./(Perimeter.*Perimeter));
                compactness(k,1) = 4*pi*nanmedian(Area./(Perimeter.*Perimeter));
                compactness(k,1) = 4*pi*mean(Area./(Perimeter.*Perimeter));
                % correct for any NaN values due to zero Perimeter -- replace with
                % compactness = 1
                if num_objects == 0
                    compactness(k,1) = 1;
                end
                %----------------------------------------------------------------------
                % compactness heterogeneity -- fractional difference between top and bottom
                % 25% quantiles
                compact_array = 4*pi*Area./(Perimeter.*Perimeter);
                compact_quantiles = quantile(compact_array,[0.25, 0.75]);
                compact_hetero(k,1) = abs(compact_quantiles(1)-compact_quantiles(2))/...
                    (compact_quantiles(1)+compact_quantiles(2));
                % correct for any NaN values due to zero Area -- replace with
                % 0, for no objects
                if isnan(compact_hetero(k,1))
                    compact_hetero(k,1) = 0;
                end
                % ---------------------------------------------------------------------
                % mean or median area
                %     meanarea(k,1) = sum(p_Area.*Area)/num_pixels; %
                %     area-weighted--really same as % pigmented
                meanarea(k,1) = mean(Area)/num_pixels;        % simple mean
                medianarea(k,1) = nanmedian(Area)/num_pixels; % median
                
                if num_objects == 0
                    medianarea(k,1) = 0;
                    meanarea(k,1) = 0;
                end
                if isnan(medianarea(k,1))
                    medianarea(k,1) = 0;
                end
                sqrtmeanarea(k,1) = meanarea(k,1);  % take sqrt to make small features count more toward distances
                %----------------------------------------------------------------------
                % size heterogeneity -- fractional difference between top and bottom
                % 25% quantiles
                Area_quantiles = quantile(Area,[0.25, 0.75]);
                size_hetero(k,1) = abs(Area_quantiles(1)-Area_quantiles(2))/...
                    (Area_quantiles(1)+Area_quantiles(2));
                % correct for any NaN values due to zero Area -- replace with
                % 0, for no objects
                if isnan(size_hetero(k,1))
                    size_hetero(k,1) = 0;
                end
                %----------------------------------------------------------------------
                % eccentricity: range [0,1] -- measure of elongation and extension
                % (i.e., an elongated object could be coiled up into a greater or
                % lesser area)
                % area-weighted average used
                Eccentricity = cat(1, s.Eccentricity);
                eccenvec(k,1) = sum(p_Area.*Eccentricity);
                if num_objects > 0
                    eccenvec(k,1) = nanmedian(Eccentricity);
                else
                    eccenvec(k,1) = 1;
                end
                %----------------------------------------------------------------------
                % area-weighted ratio of object area to its convex hull area
                % this distinguishes between extended stripes and labyrinthin stripes!
                Convex_Area = cat(1, s.ConvexArea);
                labyrinthine(k,1) = mean(Area./Convex_Area);
                % correct for any NaN values due to zero Area -- replace with
                % 0, for no objects
                if isnan(labyrinthine(k,1)) && num_objects == 0
                    labyrinthine(k,1) = 0;
                end
                %----------------------------------------------------------------------
                % how heterogeneous the labyrinthin measure is
                labyr_quantiles = quantile(Area./Convex_Area,[0.25, 0.75]);
                labyr_hetero(k,1) = abs(labyr_quantiles(1)-labyr_quantiles(2))/...
                    (labyr_quantiles(1)+Area_quantiles(2));
                if isnan(labyr_hetero(k,1))
                    labyr_hetero(k,1) = 0;
                end
                %----------------------------------------------------------------------
                ent(k,1) = entropy(img);            % Shannon entropy
                %----------------------------------------------------------------------
                image_canny = edge(BW, 'canny');
                %     figure;imshow(image_canny);
                image_edges = nnz(image_canny);
                num_pixels = numel(img);
                edge_ratio(k,1) = image_edges/num_pixels; % fraction edge pixels
                %----------------------------------------------------------------------
                % first skeletonize the image to turn all objects into lines or points
                if num_objects > 0
                    BW_skel = bwmorph(BW,'skel',Inf);
                    s_skel = regionprops(BW_skel,'Area','Centroid');
                    Area_skel = cat(1, s_skel.Area);
                    obj_width(k,1) = mean(Area./Area_skel)*sqrt(num_objects)/sqrt(num_pixels); % assuming object positions correspond, this should be object width
                else
                    obj_width(k,1) = 0;
                end
                %     centroids = cat(1, s.Centroid);
                %     centroids_skel = cat(1, s_skel.Centroid);
                % find the idx skeleton objects that agree in space with the original
                % using pdist(centroids,centroids_skel);
                % compute width of each object as
                % width = Area./Area_skel(idx);
                %----------------------------------------------------------------------
                % ideally, we should measure how "connected" different points of hte
                % image are--not sure how to do this
                % ---------------------------------------------------------------------
                %Compiles all vectors into single matrix
            end
            PatternMeasures = [sqrtmeanarea(:),fraction_pigmented(:),compactness(:),...
                size_hetero(:),labyrinthine,labyr_hetero(:),obj_num_meas(:),obj_width(:)];
            % varNames = {'sqrt area';'% Pigm';'Compact.';'area het';'labyr';...
            %    'lab het.';'1/sqrt(N)';'width'};
            
            % Combine A, B, C values, and PSM values
            
            ABCandPM = cat(2,ABCValues,PatternMeasures);
            
            % Read code into .csv file with headers.
            % Used code from https://www.mathworks.com/matlabcentral/answers/259937-csvwrite-a-matrix-with-header
            
            cHeader = {'A' 'B' 'C' 'sqrt area' '% Pigm' 'Compact.' 'area het' 'labyr' 'lab het.' '1/sqrt(N)' 'width'}; %dummy header
            commaHeader = [cHeader;repmat({','},1,numel(cHeader))]; %insert commaas
            commaHeader = commaHeader(:)';
            textHeader = cell2mat(commaHeader); %cHeader in text with commas
            %write header to file
            fid = fopen('measfileGM.csv','w');
            fprintf(fid,'%s\n',textHeader);
            fclose(fid);
            %write data to end of file
            dlmwrite('measfileGM.csv',ABCandPM,'-append');
            
            
    end
end


if ismember(2,indx)
    disp('Allen Cat Model Selected')
    disp('Allen Cat Model still in development')
end


if ismember(3,indx)
    disp('Allen Snake Model Selected')
    disp('Allen Snake Model still in development')
end
%end
% end of main function
% *************************************************************************
