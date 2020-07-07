%illfunction [A B] = make_dorsal_stripe_snake2

A = zeros(32,416, 4);
B = zeros(32,416, 4);
%add first pattern



prompt = {'min. diffusion speed of B:',...      %answer 1
    'max. diffusion speed of B:',...            %answer 2
    'diffusion speed of A:',...                 %answer 3
    'rate of reaction between A and B:',...     %answer 4
    'noise in Beta:',...                        %answer 5
    'min. concentration of A & B:',...          %answer 6
    'max. concentration of A & B:',...          %answer 7
    };
title = 'Pattern formation specifications';            
dims = [1 35];              % input field specifications
definput = {'0.02','0.07','0.3','0.001','0.1','0','10000'};   % default values
answer = inputdlg(prompt,title,dims,definput);

[Dbmin,statusflag] = str2num(answer{1});
[Dbmax,statusflag] = str2num(answer{2});
[ad.params.Da,statusflag] = str2num(answer{3});
[ad.params.speed,statusflag] = str2num(answer{4});
[ad.params.noise_beta,statusflag] = str2num(answer{5});
[ad.params.MinConst,statusflag] = str2num(answer{6});
[ad.params.MaxConst,statusflag] = str2num(answer{7});


ad.params.Db = [Dbmin  Dbmax]; %0.03 0.04 0.05 0.06
% ad.params.speed = 0.001%[0.00001 0.00025, 0.0005 0.00075 0.001, 0.0025, 0.005 0.0075]; 
 
% ad.params.noise_beta = 0.1%[0 0.1 1 5 ];

cnt=0;

for cnti = 1:3
for cntn = 1:length(ad.params.noise_beta)%Dorsal, head and no prepattern
for cnts = 1:length(ad.params.speed)
for cntb = 1:length(ad.params.Db);
      cnt = cnt+1
ad.params.Da = 0.3;%Speed of diffusion of chemical A
%ad.params.speed = 0.005;%Controls the speed of the reaction between A and B
ad.params.mean_beta = 12;%This value is used in the reaction equations.
%ad.params.noise_beta = 0.1;%Adds some noise to beta, (only used to initialize the equations)
ad.params.MinConst = 0;%Sets the limits of the concentration for stability, both A and B cannot go below this value.
ad.params.MaxConst = 10000;%Similarly this constrains the upper bound of concentration for both chemicals
ad.params.Nx = 480;%The size of the rectangle of interest. 
ad.params.Ny = 96;
ad.params.dt = 3;%The step size of the Euler integration.
ad.iter_completed = 0;
ad.params.init_chem_A = 4*ones(ad.params.Ny, ad.params.Nx); %str2num(get(ad.handles.chemA_txt, 'String'))*
ad.params.init_chem_B = 4*ones(ad.params.Ny, ad.params.Nx); %str2num(get(ad.handles.chemB_txt, 'String'))*
ad.params.chem_A = ad.params.init_chem_A;%Initial condition for chemical A
ad.params.chem_B = ad.params.init_chem_B + (randn(ad.params.Ny, ad.params.Nx)*ad.params.noise_beta(cntn));%Initial condition for chemical A

ad.params.init_beta = ones(ad.params.Ny, ad.params.Nx)*ad.params.mean_beta %+ (randn(ad.params.Ny, ad.params.Nx)*ad.params.noise_beta(cntn));  %add noise to initial conditions
%ad.params.init_beta(48,:)= 16;
a = ad.params.chem_A;
b = ad.params.chem_B;

if cnti ==1
b = b;
elseif cnti == 2
    b(48,:)=16;%dorsal_stripe
elseif cnti == 3
    b(:,33) =16; %head_stripe
end
    
beta = ad.params.init_beta;


while (ad.iter_completed<40000)
    ad.iter_completed = ad.iter_completed + ad.params.dt;
    ab = a.*b;
    da = (ad.params.speed(cnts))*(16 - ab) + (ad.params.Da)*del2(a); %del2 = discrete laplacian
    db = (ad.params.speed(cnts))*(ab - b - beta) + (ad.params.Db(cntb)).*	(b);
    a = a + ad.params.dt*da;
    b = b + ad.params.dt*db;
     
    if ~isempty(ad.params.MinConst)
        a(find(a<ad.params.MinConst)) = ad.params.MinConst;
        b(find(b<ad.params.MinConst)) = ad.params.MinConst;
    end
    if ~isempty(ad.params.MaxConst)
        a(find(a>ad.params.MaxConst)) = ad.params.MaxConst;
        b(find(b>ad.params.MaxConst)) = ad.params.MaxConst;
    end
    
    %Take out the ROI and tile it
aroi = a((ad.params.Ny/3)+1:((ad.params.Ny/3)*2), (ad.params.Nx/15)+1:((ad.params.Nx/15)*14));
broi = b(((ad.params.Ny/3)+1):((ad.params.Ny/3)*2), ((ad.params.Nx/15)+1):((ad.params.Nx/15)*14));
aroir = fliplr(aroi(:,1:32));
broir = fliplr(broi(:,1:32));
a = [aroir aroi aroir;aroir aroi aroir;aroir aroi aroir];
b = [broir broi broir;broir broi broir;broir broi broir];

end

%when reaction complete get the ROI and store it
a = a((ad.params.Ny/3)+1:((ad.params.Ny/3)*2), (ad.params.Nx/15)+1:((ad.params.Nx/15)*14));
b = b(((ad.params.Ny/3)+1):((ad.params.Ny/3)*2), ((ad.params.Nx/15)+1):((ad.params.Nx/15)*14));
B(:,:,cnt) = b;
A(:,:,cnt) = a;
save('snake_dorsal_stripe2B.mat', 'B');
save('snake_dorsal_stripe2A.mat', 'A');

%sz = size(B)
%numLoops

end
end
end

end
