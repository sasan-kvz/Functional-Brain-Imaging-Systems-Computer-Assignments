clearvars;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%
%% Part A : defining parameters realted to spheres
R0 = 7  ; sigma0 = 0.1;
R1 = 8  ; sigma1 = sigma0;
R2 = 8.5; sigma2 = sigma0/80;
R3 = 9  ; sigma3 = sigma0;
%%%%%%%%%%%%%%%%%%%%%%%%
%% Part B : defining a matrix that represents Electrods' positions
Phies = floor((0:0.25:7.75))*(pi/4); %% = [0 0 0 0 1 1 1 ... 7 7 7 7]*(pi/4)
Thetas = (3-mod(0:31,4))*(pi/8); %% = [3 2 1 0 3 2 1 0 ... 3 2 1 0]*(pi/8)
[x,y,z]=sph2cart(Phies,Thetas,R3); %% shpere to cartizian
 x = [0,x];
 y = [0,y];
 z = [R3,z];
 electrods_pos = [x;y;z]; 
%%%%%%%%%%%%%%%%%%%%%%%%
%% Part C : defining a matrix that represents sources' positions
Phies2 = floor((0:0.25:15.75))*(pi/8); %% = [0 0 0 0 1 1 1 ... 15 15 15 15]*(pi/8)
Thetas2 = (3-mod(0:63,4))*(pi/8); %% = [3 2 1 0 3 2 1 0 ... 3 2 1 0]*(pi/8)
[x2,y2,z2]=sph2cart(Phies2,Thetas2,R0); %% shpere to cartizian
 x2 = [0,x2];
 y2 = [0,y2];
 z2 = [R0,z2];
 sources_pos = [x2;y2;z2];
%%%%%%%%%%%%%%%%%%%%%%%%
%% Q 1 : Calculade Gain matrix and plot its 28th coloum 
Gain = ones(33,3*65);
for i = 1:33
    for j = 1:65
        temp = cross(electrods_pos(:,i)/norm(electrods_pos(:,i)), sources_pos(:,j));
        distance = norm(electrods_pos(:,i)-sources_pos(:,j)); 
        distance = distance^3;
        Gain(i,3*(j-1)+1)=temp(1)/distance;
        Gain(i,3*(j-1)+2)=temp(2)/distance;
        Gain(i,3*(j-1)+3)=temp(3)/distance;
    end
end
plot(Gain(:,28));
Gain = Gain*(10^(-7));  %% magnetic permeability in free space / 4*pi = 10^(-7)
%%%%%%%%%%%%%%%%%%%%%%%%
%% Q 2 : Calculade magnetic field strenght matrix for q0 = [1 0 1]^T 
q = [1 0 1 ,zeros(1,195-3)]';
b = Gain*q;
figure
scatter3(x,y,z,40,b,'filled');
%%%%%%%%%%%%%%%%%%%%%%%%
%% Q 3 : 
q_reverse = Gain'*(Gain*Gain')^(-1)*b;
q_rev_amp = zeros(65,1);
for i = 1:65
    q_rev_amp(i)=sqrt((q_reverse(3*(i-1)+1)^2+(q_reverse(3*(i-1)+2))^2+(q_reverse(3*(i-1)+3)^2)));
end
figure
scatter3(x2,y2,z2,40,q_rev_amp,'filled'),title('Qus3 Min norm method');
%%%%%%%%%%%%%%%%%%%%%%%%
%% Q 4 :
b_mean = mean(abs(b));
sigma = b_mean/(10^(2.5));
b_noisy = b + 100*sigma*randn(33,1);
figure;
plot(b_noisy-b),title('b_n_o_i_s_y-b');
q_reverse2 = Gain'*(Gain*Gain')^(-1)*b_noisy;
q_rev_amp2 = zeros(65,1);
for i = 1:65
    q_rev_amp2(i)=sqrt((q_reverse2(3*(i-1)+1)^2+(q_reverse2(3*(i-1)+2))^2+(q_reverse2(3*(i-1)+3)^2)));
end
figure
scatter3(x2,y2,z2,40,q_rev_amp2,'filled'),title('Qus4 Min norm method + noise');
%%%%%%%%%%%%%%%%%%%%%%%%
%% Q 5 :
landa = 0.01;
q_landa = Gain'*(Gain*Gain'+landa*eye(33))^(-1)*b_noisy;
q_rev_amp2 = zeros(65,1);
for i = 1:65
    q_rev_amp2(i)=sqrt((q_reverse2(3*(i-1)+1)^2+(q_reverse2(3*(i-1)+2))^2+(q_reverse2(3*(i-1)+3)^2)));
end
figure
scatter3(x2,y2,z2,40,q_rev_amp2,'fil    led'),title('Qus3 Reg Min norm method + Landa = 0.01');