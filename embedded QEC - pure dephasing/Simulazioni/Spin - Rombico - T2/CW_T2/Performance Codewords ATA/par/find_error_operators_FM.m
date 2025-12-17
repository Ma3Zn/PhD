clear all

load ../sistema_multispin/Gamma_08apr_FM.txt
Gamma = Gamma_08apr_FM*178.9817570338398;

%Gamma = Gamma*1e3;% /min(min(abs(Gamma)))/4;

for mu=1:10
    for nu=1:10
        gamma(mu,nu) = 2*Gamma(mu,nu)-Gamma(mu,mu)-Gamma(nu,nu);
    end
end

time = 0.01;

code = 1;

alpha = 1/sqrt(2);
beta = sqrt(1-alpha^2);
Ntime = 100;

idle = logspace(log10(0.03),0,Ntime);


%% 4 LEVELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim = 4;

for mu=1:dim
    for nu=1:dim
   %     gamma(mu,nu) = 2*Gamma(mu,nu)-Gamma(mu,mu)-Gamma(nu,nu);
        chi(mu,nu) = exp(gamma(mu,nu)*time);
       % chi(nu,mu) = chi(mu,nu);
    end
    chi(mu,mu) = 1;
end

[V,d] = eig(chi,'vector');
[d, ind] = sort(d, 'descend');
V = V(:, ind);
d(find(d<0)) = 0;

P = zeros(dim,dim,dim);
for j = 1: dim
    P(j,j,j) = 1;
end

%% Error operators
E = zeros(dim,dim,dim);
for j = 1: dim  
    for m = 1: dim
        E(:,:,j) = E(:,:,j) + sqrt(d(j))*V(m,j)*P(:,:,m); 
    end
end

somma = 0;
for j = 1:dim/2
    somma = somma + E(:,:,j)'*E(:,:,j);
end
acc_4 = 1-trace(somma)/dim;

% non so perchè ma serve per far tornare
% trova (-Sz) al posto di Sz)
E(:,:,2) = -E(:,:,2);

%% code words binomiali
S = 3/2;
zero = [4 1];
uno = [2 3];

zL_4 = zeros(4,1);
oL_4 = zeros(4,1);

for j = 1: floor(S+0.5)
    
    zL_4(zero(j)) = sqrt(nchoosek(2*S,2*j-2))/sqrt(2^(2*S-1));
    oL_4(uno(j)) = sqrt(nchoosek(2*S,2*j-1))/sqrt(2^(2*S-1));
        
end


% check K-L conditions
for j = 1: dim/2
    for k = j: dim/2
        
        cond(j,k) = zL_4'*E(:,:,j)'*E(:,:,k)*zL_4 - oL_4'*E(:,:,j)'*E(:,:,k)*oL_4;
        
    end
end

cond_4 = sum(sum(abs(cond)));

%%% error words
ew0 = zeros(dim,dim/2);
ew1 = zeros(dim,dim/2);
for j = 1: dim/2
    ew0(:,j) = E(:,:,j)*zL_4/norm(E(:,:,j)*zL_4);
    ew1(:,j) = E(:,:,j)*oL_4/norm(E(:,:,j)*oL_4);
end

[Q,R] = qr(ew0);
psi00 = Q(:,1);
psi10 = Q(:,2);

[Q,R] = qr(ew1);
psi01 = Q(:,1);
psi11 = Q(:,2);

%%% 
% projector on 0/1 error
P_0 = psi00*psi00' + psi01*psi01';
P_1 = psi10*psi10' + psi11*psi11';

%%%%%%%%%%%%%%%% choose if init in cw or ew0 %%%%%%%%%%%%%%
if code == 1
    
    theta0 = -acos(zL_4'*psi00);
    theta1 = -acos(oL_4'*psi01);

    v0 = psi00-(zL_4'*psi00)*zL_4;
    v0 = v0/norm(v0);
    v1 = psi01-(oL_4'*psi01)*oL_4;
    v1 = v1/norm(v1);

    C0 = (cos(theta0)*(zL_4*zL_4'+v0*v0') + sin(theta0)*(-zL_4*v0'+v0*zL_4') + ...
          cos(theta1)*(oL_4*oL_4'+v1*v1') + sin(theta1)*(-oL_4*v1'+v1*oL_4'));
      
      
    theta0 = -acos(zL_4'*psi10);
    theta1 = -acos(oL_4'*psi11);

    v0 = psi10-(zL_4'*psi10)*zL_4;
    v0 = v0/norm(v0);
    v1 = psi11-(oL_4'*psi11)*oL_4;
    v1 = v1/norm(v1);

    C1 = (cos(theta0)*(zL_4*zL_4'+v0*v0') + sin(theta0)*(-zL_4*v0'+v0*zL_4') + ...
          cos(theta1)*(oL_4*oL_4'+v1*v1') + sin(theta1)*(-oL_4*v1'+v1*oL_4'));
      
    psi_id = alpha*zL_4 + beta*oL_4;
      
else
    
    C1 = psi00*psi10' + psi01*psi11' + psi10*psi00' + psi11*psi01';
    psi_id = alpha*psi00 + beta*psi01;
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rhoin = psi_id*psi_id' ;
Zpsi = alpha*psi00 - beta*psi01;

% loop tempo di attesa
 for it = 1: length(idle)

     delta = idle(it);
     rho = deph(delta,1,rhoin,gamma);
    % t_f = time_enc + delta;
              
     % misura fittizia 0
     rho0 = P_0*rho*P_0/trace(P_0*rho);
     p0(it) = trace(P_0*rho);
     
     if code == 1
         rho0 = C0*rho0*C0';
     end
              
     e_0(it) = 1-psi_id'*rho0*psi_id;
     Le_0(it) = Zpsi'*rho0*Zpsi;
     
     % misura fittizia 1
     rho1 = P_1*rho*P_1/trace(P_1*rho);
     p1(it) = trace(P_1*rho);
     
     % correggo 1
     rho1 = C1*rho1*C1';
   
     e_1(it) = 1-psi_id'*rho1*psi_id;
     Le_1(it) = Zpsi'*rho1*Zpsi;

     % no meas err
     err(it) = e_0(it)*p0(it)+e_1(it)*p1(it);
   %  Lerr(it) = Le_0(it)*p0(it)+Le_1(it)*p1(it);
     
 end

figure(97)

loglog(idle, err, '-b', 'linewidth', 2, 'HandleVisibility','off')
hold on
loglog(idle, (1-exp(-idle))/2, '--k', 'linewidth', 2)

set(gca,'fontsize',20)
xlabel('$t/T_2$', 'interpreter', 'latex')
ylabel('$\mathcal{E}$', 'interpreter', 'latex')
box on
grid on
ylim([1e-10 1])

figure(89)
hold on
loglog(idle, (1-exp(-idle))/2./err, '--b', 'linewidth', 2)
set(gca,'fontsize',20)
xlabel('$t/T_2$', 'interpreter', 'latex')
ylabel('$\mathcal{R}$', 'interpreter', 'latex')
box on
grid on

El(1) = (1-exp(-idle(1)))/2;
El(2) = err(1);

err4 = err;
%% 6 LEVELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except err4 El gamma zL_4 oL_4 time cond_4 acc_4 code alpha beta idle

dim = 6;

for mu=1:dim
    for nu=mu+1:dim
     %   gamma(mu,nu) = 2*Gamma(mu,nu)-Gamma(mu,mu)-Gamma(nu,nu);
        chi(mu,nu) = exp(gamma(mu,nu)*time);
        chi(nu,mu) = chi(mu,nu);
    end
    chi(mu,mu) = 1;
end

[V,d] = eig(chi,'vector');
[d, ind] = sort(d, 'descend');
V = V(:, ind);
d(find(d<0)) = 0;

P = zeros(dim,dim,dim);
for j = 1: dim
    P(j,j,j) = 1;
end

%% Error operators
E = zeros(dim,dim,dim);
for j = 1: dim  
    for m = 1: dim
        E(:,:,j) = E(:,:,j) + sqrt(d(j))*V(m,j)*P(:,:,m); 
    end
end

somma = 0;
for j = 1:dim/2
    somma = somma + E(:,:,j)'*E(:,:,j);
end
acc_6 = 1-trace(somma)/dim;


Err6 = zeros(6,3);
for j = 1: dim/2
    Err6(:,j) = diag(E(:,:,j));
end

format long 
save -ASCII error_operators_6levels_08apr.txt Err6

%% Sistema per code words 
S = 5/2;
zero = [6 2 3];
uno = [4 1 5];

zL_6 = zeros(6,1);
oL_6 = zeros(6,1);

for j = 1: floor(S+0.5)
    
    zL_6(zero(j)) = sqrt(nchoosek(2*S,2*j-2))/sqrt(2^(2*S-1));
    oL_6(uno(j)) = sqrt(nchoosek(2*S,2*j-1))/sqrt(2^(2*S-1));
        
end

% check K-L conditions su TUTTI ERRORI
for j = 1: dim/2
    for k = j: dim/2
        
        cond(j,k) = zL_6'*E(:,:,j)'*E(:,:,k)*zL_6 - oL_6'*E(:,:,j)'*E(:,:,k)*oL_6;
        
    end
end

cond_6 = sum(sum(abs(cond)));


%%% error words
ew0 = zeros(dim,dim/2);
ew1 = zeros(dim,dim/2);
for j = 1: dim/2
    ew0(:,j) = E(:,:,j)*zL_6/norm(E(:,:,j)*zL_6);
    ew1(:,j) = E(:,:,j)*oL_6/norm(E(:,:,j)*oL_6);
end

psi0 = GS(ew0);
psi1 = GS(ew1);

psi01 = psi1(:,1);
psi11 = psi1(:,2);
psi21 = psi1(:,3);

psi00 = psi0(:,1);
psi10 = psi0(:,2);
psi20 = psi0(:,3);

%%% 
% projector on 0/1 error
P_0 = psi0(:,1)*psi0(:,1)' + psi1(:,1)*psi1(:,1)';
P_1 = psi0(:,2)*psi0(:,2)' + psi1(:,2)*psi1(:,2)';
P_2 = psi0(:,3)*psi0(:,3)' + psi1(:,3)*psi1(:,3)';

%%%%%%%%%%%%%%%% choose if init in cw or ew0 %%%%%%%%%%%%%%
if code == 1
    
    theta0 = -acos(zL_6'*psi00);
    theta1 = -acos(oL_6'*psi01);

    v0 = psi00-(zL_6'*psi00)*zL_6;
    v0 = v0/norm(v0);
    v1 = psi01-(oL_6'*psi01)*oL_6;
    v1 = v1/norm(v1);

    C0 = (cos(theta0)*(zL_6*zL_6'+v0*v0') + sin(theta0)*(-zL_6*v0'+v0*zL_6') + ...
          cos(theta1)*(oL_6*oL_6'+v1*v1') + sin(theta1)*(-oL_6*v1'+v1*oL_6'));
      
      
    theta0 = -acos(zL_6'*psi10);
    theta1 = -acos(oL_6'*psi11);

    v0 = psi10-(zL_6'*psi10)*zL_6;
    v0 = v0/norm(v0);
    v1 = psi11-(oL_6'*psi11)*oL_6;
    v1 = v1/norm(v1);

    C1 = (cos(theta0)*(zL_6*zL_6'+v0*v0') + sin(theta0)*(-zL_6*v0'+v0*zL_6') + ...
          cos(theta1)*(oL_6*oL_6'+v1*v1') + sin(theta1)*(-oL_6*v1'+v1*oL_6'));
 
   
    theta0 = -acos(zL_6'*psi20);
    theta1 = -acos(oL_6'*psi21);

    v0 = psi20-(zL_6'*psi20)*zL_6;
    v0 = v0/norm(v0);
    v1 = psi21-(oL_6'*psi21)*oL_6;
    v1 = v1/norm(v1);

    C2 = (cos(theta0)*(zL_6*zL_6'+v0*v0') + sin(theta0)*(-zL_6*v0'+v0*zL_6') + ...
          cos(theta1)*(oL_6*oL_6'+v1*v1') + sin(theta1)*(-oL_6*v1'+v1*oL_6'));
      
    psi_id = alpha*zL_6 + beta*oL_6;
      
else
    
    C1 = psi0(:,1)*psi0(:,2)' + psi1(:,1)*psi1(:,2)' + psi0(:,2)*psi0(:,1)' + psi1(:,2)*psi1(:,1)' + psi0(:,3)*psi0(:,3)' + psi1(:,3)*psi1(:,3)';
    C2 = psi0(:,1)*psi0(:,3)' + psi1(:,1)*psi1(:,3)' + psi0(:,3)*psi0(:,1)' + psi1(:,3)*psi1(:,1)' + psi0(:,2)*psi0(:,2)' + psi1(:,2)*psi1(:,2)';
    
    psi_id = alpha*psi0(:,1) + beta*psi1(:,1);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rhoin = psi_id*psi_id' ;
Zpsi = alpha*psi0(:,1) - beta*psi1(:,1);

% loop tempo di attesa
 for it = 1: length(idle)

     delta = idle(it);
     rho = deph(delta,1,rhoin,gamma);
    % t_f = time_enc + delta;
              
     % misura fittizia 0
     rho0 = P_0*rho*P_0/trace(P_0*rho);
     p0(it) = trace(P_0*rho);
     
     if code == 1
         rho0 = C0*rho0*C0';
     end
              
     e_0(it) = 1-psi_id'*rho0*psi_id;
     Le_0(it) = Zpsi'*rho0*Zpsi;
     
     % misura fittizia 1
     rho1 = P_1*rho*P_1/trace(P_1*rho);
     p1(it) = trace(P_1*rho);
     
     % correggo 1
     rho1 = C1*rho1*C1';
   
     e_1(it) = 1-psi_id'*rho1*psi_id;
     Le_1(it) = Zpsi'*rho1*Zpsi;

    % misura fittizia 2
     rho2 = P_2*rho*P_2/trace(P_2*rho);
     p2(it) = trace(P_2*rho);
     
     % correggo 2
     rho2 = C2*rho2*C2';
   
     e_2(it) = 1-psi_id'*rho2*psi_id;
     Le_2(it) = Zpsi'*rho2*Zpsi;
     
     
     % no meas err
     err(it) = e_0(it)*p0(it)+e_1(it)*p1(it)+e_2(it)*p2(it);
     Lerr(it) = Le_0(it)*p0(it)+Le_1(it)*p1(it)+Le_2(it)*p2(it);
     
 end

figure(97)
hold on
loglog(idle, err, '-r', 'linewidth', 2, 'HandleVisibility','off')

figure(89)
hold on
loglog(idle, (1-exp(-idle))/2./err, '--r', 'linewidth', 2)

El(3) = err(1);

err6 = err;

%% 8 LEVELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except err4 err6  El gamma zL_4 oL_4 zL_6 oL_6 time cond_4 acc_4 cond_6 acc_6 code alpha beta idle Err6

dim = 8;

for mu=1:dim
    for nu=1:dim
    %    gamma(mu,nu) = 2*Gamma(mu,nu)-Gamma(mu,mu)-Gamma(nu,nu);
        chi(mu,nu) = exp(gamma(mu,nu)*time);
        chi(nu,mu) = chi(mu,nu);
    end
    chi(mu,mu) = 1;
end
%chi = (chi'+chi)/2;

[V,d] = eig(chi,'vector');
d(find(abs(d)<1e-12)) = 0;
[d, ind] = sort(d, 'descend');
V = V(:, ind);
d(find(d<0)) = 0;

P = zeros(dim,dim,dim);
for j = 1: dim
    P(j,j,j) = 1;
end

%% Error operators
E = zeros(dim,dim,dim);
for j = 1: dim  
    for m = 1: dim
        E(:,:,j) = E(:,:,j) + sqrt(d(j))*V(m,j)*P(:,:,m); 
    end
end

somma = 0;
for j = 1:dim/2
    somma = somma + E(:,:,j)'*E(:,:,j);
end
acc_8 = 1-trace(somma)/dim;

display('Accuracy rho_reconstruction')
log10([acc_4 acc_6 acc_8])

Err8 = zeros(8,4);
for j = 1: dim/2
    Err8(:,j) = diag(E(:,:,j));
end

format long 
save -ASCII error_operators_8levels_08apr.txt Err8


%% code words binomiali
S = 7/2;
zero = [8 4 1 5];
uno = [6 2 3 7];

zL_8 = zeros(8,1);
oL_8 = zeros(8,1);

for j = 1: floor(S+0.5)
    
    zL_8(zero(j)) = sqrt(nchoosek(2*S,2*j-2))/sqrt(2^(2*S-1));
    oL_8(uno(j)) = sqrt(nchoosek(2*S,2*j-1))/sqrt(2^(2*S-1));
        
end


% check K-L conditions
for j = 1: dim/2
    for k = j: dim/2
        
        cond(j,k) = zL_8'*E(:,:,j)'*E(:,:,k)*zL_8 - oL_8'*E(:,:,j)'*E(:,:,k)*oL_8;
        
    end
end

cond_8 = sum(sum(abs(cond)));


%%% error words
ew0 = zeros(dim,dim/2);
ew1 = zeros(dim,dim/2);
for j = 1: dim/2
    ew0(:,j) = E(:,:,j)*zL_8/norm(E(:,:,j)*zL_8);
    ew1(:,j) = E(:,:,j)*oL_8/norm(E(:,:,j)*oL_8);
end

psi0 = GS(ew0);
psi1 = GS(ew1);

psi01 = psi1(:,1);
psi11 = psi1(:,2);
psi21 = psi1(:,3);
psi31 = psi1(:,4);

psi00 = psi0(:,1);
psi10 = psi0(:,2);
psi20 = psi0(:,3);
psi30 = psi0(:,4);

%%% 
% projector on 0/1 error
P_0 = psi00*psi00' + psi01*psi01';
P_1 = psi10*psi10' + psi11*psi11';
P_2 = psi20*psi20' + psi21*psi21';
P_3 = psi30*psi30' + psi31*psi31';

%%%%%%%%%%%%%%%% choose if init in cw or ew0 %%%%%%%%%%%%%%
if code == 1
    
    theta0 = -acos(zL_8'*psi00);
    theta1 = -acos(oL_8'*psi01);

    v0 = psi00-(zL_8'*psi00)*zL_8;
    v0 = v0/norm(v0);
    v1 = psi01-(oL_8'*psi01)*oL_8;
    v1 = v1/norm(v1);

    C0 = (cos(theta0)*(zL_8*zL_8'+v0*v0') + sin(theta0)*(-zL_8*v0'+v0*zL_8') + ...
          cos(theta1)*(oL_8*oL_8'+v1*v1') + sin(theta1)*(-oL_8*v1'+v1*oL_8'));
      
      
    theta0 = -acos(zL_8'*psi10);
    theta1 = -acos(oL_8'*psi11);

    v0 = psi10-(zL_8'*psi10)*zL_8;
    v0 = v0/norm(v0);
    v1 = psi11-(oL_8'*psi11)*oL_8;
    v1 = v1/norm(v1);

    C1 = (cos(theta0)*(zL_8*zL_8'+v0*v0') + sin(theta0)*(-zL_8*v0'+v0*zL_8') + ...
          cos(theta1)*(oL_8*oL_8'+v1*v1') + sin(theta1)*(-oL_8*v1'+v1*oL_8'));
 
   
    theta0 = -acos(zL_8'*psi20);
    theta1 = -acos(oL_8'*psi21);

    v0 = psi20-(zL_8'*psi20)*zL_8;
    v0 = v0/norm(v0);
    v1 = psi21-(oL_8'*psi21)*oL_8;
    v1 = v1/norm(v1);

    C2 = (cos(theta0)*(zL_8*zL_8'+v0*v0') + sin(theta0)*(-zL_8*v0'+v0*zL_8') + ...
          cos(theta1)*(oL_8*oL_8'+v1*v1') + sin(theta1)*(-oL_8*v1'+v1*oL_8'));
      
    theta0 = -acos(zL_8'*psi30);
    theta1 = -acos(oL_8'*psi31);

    v0 = psi30-(zL_8'*psi30)*zL_8;
    v0 = v0/norm(v0);
    v1 = psi31-(oL_8'*psi31)*oL_8;
    v1 = v1/norm(v1);

    C3 = (cos(theta0)*(zL_8*zL_8'+v0*v0') + sin(theta0)*(-zL_8*v0'+v0*zL_8') + ...
          cos(theta1)*(oL_8*oL_8'+v1*v1') + sin(theta1)*(-oL_8*v1'+v1*oL_8'));
            
    psi_id = alpha*zL_8 + beta*oL_8;
      
else
    
    C1 = psi00*psi10' + psi01*psi11' + psi10*psi00' + psi11*psi01' + psi20*psi20' + psi21*psi21' + psi30*psi30' + psi31*psi31';
    C2 = psi00*psi20' + psi01*psi21' + psi20*psi00' + psi21*psi01' + psi10*psi10' + psi11*psi11' + psi30*psi30' + psi31*psi31';
    C3 = psi00*psi30' + psi01*psi31' + psi30*psi00' + psi31*psi01' + psi10*psi10' + psi11*psi11' + psi20*psi20' + psi21*psi21';
    
    psi_id = alpha*psi00 + beta*psi01;
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rhoin = psi_id*psi_id' ;
Zpsi = alpha*psi00 - beta*psi01;

% loop tempo di attesa
 for it = 1: length(idle)

     delta = idle(it);
     rho = deph(delta,1,rhoin,gamma);
    % t_f = time_enc + delta;
              
     % misura fittizia 0
     rho0 = P_0*rho*P_0/trace(P_0*rho);
     p0(it) = trace(P_0*rho);
     
     if code == 1
         rho0 = C0*rho0*C0';
     end
              
     e_0(it) = 1-psi_id'*rho0*psi_id;
     Le_0(it) = Zpsi'*rho0*Zpsi;
     
     % misura fittizia 1
     rho1 = P_1*rho*P_1/trace(P_1*rho);
     p1(it) = trace(P_1*rho);
     
     % correggo 1
     rho1 = C1*rho1*C1';
   
     e_1(it) = 1-psi_id'*rho1*psi_id;
     Le_1(it) = Zpsi'*rho1*Zpsi;

    % misura fittizia 2
     rho2 = P_2*rho*P_2/trace(P_2*rho);
     p2(it) = trace(P_2*rho);
     
     % correggo 2
     rho2 = C2*rho2*C2';
   
     e_2(it) = 1-psi_id'*rho2*psi_id;
     Le_2(it) = Zpsi'*rho2*Zpsi;
     
     
     % misura fittizia 3
     rho3 = P_3*rho*P_3/trace(P_3*rho);
     p3(it) = trace(P_3*rho);
     
     % correggo 2
     rho3 = C3*rho3*C3';
   
     e_3(it) = 1-psi_id'*rho3*psi_id;
     Le_3(it) = Zpsi'*rho3*Zpsi;
     
     
     % no meas err
     err(it) = e_0(it)*p0(it)+e_1(it)*p1(it)+e_2(it)*p2(it)+e_3(it)*p3(it);
     Lerr(it) = Le_0(it)*p0(it)+Le_1(it)*p1(it)+Le_2(it)*p2(it)+Le_3(it)*p3(it);
     
 end
 
 figure(97)
hold on
loglog(idle, err, '-k', 'linewidth', 2, 'HandleVisibility','off')

figure(89)
hold on
loglog(idle, (1-exp(-idle))/2./err, '--k', 'linewidth', 2)

El(4) = err(1);

err8 = err;

%% 10 LEVELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except err4 err6 err8  El gamma zL_4 oL_4 zL_6 oL_6  zL_8 oL_8 time cond_4 acc_4 cond_6 acc_6 cond_8 acc_8 code alpha beta idle Err6 Err8

dim = 10;

for mu=1:dim
    for nu=1:dim
    %    gamma(mu,nu) = 2*Gamma(mu,nu)-Gamma(mu,mu)-Gamma(nu,nu);
        chi(mu,nu) = exp(gamma(mu,nu)*time);
        chi(nu,mu) = chi(mu,nu);
    end
    chi(mu,mu) = 1;
end
%chi = (chi'+chi)/2;

[V,d] = eig(chi,'vector');
d(find(abs(d)<1e-12)) = 0;
[d, ind] = sort(d, 'descend');
V = V(:, ind);
d(find(d<0)) = 0;

P = zeros(dim,dim,dim);
for j = 1: dim
    P(j,j,j) = 1;
end

%% Error operators
E = zeros(dim,dim,dim);
for j = 1: dim  
    for m = 1: dim
        E(:,:,j) = E(:,:,j) + sqrt(d(j))*V(m,j)*P(:,:,m); 
    end
end

somma = 0;
for j = 1:dim/2
    somma = somma + E(:,:,j)'*E(:,:,j);
end
acc_10 = 1-trace(somma)/dim;

display('Accuracy rho_reconstruction')
log10([acc_4 acc_6 acc_8 acc_10])

Err10 = zeros(10,5);
for j = 1: dim/2
    Err10(:,j) = diag(E(:,:,j));
end
format long 
save -ASCII error_operators_10levels_08apr.txt Err10


%% code words binomiali
S = 9/2;
zero = [10 6 2 3 7];
uno = [8 4 1 5 9];

zL_10 = zeros(10,1);
oL_10 = zeros(10,1);

for j = 1: floor(S+0.5)
    
    zL_10(zero(j)) = sqrt(nchoosek(2*S,2*j-2))/sqrt(2^(2*S-1));
    oL_10(uno(j)) = sqrt(nchoosek(2*S,2*j-1))/sqrt(2^(2*S-1));
        
end


% check K-L conditions
for j = 1: dim/2
    for k = j: dim/2
        
        cond(j,k) = zL_10'*E(:,:,j)'*E(:,:,k)*zL_10 - oL_10'*E(:,:,j)'*E(:,:,k)*oL_10;
        
    end
end

cond_10 = sum(sum(abs(cond)));



%%% error words
ew0 = zeros(dim,dim/2);
ew1 = zeros(dim,dim/2);
for j = 1: dim/2
    ew0(:,j) = E(:,:,j)*zL_10/norm(E(:,:,j)*zL_10);
    ew1(:,j) = E(:,:,j)*oL_10/norm(E(:,:,j)*oL_10);
end

psi0 = GS(ew0);
psi1 = GS(ew1);

psi01 = psi1(:,1);
psi11 = psi1(:,2);
psi21 = psi1(:,3);
psi31 = psi1(:,4);
psi41 = psi1(:,5);

psi00 = psi0(:,1);
psi10 = psi0(:,2);
psi20 = psi0(:,3);
psi30 = psi0(:,4);
psi40 = psi0(:,5);

%%% 
% projector on 0/1 error
P_0 = psi00*psi00' + psi01*psi01';
P_1 = psi10*psi10' + psi11*psi11';
P_2 = psi20*psi20' + psi21*psi21';
P_3 = psi30*psi30' + psi31*psi31';
P_4 = psi40*psi40' + psi41*psi41';

%%%%%%%%%%%%%%%% choose if init in cw or ew0 %%%%%%%%%%%%%%
if code == 1
    
    theta0 = -acos(zL_10'*psi00);
    theta1 = -acos(oL_10'*psi01);

    v0 = psi00-(zL_10'*psi00)*zL_10;
    v0 = v0/norm(v0);
    v1 = psi01-(oL_10'*psi01)*oL_10;
    v1 = v1/norm(v1);

    C0 = (cos(theta0)*(zL_10*zL_10'+v0*v0') + sin(theta0)*(-zL_10*v0'+v0*zL_10') + ...
          cos(theta1)*(oL_10*oL_10'+v1*v1') + sin(theta1)*(-oL_10*v1'+v1*oL_10'));
      
      
    theta0 = -acos(zL_10'*psi10);
    theta1 = -acos(oL_10'*psi11);

    v0 = psi10-(zL_10'*psi10)*zL_10;
    v0 = v0/norm(v0);
    v1 = psi11-(oL_10'*psi11)*oL_10;
    v1 = v1/norm(v1);

    C1 = (cos(theta0)*(zL_10*zL_10'+v0*v0') + sin(theta0)*(-zL_10*v0'+v0*zL_10') + ...
          cos(theta1)*(oL_10*oL_10'+v1*v1') + sin(theta1)*(-oL_10*v1'+v1*oL_10'));
 
   
    theta0 = -acos(zL_10'*psi20);
    theta1 = -acos(oL_10'*psi21);

    v0 = psi20-(zL_10'*psi20)*zL_10;
    v0 = v0/norm(v0);
    v1 = psi21-(oL_10'*psi21)*oL_10;
    v1 = v1/norm(v1);

    C2 = (cos(theta0)*(zL_10*zL_10'+v0*v0') + sin(theta0)*(-zL_10*v0'+v0*zL_10') + ...
          cos(theta1)*(oL_10*oL_10'+v1*v1') + sin(theta1)*(-oL_10*v1'+v1*oL_10'));
      
    theta0 = -acos(zL_10'*psi30);
    theta1 = -acos(oL_10'*psi31);

    v0 = psi30-(zL_10'*psi30)*zL_10;
    v0 = v0/norm(v0);
    v1 = psi31-(oL_10'*psi31)*oL_10;
    v1 = v1/norm(v1);

    C3 = (cos(theta0)*(zL_10*zL_10'+v0*v0') + sin(theta0)*(-zL_10*v0'+v0*zL_10') + ...
          cos(theta1)*(oL_10*oL_10'+v1*v1') + sin(theta1)*(-oL_10*v1'+v1*oL_10'));
            
    theta0 = -acos(zL_10'*psi40);
    theta1 = -acos(oL_10'*psi41);

    v0 = psi40-(zL_10'*psi40)*zL_10;
    v0 = v0/norm(v0);
    v1 = psi41-(oL_10'*psi41)*oL_10;
    v1 = v1/norm(v1);

    C4 = (cos(theta0)*(zL_10*zL_10'+v0*v0') + sin(theta0)*(-zL_10*v0'+v0*zL_10') + ...
          cos(theta1)*(oL_10*oL_10'+v1*v1') + sin(theta1)*(-oL_10*v1'+v1*oL_10'));
              
    psi_id = alpha*zL_10 + beta*oL_10;
      
else
    
    C1 = psi00*psi10' + psi01*psi11' + psi10*psi00' + psi11*psi01' + psi20*psi20' + psi21*psi21' + psi30*psi30' + psi31*psi31' + psi40*psi40' + psi41*psi41';
    C2 = psi00*psi20' + psi01*psi21' + psi20*psi00' + psi21*psi01' + psi10*psi10' + psi11*psi11' + psi30*psi30' + psi31*psi31' + psi40*psi40' + psi41*psi41';
    C3 = psi00*psi30' + psi01*psi31' + psi30*psi00' + psi31*psi01' + psi10*psi10' + psi11*psi11' + psi20*psi20' + psi21*psi21' + psi40*psi40' + psi41*psi41';
    C4 = psi00*psi40' + psi01*psi41' + psi40*psi00' + psi41*psi01' + psi10*psi10' + psi11*psi11' + psi20*psi20' + psi21*psi21' + psi30*psi30' + psi31*psi31';
    
    psi_id = alpha*psi00 + beta*psi01;
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rhoin = psi_id*psi_id' ;
Zpsi = alpha*psi00 - beta*psi01;

% loop tempo di attesa
 for it = 1: length(idle)

     delta = idle(it);
     rho = deph(delta,1,rhoin,gamma);
    % t_f = time_enc + delta;
              
     % misura fittizia 0
     rho0 = P_0*rho*P_0/trace(P_0*rho);
     p0(it) = trace(P_0*rho);
     
     if code == 1
         rho0 = C0*rho0*C0';
     end
              
     e_0(it) = 1-psi_id'*rho0*psi_id;
     Le_0(it) = Zpsi'*rho0*Zpsi;
     
     % misura fittizia 1
     rho1 = P_1*rho*P_1/trace(P_1*rho);
     p1(it) = trace(P_1*rho);
     
     % correggo 1
     rho1 = C1*rho1*C1';
   
     e_1(it) = 1-psi_id'*rho1*psi_id;
     Le_1(it) = Zpsi'*rho1*Zpsi;

    % misura fittizia 2
     rho2 = P_2*rho*P_2/trace(P_2*rho);
     p2(it) = trace(P_2*rho);
     
     % correggo 2
     rho2 = C2*rho2*C2';
   
     e_2(it) = 1-psi_id'*rho2*psi_id;
     Le_2(it) = Zpsi'*rho2*Zpsi;
     
     
     % misura fittizia 3
     rho3 = P_3*rho*P_3/trace(P_3*rho);
     p3(it) = trace(P_3*rho);
     
     % correggo 3
     rho3 = C3*rho3*C3';
   
     e_3(it) = 1-psi_id'*rho3*psi_id;
     Le_3(it) = Zpsi'*rho3*Zpsi;
     
     % misura fittizia 4
     rho4 = P_4*rho*P_4/trace(P_4*rho);
     p4(it) = trace(P_4*rho);
     
     % correggo 4
     rho4 = C4*rho4*C4';
   
     e_4(it) = 1-psi_id'*rho4*psi_id;
     Le_4(it) = Zpsi'*rho4*Zpsi;     
     
     
     % no meas err
     err(it) = e_0(it)*p0(it)+e_1(it)*p1(it)+e_2(it)*p2(it)+e_3(it)*p3(it)+e_4(it)*p4(it);
     Lerr(it) = Le_0(it)*p0(it)+Le_1(it)*p1(it)+Le_2(it)*p2(it)+Le_3(it)*p3(it)+Le_4(it)*p4(it);
     
 end
 
figure(97)
hold on
loglog(idle, err, '-g', 'linewidth', 2, 'HandleVisibility','off')

figure(89)
hold on
loglog(idle, (1-exp(-idle))/2./err, '--g', 'linewidth', 2)

display('Accuracy of Knill-Laflamme conditions')
log10([cond_4 cond_6 cond_8 cond_10])

El(5) = err(1);

err10 = err;
%%
figure(123)
hold on
semilogy([2 4 6 8 10], El, '-or', 'linewidth', 2, 'markersize', 10)

set(gca,'fontsize',24)
xlabel('levels', 'fontsize', 28)
ylabel('$\mathcal{E}$', 'interpreter', 'latex', 'fontsize', 32)
box on

xlim([1.5 12.5])
ylim([1e-10 0.1])
yticks([1e-10 1e-7 1e-4 1e-1])