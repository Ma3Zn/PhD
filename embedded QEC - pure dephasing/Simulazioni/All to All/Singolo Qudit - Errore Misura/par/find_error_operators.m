clear all

load ./Gamma/Gamma_08apr.txt
Gamma = Gamma_08apr*178.9817570338398;  
% fattore di scala consistente col T2 del caso FM

for mu=1:16
    for nu=1:16
        gamma(mu,nu) = 2*Gamma(mu,nu)-Gamma(mu,mu)-Gamma(nu,nu);
    end
end

time = 0.69;

code = 0;

alpha = 1/sqrt(2);
beta = 1/sqrt(2);
Ntime = 100;

idle = logspace(-3,0,Ntime);

% for j = 1: 4
%     for k = 1:4
%         
%         Gamma(2*j,2*k) = 0.25;
%         Gamma(2*j-1,2*k-1) = 0.25;
%         Gamma(2*j-1,2*k) = 0.25;
%         Gamma(2*j,2*k-1) = 0.25;
%         
%     end
% end
% time = 0.001;
% 
% gamma = randn(8)-2*ones(8);
% gamma = gamma-diag(diag(gamma));

%% 4 LEVELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim = 4;

for mu=1:dim
    for nu=mu+1:dim
   %     gamma(mu,nu) = 2*Gamma(mu,nu)-Gamma(mu,mu)-Gamma(nu,nu);
        chi(mu,nu) = exp(gamma(mu,nu)*time);
        chi(nu,mu) = chi(mu,nu);
    end
    chi(mu,mu) = 1;
end


[V,d] = eig(chi,'vector');
[d, ind] = sort(d, 'descend');
V = V(:, ind);
d(find(d<0)) = 0;
% d = d*fact;

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


E(:,:,2) = -E(:,:,2);

somma = 0;
for j = 1:dim/2
    somma = somma + E(:,:,j)'*E(:,:,j);
end
acc_4 = 1-trace(somma)/dim;


%% Sistema per code words 

zero = [1 4];
uno = [2 3];
cw = [zero uno];
z = [0 0 1 1];

for j = 1: length(cw)
    m=1;
    for k = 1: dim/2
        for l = k: dim/2
        
            A(m,j) = sqrt(d(k)*d(l))*V(cw(j),k)*V(cw(j),l)*(-1)^(z(j));
        
            m = m+1;
        end
    end
end


A(m,:) = [1 1 0 0];
A(m+1,:) = [0 0 1 1];

B = A([1 2 end-1 end],:);
b = [0 0 1 1]';


x = sqrt(pinv(B)*b)

zL_4 = zeros(dim,1);
oL_4 = zeros(dim,1);

for k = 1: dim/2
    
    zL_4(zero(k)) = x(k);
    oL_4(uno(k)) = x(k+dim/2);
    
end
% zL_4 = [1 0 0 1]'/sqrt(2);
% oL_4 = [0 1 1 0]'/sqrt(2);

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

cond_4ew = 0;
cond_4 = 0;
% check K-L conditions
for j = 1: dim/2
    for k = j: dim/2
        
        cond_4ew = cond_4ew + ...
            abs( psi00(:,1)'*E(:,:,j)'*E(:,:,k)*psi00(:,1) - psi01(:,1)'*E(:,:,j)'*E(:,:,k)*psi01(:,1) );
        cond_4 = cond_4 + ...
            abs( zL_4'*E(:,:,j)'*E(:,:,k)*zL_4 - oL_4'*E(:,:,j)'*E(:,:,k)*oL_4 );
        
    end
end


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

     err_pop(it) = sum(abs(diag(rho0*p0(it)+rho1*p1(it))-diag(psi_id*psi_id')));
     % no meas err
     err(it) = e_0(it)*p0(it)+e_1(it)*p1(it);
     Lerr(it) = Le_0(it)*p0(it)+Le_1(it)*p1(it);
     
 end

figure(78)
clf
loglog(idle, err, '-b', 'linewidth', 2)

codewords = [zL_4 oL_4];
errorwords0 = [psi00 psi10];
errorwords1 = [psi01 psi11];

% Scriviamo i dati su file in formato binario
fileID = fopen('./data/errw_0_4_lv.bin','w');
fwrite(fileID,errorwords0,'double');
fclose(fileID);

fileID = fopen('./data/errw_1_4_lv.bin','w');
fwrite(fileID,errorwords1,'double');
fclose(fileID);

save -ASCII codewords_4levs_08apr.txt codewords
% save -ASCII errorwords0_4levs_08apr.txt errorwords0
% save -ASCII errorwords1_4levs_08apr.txt errorwords1


%% 6 LEVELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except gamma zL_4 oL_4 time cond_4 cond_4ew acc_4 code alpha beta idle

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

zero = [1 3 6];
uno = [2 4 5];
cw = [zero uno];
z = [0 0 0 1 1 1];

for j = 1: length(cw)
    m=1;
    for k = 1: dim/2
        for l = k: dim/2
        
            A(m,j) = sqrt(d(k)*d(l))*V(cw(j),k)*V(cw(j),l)*(-1)^(z(j));
        
            m = m+1;
        end
    end
end

A(m,:) = [1 1 1 0 0 0];
A(m+1,:) = [0 0 0 1 1 1];


B = A([1:4 7 8],:);
b = [0 0 0 0 1 1]';

%x = sqrt(B\b);
x = sqrt(pinv(B)*b)


zL_6 = zeros(dim,1);
oL_6 = zeros(dim,1);

for k = 1: dim/2
    
    zL_6(zero(k)) = x(k);
    oL_6(uno(k)) = x(k+dim/2);
    
end


% zL_6 = [1 0 0 0.7 0 1.3].'; 
% oL_6 = [0 1 0.7 0 1.3 0].'; 
% 
% zL_6 = zL_6/norm(zL_6);
% oL_6 = oL_6/norm(oL_6);

%%% error words
ew0 = zeros(dim,dim/2);
ew1 = zeros(dim,dim/2);
for j = 1: dim/2
    ew0(:,j) = E(:,:,j)*zL_6/norm(E(:,:,j)*zL_6);
    ew1(:,j) = E(:,:,j)*oL_6/norm(E(:,:,j)*oL_6);
end

% app = ew0(:,2);
% ew0(:,2) = ew0(:,3);
% ew0(:,3) = app;
% 
% app = ew1(:,2);
% ew1(:,2) = ew1(:,3);
% ew1(:,3) = app;

psi0 = GS(ew0);
psi1 = GS(ew1);

psi01 = psi1(:,1);
psi11 = psi1(:,2);
psi21 = psi1(:,3);

psi00 = psi0(:,1);
psi10 = psi0(:,2);
psi20 = psi0(:,3);

% app = psi10;
% psi10 = psi20;
% psi20 = app;
% 
% app = psi11;
% psi11 = psi21;
% psi21 = app;

% check K-L conditions su TUTTI ERRORI
cond_6ew = 0;
cond_6 = 0;
for j = 1: dim/2
    for k = j: dim/2
        
        cond_6ew = cond_6ew + ...
            abs( psi00(:,1)'*E(:,:,j)'*E(:,:,k)*psi00(:,1) - psi01(:,1)'*E(:,:,j)'*E(:,:,k)*psi01(:,1) );
        cond_6 = cond_6 + ...
            abs( zL_6'*E(:,:,j)'*E(:,:,k)*zL_6 - oL_6'*E(:,:,j)'*E(:,:,k)*oL_6 );
        
    end
end


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

figure(78)
hold on
loglog(idle, err, '-r', 'linewidth', 2)

codewords = [zL_6 oL_6];
errorwords0 = [psi00 psi10 psi20];
errorwords1 = [psi01 psi11 psi21];

fileID = fopen('./data/errw_0_6_lv.bin','w');
fwrite(fileID,errorwords0,'double');
fclose(fileID);

fileID = fopen('./data/errw_1_6_lv.bin','w');
fwrite(fileID,errorwords1,'double');
fclose(fileID);

save -ASCII codewords_6levs_08apr.txt codewords
% save -ASCII errorwords0_6levs_08apr.txt errorwords0
% save -ASCII errorwords1_6levs_08apr.txt errorwords1


%% 8 LEVELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except gamma zL_4 oL_4 zL_6 oL_6 time cond_4  cond_4ew  cond_6ew acc_4 cond_6 acc_6 code alpha beta idle Err6

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


%% Sistema per code words 

zero = [1 4 6 7];
uno = [2 3 5 8];
cw = [zero uno];
z = [0 0 0 0 1 1 1 1];

for j = 1: length(cw)
    m=1;
    for k = 1: dim/2
        for l = k: dim/2
        
            A(m,j) = sqrt(d(k)*d(l))*V(cw(j),k)*V(cw(j),l)*(-1)^(z(j));
        
            m = m+1;
        end
    end
end

A(m,:) = [1 1 1 1 0 0 0 0];
A(m+1,:) = [0 0 0 0 1 1 1 1];

B = A([1:6 end-1 end],:);
b = [0 0 0 0 0 0 1 1]';

%x = sqrt(B\b);
x = sqrt(pinv(B)*b)

zL_8 = zeros(dim,1);
oL_8 = zeros(dim,1);

for k = 1: dim/2
    
    zL_8(zero(k)) = x(k);
    oL_8(uno(k)) = x(k+dim/2);
    
end


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


% check K-L conditions su TUTTI ERRORI
cond_8ew = 0;
cond_8 = 0;
for j = 1: dim/2
    for k = j: dim/2
        
        cond_8ew = cond_8ew + ...
            abs( psi00(:,1)'*E(:,:,j)'*E(:,:,k)*psi00(:,1) - psi01(:,1)'*E(:,:,j)'*E(:,:,k)*psi01(:,1) );
        cond_8 = cond_8 + ...
            abs( zL_8'*E(:,:,j)'*E(:,:,k)*zL_8 - oL_8'*E(:,:,j)'*E(:,:,k)*oL_8 );
        
    end
end


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
 
figure(78)
hold on
loglog(idle, err, '-k', 'linewidth', 2)

codewords = [zL_8 oL_8];
errorwords0 = [psi00 psi10 psi20 psi30];
errorwords1 = [psi01 psi11 psi21 psi31];

fileID = fopen('./data/errw_0_8_lv.bin','w');
fwrite(fileID,errorwords0,'double');
fclose(fileID);

fileID = fopen('./data/errw_1_8_lv.bin','w');
fwrite(fileID,errorwords1,'double');
fclose(fileID);

save -ASCII codewords_8levs_08apr.txt codewords
% save -ASCII errorwords0_8levs_08apr.txt errorwords0
% save -ASCII errorwords1_8levs_08apr.txt errorwords1

%% 10 LEVELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except gamma zL_4 oL_4 zL_6 oL_6  zL_8 oL_8 time cond_4 acc_4 cond_6 acc_6 cond_8 acc_8 code alpha beta idle Err6 Err8 cond_4ew cond_6ew cond_8ew

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


%% Sistema per code words 

% zero = randperm(10,5);
% uno = setdiff([1:10],zero);
zero = [1 3 6 8 9];
uno = [2 4 5 7 10];

cw = [zero uno];
z = [0 0 0 0 0 1 1 1 1 1];

for j = 1: length(cw)
    m=1;
    for k = 1: dim/2
        for l = k: dim/2
        
            A(m,j) = sqrt(d(k)*d(l))*V(cw(j),k)*V(cw(j),l)*(-1)^(z(j));
        
            m = m+1;
        end
    end
end

for kk = 1: m-1
    nA(kk) = norm(A(kk,:));
end
[gg,iso] = sort(nA,'descend');


A(m,:) = [1 1 1 1 1 0 0 0 0 0];
A(m+1,:) = [0 0 0 0 0 1 1 1 1 1];

B = (A([1:8 end-1 end],:));
%  B = A([iso(1:8) end-1 end],:);
b = [0 0 0 0 0 0 0 0 1 1]';

%x = sqrt(B\b);
x = sqrt(pinv(B)*b);

zL_10 = zeros(dim,1);
oL_10 = zeros(dim,1);

for k = 1: dim/2
    
    zL_10(zero(k)) = x(k);
    oL_10(uno(k)) = x(k+dim/2);
    
end

zL_10+oL_10

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


% check K-L conditions su TUTTI ERRORI
cond_10ew = 0;
cond_10 = 0;
for j = 1: dim/2
    for k = j: dim/2
        
        cond_10ew = cond_10ew + ...
            abs( psi00(:,1)'*E(:,:,j)'*E(:,:,k)*psi00(:,1) - psi01(:,1)'*E(:,:,j)'*E(:,:,k)*psi01(:,1) );
        cond_10 = cond_10 + ...
            abs( zL_10'*E(:,:,j)'*E(:,:,k)*zL_10 - oL_10'*E(:,:,j)'*E(:,:,k)*oL_10 );
        
    end
end


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
     
     err_pop(it) = sum(abs(diag(rho0*p0(it)+rho1*p1(it)+rho2*p2(it)+rho3*p3(it)+rho4*p4(it))-diag(psi_id*psi_id')));
     % no meas err
     err(it) = e_0(it)*p0(it)+e_1(it)*p1(it)+e_2(it)*p2(it)+e_3(it)*p3(it)+e_4(it)*p4(it);
     Lerr(it) = Le_0(it)*p0(it)+Le_1(it)*p1(it)+Le_2(it)*p2(it)+Le_3(it)*p3(it)+Le_4(it)*p4(it);
     
 end
 
figure(78)
hold on
loglog(idle, err, '-g', 'linewidth', 2)
set(gca,'fontsize',20)
xlabel('$t/T_2$', 'interpreter', 'latex')
ylabel('$\mathcal{E}$', 'interpreter', 'latex')
legend('4 levels', '6 levels', '8 levels', '10 levels','location','southeast')
legend boxoff
box on
grid on

codewords = [zL_10 oL_10];
errorwords0 = [psi00 psi10 psi20 psi30 psi40];
errorwords1 = [psi01 psi11 psi21 psi31 psi41];

fileID = fopen('./data/errw_0_10_lv.bin','w');
fwrite(fileID,errorwords0,'double');
fclose(fileID);

fileID = fopen('./data/errw_1_10_lv.bin','w');
fwrite(fileID,errorwords1,'double');
fclose(fileID);

save -ASCII codewords_10levs_08apr.txt codewords
% save -ASCII errorwords0_10levs_08apr.txt errorwords0
% save -ASCII errorwords1_10levs_08apr.txt errorwords1


display('Accuracy of Knill-Laflamme conditions')
log10([cond_4 cond_6 cond_8 cond_10])
%log10([cond_4ew cond_6ew cond_8ew cond_10ew])

%% 12 LEVELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except gamma zL_4 oL_4 zL_6 oL_6  zL_8 oL_8 zL_10 oL_10 time cond_4 acc_4 cond_6 acc_6 cond_8 acc_8 cond_10 acc_10 code alpha beta idle Err6 Err8 Err10

dim = 12;

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
acc_12 = 1-trace(somma)/dim;

display('Accuracy rho_reconstruction')
log10([acc_4 acc_6 acc_8 acc_10 acc_12])

Err12 = zeros(dim,5);
for j = 1: dim/2
    Err12(:,j) = diag(E(:,:,j));
end
format long 
save -ASCII error_operators_12levels_08apr.txt Err12


%% Sistema per code words 

% zero = randperm(12,6);
% uno = setdiff([1:12],zero);

zero = [1 3 6 8 10 11];
uno = [2 4 5 7 9 12];

cw = [zero uno];
z = [0 0 0 0 0 0 1 1 1 1 1 1];

for j = 1: length(cw)
    m=1;
    for k = 1: dim/2
        for l = k: dim/2
        
            A(m,j) = sqrt(d(k)*d(l))*V(cw(j),k)*V(cw(j),l)*(-1)^(z(j));
        
            m = m+1;
        end
    end
end

A(m,:) = [1 1 1 1 1 1 0 0 0 0 0 0];
A(m+1,:) = [0 0 0 0 0 0 1 1 1 1 1 1];

B = A([1:5 7:9 12 13 end-1 end],:);
b = [0 0 0 0 0 0 0 0 0 0 1 1]';

%x = sqrt(B\b);
x = sqrt(pinv(B)*b);

zL_12 = zeros(dim,1);
oL_12 = zeros(dim,1);

for k = 1: dim/2
    
    zL_12(zero(k)) = x(k);
    oL_12(uno(k)) = x(k+dim/2);
    
end

zL_12


% check K-L conditions
for j = 1: dim/2
    for k = j: dim/2
        
        cond(j,k) = zL_12'*E(:,:,j)'*E(:,:,k)*zL_12 - oL_12'*E(:,:,j)'*E(:,:,k)*oL_12;
        
    end
end

cond_12 = sum(sum(abs(cond)));


%%% error words
ew0 = zeros(dim,dim/2);
ew1 = zeros(dim,dim/2);
for j = 1: dim/2
    ew0(:,j) = E(:,:,j)*zL_12/norm(E(:,:,j)*zL_12);
    ew1(:,j) = E(:,:,j)*oL_12/norm(E(:,:,j)*oL_12);
end

psi0 = GS(ew0);
psi1 = GS(ew1);

psi01 = psi1(:,1);
psi11 = psi1(:,2);
psi21 = psi1(:,3);
psi31 = psi1(:,4);
psi41 = psi1(:,5);
psi51 = psi1(:,6);

psi00 = psi0(:,1);
psi10 = psi0(:,2);
psi20 = psi0(:,3);
psi30 = psi0(:,4);
psi40 = psi0(:,5);
psi50 = psi0(:,6);



% projector on 0/1 error
P_0 = psi00*psi00' + psi01*psi01';
P_1 = psi10*psi10' + psi11*psi11';
P_2 = psi20*psi20' + psi21*psi21';
P_3 = psi30*psi30' + psi31*psi31';
P_4 = psi40*psi40' + psi41*psi41';
P_5 = psi50*psi50' + psi51*psi51';

%%%%%%%%%%%%%%%% choose if init in cw or ew0 %%%%%%%%%%%%%%
if code == 1
    
    theta0 = -acos(zL_12'*psi00);
    theta1 = -acos(oL_12'*psi01);

    v0 = psi00-(zL_12'*psi00)*zL_12;
    v0 = v0/norm(v0);
    v1 = psi01-(oL_12'*psi01)*oL_12;
    v1 = v1/norm(v1);

    C0 = (cos(theta0)*(zL_12*zL_12'+v0*v0') + sin(theta0)*(-zL_12*v0'+v0*zL_12') + ...
          cos(theta1)*(oL_12*oL_12'+v1*v1') + sin(theta1)*(-oL_12*v1'+v1*oL_12'));
      
      
    theta0 = -acos(zL_12'*psi10);
    theta1 = -acos(oL_12'*psi11);

    v0 = psi10-(zL_12'*psi10)*zL_12;
    v0 = v0/norm(v0);
    v1 = psi11-(oL_12'*psi11)*oL_12;
    v1 = v1/norm(v1);

    C1 = (cos(theta0)*(zL_12*zL_12'+v0*v0') + sin(theta0)*(-zL_12*v0'+v0*zL_12') + ...
          cos(theta1)*(oL_12*oL_12'+v1*v1') + sin(theta1)*(-oL_12*v1'+v1*oL_12'));
 
   
    theta0 = -acos(zL_12'*psi20);
    theta1 = -acos(oL_12'*psi21);

    v0 = psi20-(zL_12'*psi20)*zL_12;
    v0 = v0/norm(v0);
    v1 = psi21-(oL_12'*psi21)*oL_12;
    v1 = v1/norm(v1);

    C2 = (cos(theta0)*(zL_12*zL_12'+v0*v0') + sin(theta0)*(-zL_12*v0'+v0*zL_12') + ...
          cos(theta1)*(oL_12*oL_12'+v1*v1') + sin(theta1)*(-oL_12*v1'+v1*oL_12'));
      
    theta0 = -acos(zL_12'*psi30);
    theta1 = -acos(oL_12'*psi31);

    v0 = psi30-(zL_12'*psi30)*zL_12;
    v0 = v0/norm(v0);
    v1 = psi31-(oL_12'*psi31)*oL_12;
    v1 = v1/norm(v1);

    C3 = (cos(theta0)*(zL_12*zL_12'+v0*v0') + sin(theta0)*(-zL_12*v0'+v0*zL_12') + ...
          cos(theta1)*(oL_12*oL_12'+v1*v1') + sin(theta1)*(-oL_12*v1'+v1*oL_12'));
            
    theta0 = -acos(zL_12'*psi40);
    theta1 = -acos(oL_12'*psi41);

    v0 = psi40-(zL_12'*psi40)*zL_12;
    v0 = v0/norm(v0);
    v1 = psi41-(oL_12'*psi41)*oL_12;
    v1 = v1/norm(v1);

    C4 = (cos(theta0)*(zL_12*zL_12'+v0*v0') + sin(theta0)*(-zL_12*v0'+v0*zL_12') + ...
          cos(theta1)*(oL_12*oL_12'+v1*v1') + sin(theta1)*(-oL_12*v1'+v1*oL_12'));
      
    theta0 = -acos(zL_12'*psi50);
    theta1 = -acos(oL_12'*psi51);

    v0 = psi50-(zL_12'*psi50)*zL_12;
    v0 = v0/norm(v0);
    v1 = psi51-(oL_12'*psi51)*oL_12;
    v1 = v1/norm(v1);

    C5 = (cos(theta0)*(zL_12*zL_12'+v0*v0') + sin(theta0)*(-zL_12*v0'+v0*zL_12') + ...
          cos(theta1)*(oL_12*oL_12'+v1*v1') + sin(theta1)*(-oL_12*v1'+v1*oL_12'));
      
              
    psi_id = alpha*zL_12 + beta*oL_12;
      
else
    
    C1 = psi00*psi10' + psi01*psi11' + psi10*psi00' + psi11*psi01' + psi20*psi20' + psi21*psi21' + psi30*psi30' + psi31*psi31' + psi40*psi40' + psi41*psi41' + psi50*psi50' + psi51*psi51';
    C2 = psi00*psi20' + psi01*psi21' + psi20*psi00' + psi21*psi01' + psi10*psi10' + psi11*psi11' + psi30*psi30' + psi31*psi31' + psi40*psi40' + psi41*psi41' + psi50*psi50' + psi51*psi51';
    C3 = psi00*psi30' + psi01*psi31' + psi30*psi00' + psi31*psi01' + psi10*psi10' + psi11*psi11' + psi20*psi20' + psi21*psi21' + psi40*psi40' + psi41*psi41' + psi50*psi50' + psi51*psi51';
    C4 = psi00*psi40' + psi01*psi41' + psi40*psi00' + psi41*psi01' + psi10*psi10' + psi11*psi11' + psi20*psi20' + psi21*psi21' + psi30*psi30' + psi31*psi31' + psi50*psi50' + psi51*psi51';
    C5 = psi00*psi50' + psi01*psi51' + psi50*psi00' + psi51*psi01' + psi10*psi10' + psi11*psi11' + psi20*psi20' + psi21*psi21' + psi30*psi30' + psi31*psi31' + psi40*psi40' + psi41*psi41';
    
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
     
     % misura fittizia 5
     rho5 = P_5*rho*P_5/trace(P_5*rho);
     p5(it) = trace(P_5*rho);
     
     % correggo 4
     rho5 = C5*rho5*C5';
   
     e_5(it) = 1-psi_id'*rho5*psi_id;
     
    % no meas err
     err(it) = e_0(it)*p0(it)+e_1(it)*p1(it)+e_2(it)*p2(it)+e_3(it)*p3(it)+e_4(it)*p4(it)+e_5(it)*p5(it);
     Lerr(it) = Le_0(it)*p0(it)+Le_1(it)*p1(it)+Le_2(it)*p2(it)+Le_3(it)*p3(it)+Le_4(it)*p4(it);
     
 end
 
figure(78)
hold on
loglog(idle, err, '-m', 'linewidth', 2)
set(gca,'fontsize',20)
xlabel('$t/T_2$', 'interpreter', 'latex')
ylabel('$\mathcal{E}$', 'interpreter', 'latex')
legend('4 levs', '6 levs', '8 levs', '10 levs', '12 levs', 'location','southeast')
legend boxoff
box on
grid on

codewords = [zL_12 oL_12];
errorwords0 = [psi00 psi10 psi20 psi30 psi40 psi50];
errorwords1 = [psi01 psi11 psi21 psi31 psi41 psi51];

fileID = fopen('./data/errw_0_12_lv.bin','w');
fwrite(fileID,errorwords0,'double');
fclose(fileID);

fileID = fopen('./data/errw_1_12_lv.bin','w');
fwrite(fileID,errorwords1,'double');
fclose(fileID);

save -ASCII codewords_12levs_08apr.txt codewords
% save -ASCII errorwords0_12levs_08apr.txt errorwords0
% save -ASCII errorwords1_12levs_08apr.txt errorwords1

return
%% 16 LEVELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except gamma zL_4 oL_4 zL_6 oL_6  zL_8 oL_8 zL_10 oL_10  zL_12 oL_12 time cond_4 acc_4 cond_6 acc_6 cond_8 acc_8 acc_10 acc_12 cond_12 cond_10 acc_10 code alpha beta idle Err6 Err8 Err10 Err12

dim = 16;

for mu=1:dim
    for nu=1:dim
    %    gamma(mu,nu) = 2*Gamma(mu,nu)-Gamma(mu,mu)-Gamma(nu,nu);
        chi(mu,nu) = double(vpa(exp(gamma(mu,nu)*time)));
        chi(nu,mu) = chi(mu,nu);
    end
    chi(mu,mu) = 1;
end
%chi = (chi'+chi)/2;


[V,d] = eig(chi,'vector');
d(find(abs(d)<1e-20)) = 0;
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
        E(:,:,j) = E(:,:,j) +  (sqrt(d(j))*V(m,j)*P(:,:,m)); 
    end
end

somma = 0;
for j = 1:dim/2
    somma = somma +  (E(:,:,j)'*E(:,:,j));
end
acc_16 = 1- (trace(somma)/dim);

display('Accuracy rho_reconstruction')
log10([acc_4 acc_6 acc_8 acc_10 acc_12 acc_16])

Err12 = zeros(dim,5);
for j = 1: dim/2
    Err12(:,j) = diag(E(:,:,j));
end
format long 
save -ASCII error_operators_16levels_08apr.txt Err12


%% Sistema per code words 

x = ones(dim,1)*1i;

% for reps = 1: 100000
%     
%     if sum(imag(x) ~= 0)
%         
% zero = sort([randperm(8,4)*2-1  randperm(8,4)*2]);
% %zero = randperm(16,8);
% uno = setdiff([1:16],zero);

%%%%%%%% BUONO, TESTARE LA PERFORMANCE QEC %%%%%%%%%%%%
 zero = [  3     4     5     6     9    12    14    15 ];
 uno = [ 1     2     7     8    10    11    13    16];

%%%%%%%%  PENSO BEST, DA TESTARE CON 32 DIGITS %%%%%%%%%%
% zero = [1     2     7     8    10    11    13    16];
% uno = [3     4     5     6     9    12    14    15];

cw = [zero uno];
z = [0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1];

for j = 1: length(cw)
    m=1;
    for k = 1: dim/2
        for l = k: dim/2
        
            A(m,j) =  (sqrt(d(k)*d(l))*V(cw(j),k)*V(cw(j),l)*(-1)^(z(j)));
        
            m = m+1;
        end
    end
end

A(m,:) = [1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0];
A(m+1,:) = [0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1];

for kk = 1: m-1
    nA(kk) = norm(A(kk,:));
end
[gg,iso] = sort(nA,'descend');

%B =  (A([1:14 end-1 end],:));
B =  (A([iso(1:14) end-1 end],:));
%iso = [1     2     3     4     5     6     7     8     9    10    11    12    13    19];
%B =  (A([iso(1:14) end-1 end],:));
b = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1]';

        x =  (sqrt( (pinv(B))*b))
%        reps
% 
%     end
% 
% end

disp('accuracy matrix inversion')
max(max(abs(B*pinv(B)-eye(16))))
disp('matrix determinant')
det(B)

% B(1:5,:) = A(1:5,:);
% B(15,:) = A(end-1,:);
% B(16,:) = A(end,:);
% E = [];
% G = [];
% c = 0;
% 
% for pp = 6: m-1
%     B(6,:) = A(pp,:);
%     for mm = pp+1: m-1
%         B(7,:) = A(mm,:);
%         for nn = mm+1: m-1
%             B(8,:) = A(nn,:);
%             for qq = nn+1: m-1
%                 B(9,:) = A(qq,:);
%                 for ss = qq+1: m-1
%                     B(10,:) = A(ss,:);
%                     for tt = ss+1: m-1
%                         B(11,:) = A(tt,:);
%                         for vv = tt+1: m-1
%                             B(12,:) = A(vv,:);
%                             for uu = vv+1: m-1
%                                 B(13,:) = A(uu,:);
%                                 for yy = uu+1: m-1
%                                     B(14,:) = A(yy,:);
% 
%                                     c = c+1;
%                                     E(c) = det(B);
%                                     G(c,:) = [pp,mm,nn,qq,ss,tt,vv,uu,yy];
% 
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
%     pp
% end
%                
% 
% [F,c] = max(E);
% [jj,kk,ii,ll,pp,mm,nn,qq,ss,tt,vv,uu,yy] = G(c,:)
% 

zL_16 = zeros(dim,1);
oL_16 = zeros(dim,1);

for k = 1: dim/2
    
    zL_16(zero(k)) = x(k);
    oL_16(uno(k)) = x(k+dim/2);
    
end

% check K-L conditions
for j = 1: dim/2
    for k = j: dim/2
        
        cond(j,k) =  (zL_16'*E(:,:,j)'*E(:,:,k)*zL_16) -  (oL_16'*E(:,:,j)'*E(:,:,k)*oL_16);
        
    end
end

cond_16 =  (sum(sum(abs(cond))));


%%% error words
ew0 = zeros(dim,dim/2);
ew1 = zeros(dim,dim/2);
for j = 1: dim/2
    ew0(:,j) =  (E(:,:,j)*zL_16/norm(E(:,:,j)*zL_16));
    ew1(:,j) =  (E(:,:,j)*oL_16/norm(E(:,:,j)*oL_16));
end

psi0 =  (GS(ew0));
psi1 =  (GS(ew1));

psi01 = psi1(:,1);
psi11 = psi1(:,2);
psi21 = psi1(:,3);
psi31 = psi1(:,4);
psi41 = psi1(:,5);
psi51 = psi1(:,6);
psi61 = psi1(:,7);
psi71 = psi1(:,8);

psi00 = psi0(:,1);
psi10 = psi0(:,2);
psi20 = psi0(:,3);
psi30 = psi0(:,4);
psi40 = psi0(:,5);
psi50 = psi0(:,6);
psi60 = psi0(:,7);
psi70 = psi0(:,8);



% projector on 0/1 error
P_0 =  (psi00*psi00' + psi01*psi01');
P_1 =  (psi10*psi10' + psi11*psi11');
P_2 =  (psi20*psi20' + psi21*psi21');
P_3 =  (psi30*psi30' + psi31*psi31');
P_4 =  (psi40*psi40' + psi41*psi41');
P_5 =  (psi50*psi50' + psi51*psi51');
P_6 =  (psi60*psi60' + psi61*psi61');
P_7 =  (psi70*psi70' + psi71*psi71');

%%%%%%%%%%%%%%%% choose if init in cw or ew0 %%%%%%%%%%%%%%
if code == 1
    
    theta0 = - (acos(zL_16'*psi00));
    theta1 = - (acos(oL_16'*psi01));

    v0 =  (psi00-(zL_16'*psi00)*zL_16);
    v0 =  (v0/norm(v0));
    v1 =  (psi01-(oL_16'*psi01)*oL_16);
    v1 =  (v1/norm(v1));

    C0 =  (cos(theta0)*(zL_16*zL_16'+v0*v0') + sin(theta0)*(-zL_16*v0'+v0*zL_16') + ...
             cos(theta1)*(oL_16*oL_16'+v1*v1') + sin(theta1)*(-oL_16*v1'+v1*oL_16'));
      
    theta0 = - (acos(zL_16'*psi10));
    theta1 = - (acos(oL_16'*psi11));

    v0 =  (psi10-(zL_16'*psi10)*zL_16);
    v0 =  (v0/norm(v0));
    v1 =  (psi11-(oL_16'*psi11)*oL_16);
    v1 =  (v1/norm(v1));

    C1 =  (cos(theta0)*(zL_16*zL_16'+v0*v0') + sin(theta0)*(-zL_16*v0'+v0*zL_16') + ...
             cos(theta1)*(oL_16*oL_16'+v1*v1') + sin(theta1)*(-oL_16*v1'+v1*oL_16'));
         
     
    theta0 = - (acos(zL_16'*psi20));
    theta1 = - (acos(oL_16'*psi21));

    v0 =  (psi20-(zL_16'*psi20)*zL_16);
    v0 =  (v0/norm(v0));
    v1 =  (psi21-(oL_16'*psi21)*oL_16);
    v1 =  (v1/norm(v1));

    C2 =  (cos(theta0)*(zL_16*zL_16'+v0*v0') + sin(theta0)*(-zL_16*v0'+v0*zL_16') + ...
             cos(theta1)*(oL_16*oL_16'+v1*v1') + sin(theta1)*(-oL_16*v1'+v1*oL_16'));
         
    
    theta0 = - (acos(zL_16'*psi30));
    theta1 = - (acos(oL_16'*psi31));

    v0 =  (psi30-(zL_16'*psi30)*zL_16);
    v0 =  (v0/norm(v0));
    v1 =  (psi31-(oL_16'*psi31)*oL_16);
    v1 =  (v1/norm(v1));

    C3 =  (cos(theta0)*(zL_16*zL_16'+v0*v0') + sin(theta0)*(-zL_16*v0'+v0*zL_16') + ...
             cos(theta1)*(oL_16*oL_16'+v1*v1') + sin(theta1)*(-oL_16*v1'+v1*oL_16'));
         
         
         
    theta0 = - (acos(zL_16'*psi40));
    theta1 = - (acos(oL_16'*psi41));

    v0 =  (psi40-(zL_16'*psi40)*zL_16);
    v0 =  (v0/norm(v0));
    v1 =  (psi41-(oL_16'*psi41)*oL_16);
    v1 =  (v1/norm(v1));

    C4 =  (cos(theta0)*(zL_16*zL_16'+v0*v0') + sin(theta0)*(-zL_16*v0'+v0*zL_16') + ...
             cos(theta1)*(oL_16*oL_16'+v1*v1') + sin(theta1)*(-oL_16*v1'+v1*oL_16'));
         
         
    theta0 = - (acos(zL_16'*psi50));
    theta1 = - (acos(oL_16'*psi51));

    v0 =  (psi50-(zL_16'*psi50)*zL_16);
    v0 =  (v0/norm(v0));
    v1 =  (psi51-(oL_16'*psi51)*oL_16);
    v1 =  (v1/norm(v1));

    C5 =  (cos(theta0)*(zL_16*zL_16'+v0*v0') + sin(theta0)*(-zL_16*v0'+v0*zL_16') + ...
             cos(theta1)*(oL_16*oL_16'+v1*v1') + sin(theta1)*(-oL_16*v1'+v1*oL_16'));     
         
      
         
    theta0 = - (acos(zL_16'*psi60));
    theta1 = - (acos(oL_16'*psi61));

    v0 =  (psi60-(zL_16'*psi60)*zL_16);
    v0 =  (v0/norm(v0));
    v1 =  (psi61-(oL_16'*psi61)*oL_16);
    v1 =  (v1/norm(v1));

    C6 =  (cos(theta0)*(zL_16*zL_16'+v0*v0') + sin(theta0)*(-zL_16*v0'+v0*zL_16') + ...
             cos(theta1)*(oL_16*oL_16'+v1*v1') + sin(theta1)*(-oL_16*v1'+v1*oL_16'));
         
         
    
    theta0 = - (acos(zL_16'*psi70));
    theta1 = - (acos(oL_16'*psi71));

    v0 =  (psi70-(zL_16'*psi70)*zL_16);
    v0 =  (v0/norm(v0));
    v1 =  (psi71-(oL_16'*psi71)*oL_16);
    v1 =  (v1/norm(v1));

    C7 =  (cos(theta0)*(zL_16*zL_16'+v0*v0') + sin(theta0)*(-zL_16*v0'+v0*zL_16') + ...
             cos(theta1)*(oL_16*oL_16'+v1*v1') + sin(theta1)*(-oL_16*v1'+v1*oL_16'));
         
         
    psi_id =  (alpha*zL_16 + beta*oL_16);
      
else
    
    C1 =  (psi00*psi10' + psi01*psi11' + psi10*psi00' + psi11*psi01' + psi20*psi20' + psi21*psi21' + psi30*psi30' + psi31*psi31' + psi40*psi40' + psi41*psi41' + psi50*psi50' + psi51*psi51' + psi60*psi60' + psi61*psi61' + psi70*psi70' + psi71*psi71');
    C2 =  (psi00*psi20' + psi01*psi21' + psi20*psi00' + psi21*psi01' + psi10*psi10' + psi11*psi11' + psi30*psi30' + psi31*psi31' + psi40*psi40' + psi41*psi41' + psi50*psi50' + psi51*psi51' + psi60*psi60' + psi61*psi61' + psi70*psi70' + psi71*psi71');
    C3 =  (psi00*psi30' + psi01*psi31' + psi30*psi00' + psi31*psi01' + psi10*psi10' + psi11*psi11' + psi20*psi20' + psi21*psi21' + psi40*psi40' + psi41*psi41' + psi50*psi50' + psi51*psi51' + psi60*psi60' + psi61*psi61' + psi70*psi70' + psi71*psi71');
    C4 =  (psi00*psi40' + psi01*psi41' + psi40*psi00' + psi41*psi01' + psi10*psi10' + psi11*psi11' + psi20*psi20' + psi21*psi21' + psi30*psi30' + psi31*psi31' + psi50*psi50' + psi51*psi51' + psi60*psi60' + psi61*psi61' + psi70*psi70' + psi71*psi71');
    C5 =  (psi00*psi50' + psi01*psi51' + psi50*psi00' + psi51*psi01' + psi10*psi10' + psi11*psi11' + psi20*psi20' + psi21*psi21' + psi30*psi30' + psi31*psi31' + psi40*psi40' + psi41*psi41' + psi60*psi60' + psi61*psi61' + psi70*psi70' + psi71*psi71');
    C6 =  (psi00*psi60' + psi01*psi61' + psi60*psi00' + psi61*psi01' + psi10*psi10' + psi11*psi11' + psi20*psi20' + psi21*psi21' + psi30*psi30' + psi31*psi31' + psi40*psi40' + psi41*psi41' + psi50*psi50' + psi51*psi51' + psi70*psi70' + psi71*psi71');
    C7 =  (psi00*psi70' + psi01*psi71' + psi70*psi00' + psi71*psi01' + psi10*psi10' + psi11*psi11' + psi20*psi20' + psi21*psi21' + psi30*psi30' + psi31*psi31' + psi40*psi40' + psi41*psi41' + psi60*psi60' + psi61*psi61' + psi50*psi50' + psi51*psi51');
    
    psi_id =  (alpha*psi00 + beta*psi01);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rhoin =  (psi_id*psi_id' );

% loop tempo di attesa
 for it = 1: length(idle)

     delta = idle(it);
     rho =  (deph(delta,1,rhoin,gamma));
    % t_f = time_enc + delta;
              
     % misura fittizia 0
     rho0 =  (P_0*rho*P_0);
     p0(it) =  (trace(P_0*rho));
     
     if code == 1
         rho0 =  (C0*rho0*C0');
     end
              
     f_0(it) =  (psi_id'*rho0*psi_id);
     
     % misura fittizia 1
     rho1 =  (P_1*rho*P_1);
     p1(it) =  (trace(P_1*rho));
     
     % correggo 1
     rho1 =  (C1*rho1*C1');
   
     f_1(it) =  (psi_id'*rho1*psi_id);

    % misura fittizia 2
     rho2 =  (P_2*rho*P_2);
     p2(it) =  (trace(P_2*rho));
     
     % correggo 2
     rho2 =  (C2*rho2*C2');
   
     f_2(it) =  (psi_id'*rho2*psi_id);
     
     
     % misura fittizia 3
     rho3 =  (P_3*rho*P_3);
     p3(it) =  (trace(P_3*rho));    
     % correggo 3
     rho3 =  (C3*rho3*C3');
     f_3(it) =  (psi_id'*rho3*psi_id);
     
     % misura fittizia 4
     rho4 =  (P_4*rho*P_4);
     p4(it) =  (trace(P_4*rho));
     % correggo 4
     rho4 =  (C4*rho4*C4');
     f_4(it) =  (psi_id'*rho4*psi_id);  
     
     % misura fittizia 5
     rho5 =  (P_5*rho*P_5);
     p5(it) =  (trace(P_5*rho));
     % correggo 5
     rho5 =  (C5*rho5*C5');
     f_5(it) =  (psi_id'*rho5*psi_id);  
     
     
     % misura fittizia 6
     rho6 =  (P_6*rho*P_6);
     p6(it) =  (trace(P_6*rho));
     % correggo 6
     rho6 =  (C6*rho6*C6');
     f_6(it) =  (psi_id'*rho6*psi_id);  
     
    % no meas err
     err(it) =  (1-f_0(it)-f_1(it)-f_2(it)-f_3(it)-f_4(it)-f_5(it)-f_6(it));
     
 end
 
figure(78)
hold on
loglog(idle, err, '-c', 'linewidth', 2)
set(gca,'fontsize',20)
xlabel('$t/T_2$', 'interpreter', 'latex')
ylabel('$\mathcal{E}$', 'interpreter', 'latex')
legend('4 levels', '6 levs', '8 levs', '10 levs', '12 levs', '16 levs', 'location','southeast')
legend boxoff
box on
grid on

codewords = [zL_12 oL_12];
errorwords0 = [psi00 psi10 psi20 psi30 psi40 psi50];
errorwords1 = [psi01 psi11 psi21 psi31 psi41 psi51];

fileID = fopen('./data/errw_0_16_lv.bin','w');
fwrite(fileID,errorwords0,'double');
fclose(fileID);

fileID = fopen('./data/errw_1_16_lv.bin','w');
fwrite(fileID,errorwords1,'double');
fclose(fileID);

save -ASCII codewords_16levs_08apr.txt codewords
% save -ASCII errorwords0_16levs_08apr.txt errorwords0
% save -ASCII errorwords1_16levs_08apr.txt errorwords1



display('Accuracy of Knill-Laflamme conditions')
log10([cond_4 cond_6 cond_8 cond_10 cond_12 cond_16])

