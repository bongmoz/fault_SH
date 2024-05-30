theta_i= 0;

A_i = 1;

omega=pi*2*10;
h=10;

E=5000*10^6;
nu=0.2;
rho=2100;

E_fault=2000*10^6;
nu_fault=0.2;
rho_fault=2100;

E_liner = 30000*10^6;
nu_liner =0.2;
rho_liner = 2500;

r0 = 3;
r1 = 3.5;
d = 30;



lambda=E*nu/(1+nu)/(1-2*nu);
mu=E/(1+nu)/2;

lambda_liner = E_liner*nu_liner/(1+nu_liner)/(1-2*nu_liner);
mu_liner = E_liner/(1+nu_liner)/2;

lambda_fault = E_fault*nu_fault/(1+nu_fault)/(1-2*nu_fault);
mu_fault = E_fault/(1+nu_fault)/2;



c_ground = sqrt(mu/rho);
k_ground = omega/c_ground;

c_liner = sqrt(mu_liner/rho_liner);
k_liner = omega/c_liner;

c_fault = sqrt(mu_fault/rho_fault);
k_fault = omega/c_fault;

n=10;


if theta_i<pi/2 && theta_i>-pi/2
    
theta_f = asin(k_ground*sin(theta_i)/k_fault);

A1 = exp(1i*k_ground*d*cos(theta_i));
A2 = exp(1i*k_fault*d*cos(theta_f));
A4 = exp(1i*k_ground*(d+h)*cos(theta_i));
A3 = exp(1i*k_fault*(d+h)*cos(theta_f));

lamdazzz = (mu*k_ground*cos(theta_i)-mu_fault*k_fault*cos(theta_f))...
    /(mu*k_ground*cos(theta_i)+mu_fault*k_fault*cos(theta_f));

A_r = A_i*(lamdazzz*A1^2*(A2^2-A3^2)/(A2^2-lamdazzz^2*A3^2));

a_ff = A_i*exp(1i*(pi/2-theta_i)*(-n:n)).'+A_r./exp(1i*(pi/2-theta_i)*(-n:n)).';
end

if theta_i<3*pi/2 && theta_i>pi/2
    
theta_f = pi-asin(k_ground*sin(theta_i)/k_fault);

A1 = exp(1i*k_ground*d*cos(theta_i));
A2 = exp(1i*k_fault*d*cos(theta_f));
A4 = exp(1i*k_ground*(d+h)*cos(theta_i));
A3 = exp(1i*k_fault*(d+h)*cos(theta_f));

lamdazzz = (mu*k_ground*cos(theta_i)-mu_fault*k_fault*cos(theta_f))...
    /(mu*k_ground*cos(theta_i)+mu_fault*k_fault*cos(theta_f));
A_t = (lamdazzz^2-1)*A2*A3*A4/(-A1*A3*A3+lamdazzz^2*A1*A2*A2)*A_i;

a_ff = A_t*exp(1i*(pi/2-theta_i)*(-n:n)).';


end

J_g_r1 = besselj((-n-1:n+1).',k_ground*r1);
H1_g_r1 = besselh((-n-1:n+1).',k_ground*r1);
H1_l_r1 = besselh((-n-1:n+1).',k_liner*r1);
H1_l_r0 = besselh((-n-1:n+1).',k_liner*r0);
H2_l_r1 = besselh((-n-1:n+1).',2,k_liner*r1);
H2_l_r0 = besselh((-n-1:n+1).',2,k_liner*r0);

R1 = (H1_g_r1(1:end-2)-H1_g_r1(3:end))./H1_g_r1(2:end-1);
R2 = (H2_l_r0(1:end-2)-H2_l_r0(3:end))./(H1_l_r0(1:end-2)-H1_l_r0(3:end));

C1 = mu*k_ground*(J_g_r1(1:end-2)-J_g_r1(3:end)-J_g_r1(2:end-1).*R1);
C2 = mu_liner*k_liner*(H1_l_r1(1:end-2)-H1_l_r1(3:end))-mu*k_ground*H1_l_r1(2:end-1).*R1;
C3 = mu_liner*k_liner*(H2_l_r1(1:end-2)-H2_l_r1(3:end))-mu*k_ground*H2_l_r1(2:end-1).*R1;

b2 = (a_ff.*C1)./(-R2.*C2+C3);
b1 = -R2.*b2;
for i=1:360
    theta=pi/180*i;
    sigma_coe(i) = mu_liner/mu/k_ground/r0*abs(sum(1i*(-n:n).*(b1.'.*(H1_l_r0(2:end-1)).'...
        +b2.'.*(H2_l_r0(2:end-1)).').*exp(1i*(-n:n)*theta)));
end

sigmalen=mu_liner/mu/k_ground/r0*(1i*(-n:n).*(b1.'.*(H1_l_r0(2:end-1)).'...
        +b2.'.*(H2_l_r0(2:end-1)).'));
theta = exp(1i*2*pi/360*(1:360)'*(-n:n))';


sig=abs(sigmalen*theta);
plot(sig)
% sig=abs(real(sigmalen*theta*exp(-1i*omega*0.25)));
