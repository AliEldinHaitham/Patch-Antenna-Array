function[]=Tapering_Distributions_graphs
close all;
clc;

theta_low=0;
theta_up=180;
disc=181;
M=1800;
k=2*pi;
theta=linspace(0,pi,M+1);
dtheta=pi/M;
% CHOICE OF UNIFORM OR NONUNIFORM
option_c=0;
while ((option_c~=1)&&(option_c~=2))
    disp(strvcat('UNIFORM OR NONUNIFORM ARRAY','OPTION (1):UNIFORM ARRAY','OPTION (2):NONUNIFORM ARRAY'));
    option_c=input('OPTION NUMBER =');
end
if option_c==1 % UNIFORM ARRAY
    % CHOICE OF ARRAY TYPE
    option_d=0;
    while ((option_d~=1)&&(option_d~=2))
        disp(strvcat('ARRAY NAMES','OPTION (1):BROADSIDE ARRAY(MAXIMUM ALONG THETA = 90 DEGREES)',...
            'OPTION (2):ORDINARY END-FIRE ARRAY'));
        option_d=input('OPTION NUMBER =');
    end
    if option_d==1 % BROADSIDE
        Nelem=0;
        while (Nelem<1)
            Nelem=floor(input('NUMBER OF ELEMENTS ='));
        end
        d=input('SPACING d BETWEEN THE ELEMENTS (IN WAVELENGTHS) =');
        beta=0;
        psi=k.*d.*cos(theta)+beta;
        AF=sinc((Nelem.*psi./2)/pi)./sinc((psi./2)/pi);
    elseif option_d==2 % ORDINARY END-FIRE
        Nelem=0;
        while (Nelem<1)
            Nelem=floor(input('NUMBER OF ELEMENTS ='));
        end
        d=input('SPACING d BETWEEN THE ELEMENTS (IN WAVELENGTHS) =');
        thmax=90;
        while((thmax~=0)&&(thmax~=180))
            thmax=input('ANGLE WHERE MAXIMUM OCCURS (THETA = 0 OR 180 DEG.)=');
        end
        if abs(thmax)<eps,beta=-k*d;
        elseif abs(thmax-180)<eps,beta=k*d;
        end
        psi=k*d*cos(theta)+beta;
        AF=sinc((Nelem.*psi./2)/pi)./sinc((psi./2)/pi);
    end
elseif option_c==2 % NONUNIFORM ARRAY
    % CHOICE OF ARRAY TYPE
    option_e=0;
    while ((option_e~=1)&&(option_e~=2)&&(option_e~=3))
        disp(strvcat('ARRAY NAMES',...
            'OPTION (1):BINOMIAL',...
            'OPTION (2):DOLPH-TSCHEBYSCHEFF',...
            'OPTION (3):TAYLOR'));
        option_e=input('OPTION NUMBER =');
    end
    if option_e==1 % BINOMIAL
        Nelem=0;
        while (Nelem<1)
            Nelem=floor(input('NUMBER OF ELEMENTS ='));
        end
        d=input('SPACING d BETWEEN THE ELEMENTS (IN WAVELENGTHS) =');
        beta=0;
        [AF,Ncoef,Coef]=bin(cos(theta),Nelem,d,beta);
    elseif option_e==2 % DOLPH-TSCHEBYSCHEFF
        Nelem=0;
        while (Nelem<1)
            Nelem=floor(input('NUMBER OF ELEMENTS ='));
        end
        d=input('SPACING d BETWEEN THE ELEMENTS (IN WAVELENGTHS) =');
        beta=0;
        RdB=input('SIDE LOBE LEVEL (IN dB) =');
        [AF,Ncoef,Coef,Zo]=tscheby(cos(theta),Nelem,d,RdB,beta);
        R=10^(RdB/20);
        f=1+0.636*(2/R*cosh(sqrt((acosh(R))^2-pi^2)))^2;
        Dir=2*R^2/(1+(R^2-1)*f/(d*Nelem));
        HP_uni=acos(-0.443/(Nelem*d))-acos(0.443/(Nelem*d));
        HPT=HP_uni*f*180/pi;
    elseif option_e==3 % TAYLOR
        Nelem=0;
        while (Nelem < 2)
            Nelem=floor(input('NUMBER OF ELEMENTS ='));
        end
        d=input('SPACING d BETWEEN THE ELEMENTS (IN WAVELENGTHS) =');

        RdB = 0;
        while RdB >= 0
            RdB = -input('SIDE LOBE LEVEL (IN dB, e.g., 20) = ');
        end

        n_bar = 0;
        while n_bar < 1 || mod(n_bar,1) ~= 0
            n_bar = floor(input('TAYLOR n_bar: NO. OF NEARLY EQUAL SIDELOBES (>=1) = '));
        end

        beta = 0; % Broadside by default

        % Generate Taylor weights
        I_taylor = taylorwin(Nelem, n_bar, RdB); % Use built-in function
        I_taylor = I_taylor / max(I_taylor);   % Normalize

        % Compute AF over theta
        AF = zeros(size(theta));
        for n = 1:Nelem
            phase = 2 * pi * d * (n - (Nelem + 1)/2) * cos(theta) + beta;
            AF = AF + I_taylor(n) * exp(1j * phase);
        end
        AF = abs(AF);
        AF = AF / max(AF); % Normalize

        U=(abs(AF)).^2;
        Prad=2*pi*sum(U.*sin(theta).*dtheta);
        D=4*pi*U./Prad;
        DdB=10*log10(D + eps);
        Do=max(D);
        DodB=max(DdB);

        % HPBW calculation
        if Nelem == 1 || max(AF) < 2*min(AF)
            hp = 0;
            thmax = 0;
        else
            [hp,thmax] = hpbw(AF);
        end
        No_maxima=length(thmax);

        Coef = I_taylor;
        Ncoef = Coef / max(Coef);
    end
    Coef=Coef(1:Ncoef);
    Ncoef=Coef(1:Ncoef)/Coef(Ncoef);
end
U=(abs(AF)./max(abs(AF))).^2;
Prad=2*pi*sum(U.*sin(theta).*dtheta);
D=4*pi*U/Prad;
DdB=10.*log10(D+eps);
Do=max(D);
DodB=max(DdB);
if Nelem==1||max(AF)<2*min(AF)
    hp=0;
    thmax=0;
else
    [hp,thmax]=hpbw(U,M);
end
No_maxima=length(thmax);

% % OUTPUT
% disp(strvcat('********************************************************'));
% disp(strvcat('PROGRAM OUTPUT'));
% disp(strvcat('********************************************************'));
% disp(strvcat('INPUT SPECIFICATION'));
% disp(strvcat('--------------------------------------------------------'));
% if option_c==1&&option_d==1
%     disp(strvcat('UNIFORM BROADSIDE ARRAY'));
% end
% if option_c==1&&option_d==2
%     disp(strvcat('UNIFORM ORDINARY END-FIRE ARRAY'));
% end
% if option_c==2&&option_e==1
%     disp(strvcat('NONUNIFORM BINOMIAL ARRAY'));
% end
% if option_c==2&&option_e==2
%     disp(strvcat('NONUNIFORM DOLPH-TSCHEBYSCHEFF ARRAY'));
% end
% if option_c==2&&option_e==3
%     disp(strvcat('NONUNIFORM TAYLOR ARRAY'));
% end
% disp(['NUMBER OF ARRAY ELEMENTS = ',num2str(Nelem)]);
% disp(['SPACING BETWEEN THE ELEMENTS (IN WAVELENGTHS) = ',num2str(d)]);
% if option_c==1&&option_d~=1
%     disp(['MAXIMUM NEEDS TO OCCUR AT = ',num2str(thmax)]);
% end
% if option_c==2&&option_e==2
%     disp(['SIDE LOBE LEVEL (IN dB) = ',num2str(RdB)]);
% end
% if option_c==2&&option_e==3
%     disp(['SIDE LOBE LEVEL (IN dB) = ',num2str(-RdB)]);
%     disp(['TAYLOR n_bar = ',num2str(n_bar)]);
% end
% disp('OUTPUT CHARACTERISTICS OF THE ARRAY');
% disp('--------------------------------------------------------');
% disp(['DIRECTIVITY = ',num2str(DodB),' dB']);
% disp(['DIRECTIVITY = ',num2str(Do),' dimensionless']);
% disp(['NUMBER OF MAXIMA BETWEEN 0 AND 180 DEGREES = ',num2str(No_maxima)]);
% for i=1:No_maxima;
%     disp(['HPBW using 0.5*AF^2 FOR MAXIMUM # =',num2str(i),'    ',num2str(hp(i,:)),' degrees     THMAX = ',num2str(thmax(i)),' degrees']);
% end
% if option_c==2  
%     if option_e==2
%         fprintf('\nBEAM BROADENING FACTOR (BBF), f = %6.4f \n', f);
%         fprintf('DIRECTIVITY using (6-79) = %6.4f dB \n', 10*log10(Dir));
%         fprintf('DIRECTIVITY using (6-79) = %6.4f dimensionless\n', Dir);
%         fprintf('HALF POWER BEAMWIDTH for uniform array, using (6-22a)= %6.4f degrees\n', HP_uni*180/pi);
%         fprintf('HALF POWER BEAMWIDTH for Dolph-Tschebyscheff array computed using\n')
%         fprintf('      HPBW(uniform)*f(BBF) = %6.4f degrees\n', HPT);
%     end
% end
% if option_c==2            
%     disp('TOTAL EXCITATION COEFFICIENTS FOR THE ARRAY DESIGN');
%     disp(Coef);
%     disp('NORMALIZED TOTAL EXCITATION COEFFICIENTS (RELATIVE TO EDGE)');
%     Ncoef=Coef/Coef(round(Nelem/2));
%     disp(Ncoef);
%     disp('NORMALIZED TOTAL EXCITATION COEFFICIENTS (RELATIVE TO CENTER)');
%     Ncoef=Ncoef./Ncoef(1);
%     disp(Ncoef);
% end
% diary off;
AFdB=10.*log10(U);
for i=1:901
    thetarec(i)=theta(i*2-1);
    AFdBrec(i)=AFdB(i*2-1);
end
for i=902:1801
    thetarec(i)=2*pi-theta(3601-i*2+3);
    AFdBrec(i)=AFdB(3601-i*2+3);
end
fidaf=fopen('ArrFac.dat','wt');
fprintf(fidaf,'%7.3f        %9.5f\n',[thetarec.*180/pi; AFdBrec]);
fclose(fidaf);
% PLOT THE GRAPHS
% ARRAY FACTOR
figure;
plot(theta*180/pi,AFdB,'m','linewidth',2);
axis([0 180 max(min(AFdB)-1,-60) 1]);
xlabel(['\theta',' (degrees)']),ylabel('ARRAY FACTOR(dB)')
grid on;

if option_c==1&&option_d==1
    s1=title('LINEAR UNIFORM BROADSIDE','Fontsize',15);
    set(gca,'units','normalized');
    set(s1,'position',[0 1],'horizontalalign','left');
end
if option_c==1&&option_d==2
    s2=title('LINEAR UNIFORM ORDINARY END-FIRE','Fontsize',15);
    set(gca,'units','normalized');
    set(s2,'position',[0 1],'horizontalalign','left');
end
if option_c==2&&option_e==1
    s5=title('LINEAR NONUNIFORM BINOMIAL','Fontsize',15);
    set(gca,'units','normalized');
    set(s5,'position',[0 1],'horizontalalign','left');
end
if option_c==2&&option_e==2
    s6=title('LINEAR NONUNIFORM DOLPH-TSCHEBYSCHEFF','Fontsize',10);
    set(gca,'units','normalized');
    set(s6,'position',[10 1],'horizontalalign','left');
end
if option_c==2&&option_e==3
    s7=title('LINEAR NONUNIFORM TAYLOR','Fontsize',15);
    set(gca,'units','normalized');
    set(s7,'position',[0 1],'horizontalalign','left');
end
% =============================
% ADD MARKERS AT Y = -3 dB
% =============================
hold on;

threshold = -3;

% Remove any variable named diff that may override the function
try
    dummy = diff([1 2]);  % Test if diff is functional
catch
    clear diff;           % Clear variable shadowing
end

% Find crossing indices
dAF = sign(AFdB - threshold);
cross_indices = find(abs(dAF(2:end) - dAF(1:end-1)) > 1e-6);

for k = 1:length(cross_indices)
    idx = cross_indices(k);
    
    % Interpolate between theta(idx) and theta(idx+1)
    x_before = theta(idx) * 180/pi;
    x_after = theta(idx+1) * 180/pi;
    y_before = AFdB(idx);
    y_after = AFdB(idx+1);
    
    slope = (y_after - y_before) / (x_after - x_before + eps);
    if abs(slope) < 1e-6
        continue;
    end
    
    x_cross = x_before + (threshold - y_before)/slope;
    
    % Plot a vertical dashed line
    plot([x_cross x_cross], [min(ylim) threshold], 'k--', 'LineWidth', 1);
    
    % Plot a red circle at the crossing point
    plot(x_cross, threshold, 'ro', 'MarkerFaceColor', 'r');
end

hold off;
%Array Factor linear
figure;
plot(theta*180/pi,abs(AF)/max(abs(AF)),'m','linewidth',2);
axis([0 180 0  1]);
xlabel(['\theta',' (degrees)']),ylabel('NORMALIZED ARRAY FACTOR')
grid on;

if option_c==1&&option_d==1
    s1=title('LINEAR UNIFORM BROADSIDE','Fontsize',15);
    set(gca,'units','normalized');
    set(s1,'position',[0 1],'horizontalalign','left');
end
if option_c==1&&option_d==2
    s2=title('LINEAR UNIFORM ORDINARY END-FIRE','Fontsize',15);
    set(gca,'units','normalized');
    set(s2,'position',[0 1],'horizontalalign','left');
end
if option_c==2&&option_e==1
    s5=title('LINEAR NONUNIFORM BINOMIAL','Fontsize',15);
    set(gca,'units','normalized');
    set(s5,'position',[0 1],'horizontalalign','left');
end
if option_c==2&&option_e==2
    s6=title('LINEAR NONUNIFORM DOLPH-TSCHEBYSCHEFF','Fontsize',15);
    set(gca,'units','normalized');
    set(s6,'position',[0 1],'horizontalalign','left');
end
if option_c==2&&option_e==3
    s7=title('LINEAR NONUNIFORM TAYLOR','Fontsize',15);
    set(gca,'units','normalized');
    set(s7,'position',[0 1],'horizontalalign','left');
end
% DIRECTIVITY
figure;
diff=Do-min(D);
subplot(2,1,1)
plot(theta*180/pi,D,'r','linewidth',2);
xlabel(['\theta',' (degrees)']),ylabel('DIRECTIVITY(dimensionless)')
grid on;
axis([0 180 floor(min(D)-0.1*diff-.1) ceil(Do+0.1*diff+.1)]);
t2=text(1,1,['D_0 = ',num2str(Do),' (dimensionless)']);
set(t2,'units','normalized','position',[1 1.05],'horizontalalign','right');
if option_c==1&&option_d==1
    s1=title('LINEAR UNIFORM BROADSIDE','Fontsize',15);
    set(gca,'units','normalized');
    set(s1,'units','normalized','position',[0 1],'horizontalalign','left');
end
if option_c==1&&option_d==2
    s2=title('LINEAR UNIFORM ORDINARY END-FIRE','Fontsize',15);
    set(gca,'units','normalized');
    set(s2,'units','normalized','position',[0 1],'horizontalalign','left');
end
if option_c==2&&option_e==1
    s5=title('LINEAR NONUNIFORM BINOMIAL','Fontsize',15);
    set(gca,'units','normalized');
    set(s5,'units','normalized','position',[0 1],'horizontalalign','left');
end
if option_c==2&&option_e==2
    s6=title('LINEAR NONUNIFORM DOLPH-TSCHEBYSCHEFF','Fontsize',15);
    set(gca,'units','normalized');
    set(s6,'units','normalized','position',[0 1],'horizontalalign','left');
end
if option_c==2&&option_e==3
    s7=title('LINEAR NONUNIFORM TAYLOR','Fontsize',15);
    set(gca,'units','normalized');
    set(s7,'units','normalized','position',[0 1],'horizontalalign','left');
end
diffdB=DodB-min(DdB);
subplot(2,1,2)
plot(theta*180/pi,DdB,'b','linewidth',2);
t3=text(1,1,['D_0 = ',num2str(DodB),' (dB)']);
set(t3,'units','normalized','position',[1 1.05],'horizontalalign','right');
xlabel(['\theta',' (degrees)']),ylabel('DIRECTIVITY(dB)')
grid on;
axis([0 180 max(-50,10*floor(min(DdB)/10)) 10*ceil(DodB/10)]);
% AMPLITUDE AND PHASE OF EXCITATION COEFFICIENTS
if option_c~=2
    for i=1:Nelem
        Ncoef(i)=1;
    end
end
figure;
x=(1-Nelem)*d/2:d:(Nelem-1)*d/2;
y=(1-Nelem)/2:1:(Nelem-1)/2;
[AX,H1,H2]=plotyy(x,Ncoef(ceil(abs(y)+0.1)),x,beta.*y.*180./pi);
set(get(AX(1),'Ylabel'),'String','AMPLITUDE','color','r');
set(get(AX(2),'Ylabel'),'String','PHASE (degrees)','color','b');
set(AX(1),'ycolor','r');
set(AX(2),'ycolor','b');
set(H1,'Linestyle','-','color','r','linewidth',2,'marker','s');
set(H2,'Linestyle',':','color','b','linewidth',2,'marker','o');
xlabel(['ARRAY LENGTH',' (\lambda)']);
grid on;
hle = legend('AMPLITUDE','PHASE');
set(hle,'color',[1 1 1], 'Location', 'best');
if option_c==1&&option_d==1
    s1=title('LINEAR UNIFORM BROADSIDE','Fontsize',15);
    set(s1,'units','normalized','position',[0 1],'horizontalalign','left');
end
if option_c==1&&option_d==2
    s2=title('LINEAR UNIFORM ORDINARY END-FIRE','Fontsize',15);
    set(gca,'units','normalized');
    set(s2,'units','normalized','position',[0 1],'horizontalalign','left');
end
if option_c==2&&option_e==1
    s5=title('LINEAR NONUNIFORM BINOMIAL','Fontsize',15);
    set(gca,'units','normalized');
    set(s5,'units','normalized','position',[0 1],'horizontalalign','left');
end
if option_c==2&&option_e==2
    s6=title('LINEAR NONUNIFORM DOLPH-TSCHEBYSCHEFF','Fontsize',15);
    set(gca,'units','normalized');
    set(s6,'units','normalized','position',[0 1],'horizontalalign','left');
end
if option_c==2&&option_e==3
    s7=title('LINEAR NONUNIFORM TAYLOR','Fontsize',15);
    set(gca,'units','normalized');
    set(s7,'units','normalized','position',[0 1],'horizontalalign','left');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Subroutines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HPBWCALC
function[hp,thmax]=hpbw(U,M)
tol=0.001;
imax=0;
j=0;
for i=1:M+1
    if abs(U(i)-1)<tol && floor((j+2)/2)==imax+1
        imax=imax+1;
        thmax(imax)=(i-1)/10;
    end
    if i>1 && abs(U(i)-1)<tol && U(i)>U(i-1) && j~=0
        thmax(imax)=(i-1)/10;
    end
    if i>1
        y(1)=U(i)-0.5;
        y(2)=U(i-1)-0.5;
        x(1)=(i-1)/10;
        x(2)=(i-2)/10;
        sign=y(1)*y(2);
        if sign<0
            j=j+1;
            root(j)=x(2)-y(2)*(x(2)-x(1))/(y(2)-y(1));
            if j>=2 && y(2)>y(1)
                hp(imax,1)=root(j)-root(j-1);
            elseif j==1 && y(2)>y(1)
                hp(imax,1)=2.*root(j);
            end
        end
    end
end
if j==0
    hp(imax,:)='N/A';
elseif thmax(imax)>root(j)
    hp(imax,1)=2.*(180-root(j));
else
end
% BINOMIAL
function[AF,Ncoef,Coef]=bin(X,Nelem,d,beta)
if 2*floor(Nelem/2)==Nelem,Ncoef=Nelem/2;
else 
    Ncoef=(Nelem+1)/2;
end
Coef=zeros(1,Ncoef);
for i=1:Ncoef
    Coef(i)=1;
    for j=1:Ncoef-i
        Coef(i)=Coef(i).*(Nelem-j)./j;
    end
end
if 2*floor(Nelem/2)~=Nelem,Coef(1)=Coef(1)/2;
end
u=2*pi*d*X+beta;
if 2*floor(Nelem/2)==Nelem
    AF=0;
    for i=1:Ncoef
        AF=AF+Coef(i).*cos((2.*i-1)/2.*u);
    end
else 
    AF=0;
    for i=1:Ncoef
        AF=AF+Coef(i).*cos((i-1).*u);
    end
end
if 2*round(Nelem/2)~=Nelem
    Coef(1)=2*Coef(1);
end
% TSCEBY(THETA,NELEM,D,RDB)
function[AF,Ncoef,Coef,Zo]=tscheby(X,Nelem,d,RdB,beta)
Ro=10^(RdB/20);
P=Nelem-1;
Zo=0.5*((Ro+sqrt(Ro^2-1))^(1/P)+(Ro-sqrt(Ro^2-1))^(1/P));
if 2*floor(Nelem/2)==Nelem
    Ncoef=Nelem/2;
    M=Ncoef;
    Coef=zeros(1,M);
    for i=1:M
        Coef(i)=0;
        for j=i:M
            Coef(i)=Coef(i)+(-1)^(M-j)*Zo^(2*j-1)*fact(j+M-2)*(2*M)/(fact(j-i)*fact(j+i-1)*fact(M-j));
        end
    end
elseif 2*floor((Nelem+1)/2)==Nelem+1
    Ncoef=(Nelem+1)/2;
    M=Ncoef-1;
    Coef=zeros(1,M);
    for i=1:M+1
        Coef(i)=0;
        for j=i:M+1
            if i==1
                EN=2;
            else
                EN=1;
            end
            Coef(i)=Coef(i)+(-1)^(M-j+1)*Zo^(2*(j-1))*fact(j+M-2)*2*M/(EN*fact(j-i)*fact(j+i-2)*fact(M-j+1));
        end
    end
end
u=2*pi*d*X+beta;
if 2*floor(Nelem/2)==Nelem
    AF=0;
    for i=1:Ncoef
        AF=AF+Coef(i)*cos((2*i-1)/2*u);
    end
elseif 2*floor((Nelem+1)/2)==Nelem+1
    AF=0;
    for i=1:Ncoef
        AF=AF+Coef(i)*cos((i-1)*u);
    end
end
if 2*round(Nelem/2)~=Nelem
    Coef(1)=2*Coef(1);
end
% FACT(IARG)
function[f7]=fact(iarg)
f7=1;
for j=1:iarg
    f7=j*f7;
end
function y=sinc(x)
hh=find(x==0);
x(hh)= 1;
y = sin(pi*x)./(pi*x);
y(hh) = 1;