%%
Px=[150.000;279.904;334.808;300.000;184.808];
Py=[300.000;184.808;20.0962;-150.000;-279.904];
Qx=[900.000;977.846;934.798;797.690;670.839];
Qy=[500.000;524.479;512.551;445.655;325.301];
alpha=[degtorad(14.9314);degtorad(25.9511);degtorad(39.3782);degtorad(50.1201);degtorad(51.2326)];
f=[alpha(1)-alpha(1);alpha(2)-alpha(1);alpha(3)-alpha(1);alpha(4)-alpha(1);alpha(5)-alpha(1)];
P=[Px(1)+Py(1)*1i;Px(2)+Py(2)*1i;Px(3)+Py(3)*1i;Px(4)+Py(4)*1i;Px(5)+Py(5)*1i];
% P=[sqrt(Px(1)^2+Py(1)^2)*exp(atand(Py(1)/Px(1))*1i);sqrt(Px(2)^2+Py(2)^2)*exp(atand(Py(2)/Px(2))*1i);sqrt(Px(3)^2+Py(3)^2)*exp(atand(Py(3)/Px(3))*1i);sqrt(Px(4)^2+Py(4)^2)*exp(atand(Py(4)/Px(4))*1i);sqrt(Px(5)^2+Py(5)^2)*exp(atand(Py(5)/Px(5))*1i)];
R=[P(1)-P(1);P(2)-P(1);P(3)-P(1);P(4)-P(1);P(5)-P(1)];
% Q=[sqrt(Qx(1)^2+Qy(1)^2)*exp(atand((Qy(1)/Qx(1)))*1i);sqrt(Qx(2)^2+Qy(2)^2)*exp(atand((Qy(2)/Qx(2)))*1i);sqrt(Qx(3)^2+Qy(3)^2)*exp(atand((Qy(3)/Qx(3)))*1i);sqrt(Qx(4)^2+Qy(4)^2)*exp(atand((Qy(4)/Qx(4)))*1i);sqrt(Qx(5)^2+Qy(5)^2)*exp(atand((Qy(5)/Qx(5)))*1i)];
% S=[Q(1)-Q(1);Q(2)-Q(1);Q(3)-Q(1);Q(4)-Q(1);Q(5)-Q(1)];
E2=det([exp(1i*f(3))-1 R(3);exp(1i*f(4))-1 R(4)]);
E3=-det([exp(f(2)*1i)-1 R(2);exp(f(4)*1i)-1 R(4)]);
E4=det([exp(f(2)*1i)-1 R(2);exp(f(3)*1i)-1 R(3)]);
E1=-(E2+E3+E4);
EP2=det([exp(f(3)*1i)-1 R(3);exp(f(5)*1i)-1 R(5)]);
EP3=-det([exp(f(2)*1i)-1 R(2);exp(f(5)*1i)-1 R(5)]);
EP1=-(EP2+EP3+E4);
%%
a1=conj(E1)*conj(EP1);
a2=EP2*E3;
a3=E2*EP3;
a4=a2-a3;
a=a1*a4;
b1=conj(EP1)*conj(E2);
b2=conj(E1)*conj(EP2);
b3=conj(E1)*E3;
b4=conj(EP1)*EP3;
n=(abs(E1)^2)+(abs(E2)^2)+(abs(E3)^2)-(abs(E4)^2);
np=-abs(E4)^2+(abs(EP1)^2+abs(EP2)^2+abs(EP3)^2);
b=b1*a2-b2*a3+b3*np-b4*n;
c1=EP1*E3;
c2=E1*EP3;
c3=conj(E2)*E3;
c4=conj(EP2)*EP3;
c=b2*c1-b1*c2+c3*np-c4*n;
d1=conj(E2)*conj(EP2);
d2=c1-c2;
d=d1*d2;
u=E3*conj(EP3);
f1=conj(E1)*EP2;
fp=f1*u;
h1=EP1*conj(E2);
h=h1*u;
k=fp-conj(h);
K=norm(k);
g1=EP1*conj(E1);
g2=EP2*conj(E2);
g3=g1+g2;
gy=real(u)*imag(g3)+imag(u)*real(g3);
v=gy*4*k*1i;
m=-4*gy^2-2*K^2;
p=a*conj(d);
q=a*conj(c)+b*conj(d)+k^2;
s=a*conj(b)+b*conj(c)+c*conj(d)+v;
t=0.5*(abs(a)^2+abs(b)^2+abs(c)^2+abs(d)^2+m);
A0=real(p)+real(q)+real(s)+t;
A1=-6*imag(p)-4*imag(q)-2*imag(s);
A2=-15*real(p)-5*real(q)+real(s)+3*t;
A3=20*imag(p)-4*imag(s);
A4=15*real(p)-5*real(q)-real(s)+3*t;
A5=-6*imag(p)+4*imag(q)-2*imag(s);
A6=-real(p)+real(q)-real(s)+t;
l0=A1/A6;
l1=A2/A6;
l2=A3/A6;
l3=A4/A6;
l4=A5/A6;
P=[1 l4 l3 l2 l1 l0];
tow=roots(P); 
for i=1:1:5
    if ((1-cos(f(2)))/sin(f(2)))-0.0001<tow(i)&&tow(i)<((1-cos(f(2)))/sin(f(2)))+0.0001
            tow=removerows(tow,i);
            return;
    end
end
%% RIGHT SIDE CALCULATION
st=randi([1,4]);
beta2=angle(1-tow(st)^2+1i*2*tow(st));
tow=removerows(tow,st);
E=E1+E2*exp(1i*beta2);
x=(abs(E4)^2-abs(E3)^2-abs(E)^2)/(2*abs(E3)*abs(E));
beta3=angle(E)+atan2(abs((1-x^2)^0.5),x)-angle(E3);
x=(abs(E3)^2-abs(E4)^2-abs(E)^2)/(2*abs(E4)*abs(E));
beta4=angle(E)+atan2(abs((1-x^2)^0.5),x)-angle(E4);
D1=det([R(2) exp(1i*f(2))-1;R(3) exp(1i*f(3))-1]);
D2=det([exp(1i*beta2)-1 R(2);exp(1i*beta3)-1 R(3)]);
D3=det([exp(1i*beta2)-1 exp(1i*f(2))-1;exp(1i*beta3)-1 exp(1i*f(3))-1]); 
W=D1/D3
Z=D2/D3
%% LEFT SIDE CALCULATION
st=randi([1,3]);
gama2=angle(1-tow(st)^2+1i*2*tow(st));
x=(abs(E4)^2-abs(E3)^2-abs(E)^2)/(2*abs(E3)*abs(E));
gama3=angle(E)+atan2(abs((1-x^2)^0.5),x)-angle(E);
x=(abs(E3)^2-abs(E4)^2-abs(E)^2)/(2*abs(E4)*abs(E));
gama4=angle(E)+atan2(abs((1-x^2)^0.5),x)-angle(E4);
D4=det([R(2) exp(1i*f(2))-1;R(4) exp(1i*f(4))-1]);
D5=det([exp(1i*gama2)-1 R(2);exp(1i*gama4)-1 R(4)]);
D6=det([exp(1i*gama2)-1 exp(1i*f(2))-1;exp(1i*gama4)-1 exp(1i*f(4))-1]);        
U=D4/D6
S=D5/D6



