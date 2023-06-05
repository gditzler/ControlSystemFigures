%% chapter 1 figures 
clc; 
clear; 
close all; 

% -------------------------------------------------------------------------
t = linspace(0, 1.5, 200);
Tp = .15;
PO = .5;
yfv= 1.5;
wd = pi/Tp;
zeta = sqrt(log(PO)^2/(pi^2+log(PO)^2));
wn = wd/(sqrt(1-zeta^2));
sigma = wn*zeta;
beta = sqrt(1-zeta^2);
theta = acos(zeta);
alpha = wn^2*yfv;
f = (1 - 1/beta*exp(-sigma*t).*sin(wd*t+theta))*yfv;

h = figure();
hold on; 
box on; 
grid on; 
plot([.15,.15], [0,2.25], 'c', 'LineWidth', 3)
plot([1.25,1.25], [0,1.5], 'c', 'LineWidth', 3)
plot(t, f, 'b', 'LineWidth', 3)
plot(1.25, yfv, 'r-x', 'MarkerSize', 15, 'LineWidth', 3)  
plot(.15, max(f), 'r-x', 'MarkerSize', 15, 'LineWidth', 3)  
xlabel('time (s)','Interpreter', 'latex', 'FontSize', 18)
ylabel('$y(t)$','Interpreter', 'latex', 'FontSize', 18)
saveas(h, 'pdf/system.eps', 'eps2c')
%% chapter 2 figures 
clc; clear
close all;
t = linspace(-3, 3, 200);  
f = t.^2;
t2 = linspace(0, 4, 200);  

h = figure();
box on; 
grid on; 
hold on; 
plot(t, f, 'b', 'LineWidth', 3)
plot(t2, (2*1.5)*(t2-1.5)+1.5^2, 'c', 'LineWidth', 3)
plot(1.5, 1.5^2, 'ro', 'MarkerSize', 15, 'LineWidth', 3)
text(1, 2.5, '$x_0$','Interpreter', 'latex', 'FontSize', 18)
xlabel('$x$','Interpreter', 'latex', 'FontSize', 18)
ylabel('$f(x)$','Interpreter', 'latex', 'FontSize', 18)
xlim([-3,3])
ylim([-1, 9])
saveas(h, 'pdf/linearization.eps', 'eps2c')


% -------------------------------------------------------------------------
t = linspace(-1, 4.3, 200);  
f = cos(t)+1.5;
fig = figure();
box on; 
grid on; 
hold on; 
plot(t, f, 'b', 'LineWidth', 3)
plot(2, cos(2)+1.5, 'ro', 'LineWidth', 3, 'MarkerSize', 15)
plot(1.75, cos(1.75)+1.5, 'ro', 'LineWidth', 3, 'MarkerSize', 15)
plot(2.25, cos(2.25)+1.5, 'ro', 'LineWidth', 3, 'MarkerSize', 15)
text(2.2, 1, '$f(x_0)$', 'Interpreter', 'latex', 'FontSize', 18)
text(.3, 1.3, '$f(x_0+\delta)$', 'Interpreter', 'latex', 'FontSize', 18)
text(.8, .8, '$f(x_0+\delta)$', 'Interpreter', 'latex', 'FontSize', 18)
xlabel('$x$','Interpreter', 'latex', 'FontSize', 18)
ylabel('$f(x)$','Interpreter', 'latex', 'FontSize', 18)

saveas(fig, 'pdf/linearization_2.eps', 'eps2c')

% -------------------------------------------------------------------------
z = 10.;
h = figure();
box on; 
grid on; 
hold on; 
for i = 1:30
    t = linspace(-5, 0, 200);
    plot(t, 1.2*t+i-10, 'c',  'LineWidth', 3)
end    
plot(4, 0, 'rx', 'MarkerSize', 15, 'LineWidth', 3)  
xlim([-5, 5])
ylim([-5, 5])
xlabel('$Real(s)$', 'Interpreter', 'latex', 'FontSize', 18)
ylabel('$Imag(s)$', 'Interpreter', 'latex', 'FontSize', 18)
saveas(h, 'pdf/poles.eps', 'eps2c')



% -------------------------------------------------------------------------
t_max = 7;
tau = 1.5;

t = linspace(0, t_max, 200);
y = 1./2 - exp(-t/tau)/2;

h = figure();
box on; 
grid on; 
hold on; 
plot(t, .5*ones(numel(t), 1), 'r--', 'LineWidth', 3)
plot([4*tau, 4*tau], [0, 1./2], 'c--', 'LineWidth', 3)
plot(4*tau, 1./2, 'ro', 'MarkerSize', 15)
plot(t, y, 'b', 'LineWidth', 3)
text(4., .52, '$A_1 = y_{FV}$', 'FontSize', 18, 'Interpreter', 'latex')
text(6.2, .05, '$T_s$', 'FontSize', 18, 'Interpreter', 'latex')
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 18)
ylabel('$y(t)$', 'Interpreter', 'latex', 'FontSize', 18)
ylim([0, .555])
saveas(h, 'pdf/first_order_sys.eps', 'eps2c')



% -------------------------------------------------------------------------
num = [1.];
den = [1.,2.];
sys1 = tf(num, den);
[f_step, t_step]  = impulse(sys1);
fig = figure();
box on; 
grid on; 
hold on; 
plot(t_step, f_step, 'b', 'LineWidth', 3)
% legend()
xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 18)
ylabel('y(t)', 'Interpreter', 'latex', 'FontSize', 18)
saveas(fig, 'pdf/impluse_response_example.eps', 'eps2c')


% -------------------------------------------------------------------------
num = [1./2];
den = [1./2,1.];
tau = 1./2;

sys = tf(num, den);
[f_step, t_step]  = step(sys);

fig = figure();
box on; 
grid on; 
hold on; 

plot(t_step, f_step, 'b',  'LineWidth', 3)

text(5*tau+.5, .44, '$A_1 = y_{FV}$', 'FontSize', 18, 'Interpreter', 'latex')
% text(5*tau, .1, '$T_s$', 'FontSize', 18, 'Interpreter', 'latex')

xlabel('$t$', 'Interpreter', 'latex', 'FontSize', 18)
ylabel('$y(t)$', 'Interpreter', 'latex', 'FontSize', 18)
saveas(fig, 'pdf/step_response_example.eps', 'eps2c')


% -------------------------------------------------------------------------
t = linspace(0, 1.5*pi, 100);
w = 5.;
s = 0.75;
f1 = cos(w*t);
f2 = exp(-s*t);

fig = figure();
box on; 
grid on; 
hold on; 
plot(t, f1, 'b',  'LineWidth', 3)
plot(t, f2, 'r',  'LineWidth', 3)
plot(t, -f2, 'r',  'LineWidth', 3)
plot(t, f1.*f2, 'g',  'LineWidth', 3)
xlabel('$t$', 'FontSize', 18, 'Interpreter', 'latex')
saveas(fig, 'pdf/dampingenv.eps', 'eps2c')


% -------------------------------------------------------------------------
zeta = linspace(0.001, .99, 1000);
po = exp(-zeta.*pi./sqrt(1-zeta.^2));

po1 = (po>=.1)*.1;
fig = figure();
box on; 
grid on; 
hold on; 
plot(zeta, po,'b', 'LineWidth', 3)
plot(zeta, po1,'r', 'LineWidth', 3)
xlabel('$\zeta$', 'FontSize', 18, 'Interpreter', 'latex')
ylabel('P.O.', 'FontSize', 18, 'Interpreter', 'latex')
% plt.title('Trade off between $\zeta$ and P.O.')
text(5*tau+.5, .44, '$A_1 = y_{FV}$', 'FontSize', 18, 'Interpreter', 'latex')
saveas(fig, 'pdf/damping.eps', 'eps2c')

%% chapter 8 figures 
clc;
clear; 
close all; 

routh(poly([-2, 1, -1, 1i*5, -1i*5]), .01)
roots([2, 0, 48, 0, -50])

roots([1 0 6 5 8 20])

F = @(K)(3*K.*(K + 2) - 4);
k = linspace(0, 5, 100); 
% D = sqrt(21)/3-1; 
D = sqrt(7/3)-1;

h = figure(); 
box on; 
grid on; 
hold on; 
plot(k, F(k), 'b', 'LineWidth', 3)
plot([D, D], [-20, 400], 'c--', 'LineWidth', 3)
xlabel('$K$', 'FontSize', 18, 'Interpreter', 'latex')
ylabel('$f(K)$', 'FontSize', 18, 'Interpreter', 'latex')
xlim([0, 5])
ylim([-20, 120])
saveas(h, 'pdf/rh_plot.eps', 'eps2c')

%% chapter 9 figures 
clc;
clear; 
close all; 


% -------------------------------------------------------------------------
h = figure; 
box on; 
grid on; 
hold on; 

plot(3, 0, 'ro', 'MarkerSize', 15, 'LineWidth', 3)
plot(5, 0, 'ro', 'MarkerSize', 15, 'LineWidth', 3)
plot(-1, 0, 'rx', 'MarkerSize', 15, 'LineWidth', 3)
plot(-2, 0, 'rx', 'MarkerSize', 15, 'LineWidth', 3)
plot([3, 5], [0, 0], 'c', 'LineWidth', 3)
plot([-2, -1], [0, 0], 'c', 'LineWidth', 3)
saveas(h, 'pdf/rule2_01.eps', 'eps2c')


% -------------------------------------------------------------------------
h = figure; 
box on; 
grid on; 
hold on; 

plot(3, 0, 'ro', 'MarkerSize', 15, 'LineWidth', 3)
plot(5, 0, 'ro', 'MarkerSize', 15, 'LineWidth', 3)
plot(-1, 0, 'rx', 'MarkerSize', 15, 'LineWidth', 3)
plot(-2, 0, 'rx', 'MarkerSize', 15, 'LineWidth', 3)
plot(1, .5, 'rx', 'MarkerSize', 15, 'LineWidth', 3)
plot(1, -.5, 'rx', 'MarkerSize', 15, 'LineWidth', 3)
plot([3, 5], [0, 0], 'c', 'LineWidth', 3)
plot([-2, -1], [0, 0], 'c', 'LineWidth', 3)
saveas(h, 'pdf/rule2_02.eps', 'eps2c')


% -------------------------------------------------------------------------
h = figure; 
box on; 
grid on; 
hold on; 

plot(3, 0, 'ro', 'MarkerSize', 15, 'LineWidth', 3)
plot(5, 0, 'ro', 'MarkerSize', 15, 'LineWidth', 3)
plot(-1, 0, 'rx', 'MarkerSize', 15, 'LineWidth', 3)
plot(-2, 0, 'rx', 'MarkerSize', 15, 'LineWidth', 3)
plot(1, .5, 'rx', 'MarkerSize', 15, 'LineWidth', 3)
plot(1, -.5, 'rx', 'MarkerSize', 15, 'LineWidth', 3)

plot(-3, .5, 'rx', 'MarkerSize', 15, 'LineWidth', 3)
plot(-3, -.5, 'rx', 'MarkerSize', 15, 'LineWidth', 3)
plot(-5, 0, 'ro', 'MarkerSize', 15, 'LineWidth', 3)
plot(-10, 0, 'ro', 'MarkerSize', 15, 'LineWidth', 3)

plot([3, 5], [0, 0], 'c', 'LineWidth', 3)
plot([-2, -1], [0, 0], 'c', 'LineWidth', 3)
plot([-10, -5], [0, 0], 'c', 'LineWidth', 3)
saveas(h, 'pdf/rule2_03.eps', 'eps2c')

k = @(s) abs(s+4).*abs(s+6).*abs(s+10)./abs(s+1);
s = linspace(-10, -6, 200);
K = k(s);
[~,i] = sort(K);

h = figure; 
box on; 
grid on; 
hold on; 
plot(s, K,'LineWidth', 3)
plot([s(i(end)), s(i(end))], [0, K(i(end))], 'k--', 'MarkerSize', 15, 'LineWidth', 3)
plot(s(i(end)), K(i(end)), 'ro', 'MarkerSize', 15, 'LineWidth', 3)
xlabel('$s$', 'FontSize', 18, 'Interpreter', 'latex')
ylabel('$|K(s)|$', 'FontSize', 18, 'Interpreter', 'latex')
saveas(h, 'pdf/rule4_01.eps', 'eps2c')

