% Klausur 2017
% Ex2) 
format long g
x = [0;1;-1].';% coefficients of (delayed) input x
y = [4;-2;2;-1].';% coefficients of (delayed) output y
mySys = tf(x,y,-1,'variable','z^-1');
[poles,zeros] = pzmap(mySys);% plot pzmap with "pzmap(mySys)"
linearSystemAnalyzer({'pzmap';'impulse'},mySys);% new
h = impulse(mySys,4);
disp('Impulse response h[n], 0 <= n <= 4')
disp(h)
% Stable + IIR

% plot output of x[n]=(0.5)^n sin(2n)
figure;
N = 100;
h = impulse(mySys,20);
% n: time axis
% x: input signal
% y: output signal
% p: handle for plot
nn = (0:1:N);
p1 = subplot(2,1,1); plot(h);

x = cos(2.*nn);
y = conv(x,h);
p2 = subplot(2,1,2); plot(y);

linkaxes([p1, p2],'y')

%% Question 4
x = [1;2;];
w = [1;-1];
y_lin = conv(x,w);
y_cc = cconv(x,w,2);

disp('Linear Convolution')
disp(y_lin)
disp('Circular Convolution')
disp(y_cc)

%%
x = [1;2;-1;4];
w = [0,0,1,-1];
y_lin = conv(x,w);
y_cc = cconv(x,w,4);

disp('Linear Convolution')
disp(y_lin)
disp('Circular Convolution')
disp(y_cc)