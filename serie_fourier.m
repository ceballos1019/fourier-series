clc
clear all
close all

%Graficar funcion original
%%z = 0:0.01:2;
%f1 = (z < 1) .* (z) + (z >= 1);
%plot(z, f1);
%axis([-1 2 -1 2]);
%xlabel('Eje X');
%ylabel('Eje Y');
%ax = gca;
%ax.XAxisLocation = 'bottom';
%ax.YAxisLocation = 'left';

%Calcular serie de fourier
syms t
n = 6;
A = [0 1 2];
f = [t 1];
f = sym(f);
T = max(A) - min(A);
w0 = 2*pi /(T);

A0 = 0;
for i = 1:length(f)
    A0 = A0 + int(f(i), 't', A(i), A(i+1));
end
A0 = simplify(A0/T);

An = 0;
for i = 1:length(f)
    An = An + int(f(i) * cos(n * w0 * t), A(i), A(i + 1));
end
An = simplify(2 * An / T);

Bn = 0;
for i = 1:length(f)
    Bn = Bn + int(f(i) * sin(n * w0 * t), A(i), A(i + 1));
end
Bn = simplify(2 * Bn / T);

%Graphique
%x = linspace(min(A), max(A), 1000);
%fx = (x < 1) .* (x) + (x >= 1);
%for i = 1:length(A) -1
    %{
    if mod(i, 2) == 1
        fx = fx + ((x >= A(i)) & (x <= A(i + 1))) .* subs(f(i), x);
    else
         fx = fx + ((x > A(i)) & (x < A(i + 1))) .* subs(f(i), x);
    end
    %}
     %fx = fx + ((x >= A(i)) & (x < A(i + 1))) .* x;
%end

sum = 0;
for k = 1: n
    ank = subs(An, n, k);
    bnk = subs(Bn, n, k);
    sum = sum + ank * cos(k * w0 * t) + bnk * sin(k * w0 * t);
end

suma = (A0 / 2) + sum;

tv = linspace(0, 2);
y = subs(suma, t, tv);
plot(tv, y)
grid on
hold on

x = linspace(min(A), max(A), 1000);
fx = (x < 1) .* (x) + (x >= 1);
plot(x, fx,'b')
hold off

%{
plot(x, fx, 'b', 'Linewidth', 2); hold on

plot(x + max(x) - min(x), fx, 'b', 'Linewidth', 2)
plot(x - max(x) + min(x), fx, 'b', 'Linewidth', 2)
%plot([max(x) max(x)], [fx(1) fx(end)], 'b', 'Linewidth', 2)
%plot([min(x) min(x)], [fx(end) fx(1)], 'b', 'Linewidth', 2)

grid on
xlabel('\bfTIEMPO');
ylabel('\bfAMPLITUD');
title('\bfGRAFICA DE LA FUNCIÓN');
%}
