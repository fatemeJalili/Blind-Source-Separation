function X = templateMatrix(f)
t = (1:1250) /250;
k = 1:floor(40/f);
X = [sin(2*pi*f*k.'*t);cos(2*pi*f*k.'*t)];

   