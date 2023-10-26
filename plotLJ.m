function plotLJ(a)
  b=a';
  plot(b(1,1:end),b(2,1:end),"-;KE;",b(1,1:end),b(3,1:end),"-;PE;",b(1,1:end),b(4,1:end),"-;KE+KP;",b(1,1:end),b(5,1:end),"-;T;");
  xlabel ("t");
endfunction