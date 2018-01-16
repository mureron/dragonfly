clear all
disp(' ------------four cores------------')
tic;
ParallelConfiguration qq4.txt console slave_open
Primelessthan(3000000);
a(1)=toc;
clear ans
tic;
Primelessthan(3000000);
a(2)=toc;
clear all
disp(' ------------TWO cores------------')
ParallelClear
ParallelConfiguration qq4.txt console slave_open
tic;
Primelessthan(3000000);
a(3)=toc;
clear ans
tic;
Primelessthan(3000000);
a(4)=toc;


