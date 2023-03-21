%{
Run the lookup table simulation

Benjamin Kao
Last update: 2020/10/26
%}

clc;clear;close all;

%% clear the stop flag
to_save=0;
save('stop_flag.txt','to_save','-ascii','-tabs');

fun_MCX_run_lookup('CS');
% fun_MCX_run_lookup('BY');
