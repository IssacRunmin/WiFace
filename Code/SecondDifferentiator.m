function Re = SecondDifferentiator(x, h)
% used in app_data_processing_v8.m
% get the accleration (second differentiator) of the signal
if (length(x) < 8)
    Re = [];
    return;
end
% f''_0 = 4f_0+(f_1+f_(-1))-2(f_2+f_(-2))-(f_3+f_(-3))
Re = 4.*x(4:end-3)+(x(5:end-2)+x(3:end-4))...
    -2.*(x(6:end-1)+x(2:end-5))-(x(7:end)+x(1:end-6));
Re = Re ./ (16*h*h);
if isrow(Re)
    Re = Re';
end
Re = cat(1,Re(1:3),Re,Re(end-2:end));