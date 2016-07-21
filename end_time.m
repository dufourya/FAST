function t = end_time(x)
if size(x,1)>4
t = x(end,11);
else
    t=0;
end
end