function plotsampleald(x)
[~,n] = size(x);
realX = x;
if n>=8
    warning('There more than 7 datasets. The colors may be confused!')
end
color = [0.850,0.325,0.0980;0.929,0.694,0.125;0.494,0.184,0.556;
    0.466,0.674,0.188;0.301,0.745,0.933;0.635,0.0780,0.184;0,0.447,0.741];
xmin = min(x(:)); xmin = xmin-abs(xmin)*0.1;
xmax = max(x(:)); xmax = xmax+abs(xmax)*0.1;
figure
for i = 1:n
    ci = mod(i-1,7)+1;
    x = realX(:,i);
    opt.display = 0;
    [alpha,lambda,kappa] = mleald(x,opt);
    hold on
    if n==1;
        ci = 7;
    end
     histogram(x,'BinMethod','fd','EdgeAlpha',0.3,...
         'Normalization','pdf','FaceColor',color(ci,:),'NumBins',80)
%     histogram(x,60,'EdgeAlpha',0.3,'Normalization','pdf','FaceColor',color(ci,:))
    h = plotpdfald(alpha,lambda,kappa,xmin,xmax);
    if n>=2
        h.Color = color(ci,:);
    end
    hold off
end