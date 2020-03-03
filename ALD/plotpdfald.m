function h = plotpdfald(alpha,lambda,kappa,xmin,xmax)
if nargin == 3
    xmin = -5;
    xmax = 5;
end
series = linspace(xmin,xmax,500);
y = pdfald(series,alpha,lambda,kappa);


h = plot(series,y,'LineWidth',2);