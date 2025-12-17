% The script is an example for the function "inset"
% figure 1 is a plot of 'zero order bessel function of the first kind' in the range of x=0 to x=50.
% figure 2 is a close-up of this function around its first zero.
% figure 3 shows the function in the main plot, and the close-up in the
% inset plot.
% 
% Moshe Lindner, August 2010 (C).

close all

fig1=open('images/pt_ef_ripetizione_EC.fig');
fig2=open('images/pt_ef_ripetizione_EC_zoom.fig');
[h_m h_i]=inset(fig1,fig2);

set(h_i,'xlim', [2e-5,8e-5], 'ylim',[10e-11,10e-10])

main_fig = findobj(h_m,'Type','axes');
saveas(main_fig, 'images/pt_ef_ripetizione_EC_inset.fig');