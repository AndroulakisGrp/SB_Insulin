% This should be turned into a function, and probably moved to the +plot
% folder.

fig=figure;
set(fig,'DoubleBuffer','on');
n=60;
max_offset = find(min(t_end(t_end>71*60*60))==t_end); % 48 hours

%%
% set(gca,'NextPlot','replacechildren')
% axis equal % fix the axes
% mov = avifile('poincare_lps.avi', 'FPS', 2);
% for i=0:n
%     offset = round(max_offset*i/n);
%     plot_poincare(RRintervals(1+offset:offset+1000));
%     axis([0, 1.5, 0, 1.5]);
%     F=getframe(gcf);
%     mov = addframe(mov,F);
% end
% mov = close(mov);
% % movie(A,10,3) % Play the MATLAB movie 

%%
plot_poincare(RRintervals(1:1000));
axis([0, 1.5, 0, 1.5]);
f = getframe;
[im,map] = rgb2ind(f.cdata,256,'nodither');
im(1,1,1,n) = 0;
for i = 1:n
    offset = round(max_offset*i/n);
    plot_poincare(RRintervals(1+offset:offset+1000));
    axis([0, 1.5, 0, 1.5]);
    f = getframe;
    im(:,:,1,i) = rgb2ind(f.cdata,map,'nodither');
end
imwrite(im,map,'poincare_baseline.gif','DelayTime',0,'LoopCount',inf) %g443800