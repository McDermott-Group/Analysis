function save2PDF(h, filename)
%SAVE2PDF Saves a figure specified by handle h to a pdf file named filename.

h.PaperPositionMode = 'auto';
fig_pos = h.PaperPosition;
h.PaperSize = [fig_pos(3) fig_pos(4)];
print(h, filename, '-dpdf', '-r0')

end