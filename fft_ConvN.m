function [l_grid_padding,fLConvFFT_array] = fft_ConvN(l_grid, fL_array,composition)

l_grid_padding = composition*min(l_grid):(l_grid(2)-l_grid(1)):composition*max(l_grid);

% We can actually avoid all FFT shifting (and normalization) if we zero-pad at the end of 
% the vector, rather than on both sides. In fact the fft function will automatically do this
% if we just add a second argument with the desired length
fFFT = fft(fL_array,length(l_grid_padding));
cfFFT_array = fFFT.^composition;
fLConvFFT_array = ifft(cfFFT_array);
fLConvFFT_array = max(0,fLConvFFT_array);
end