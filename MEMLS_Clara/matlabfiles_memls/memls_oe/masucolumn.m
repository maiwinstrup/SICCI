function a=masucolumn(a)
 [ya,xa]=size(a);
 if (xa>1)
  a=a';
 end %endif

end %endfunction
