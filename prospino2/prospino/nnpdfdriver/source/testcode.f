***************************************
*     Example of fortran code
***************************************
      function NNPDF(x,Q,i)
      Implicit Double Precision (A-H,O-Z)
      
      NNPDF = nnxfx(x,Q,i)
c~       write(6,*) nnxfx(0.1d0,10d0,0)
      
      end
