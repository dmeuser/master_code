***************************************
*     Example of fortran code
***************************************
      program main
      implicit none
      double precision nnxfx
      
      call initnnset("PDF4LHC15_nlo_mc_pdfas"//char(0), 101)
      write(6,*) nnxfx(0.1d0,100.0d0,5)
      
      end
