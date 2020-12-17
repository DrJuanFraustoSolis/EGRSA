c Evolutionary Golden Ratio Simulated annealing (EGRSA)
c Copyright (C) 2020  Dr. Juan Paulo Sánchez Hernández, Fanny G. Maldonado Nava, and Dr. Juan Frausto Solis
c Copyright (C) 2005 Frank Eisenmenger, U.H.E. Hansmann, Shura Hayryan, Chin-Ku Hu

c This program is free software; you can redistribute it and/or modify
c it under the terms of the GNU General Public License as published by
c the Free Software Foundation; either version 2 of the License, or (at
c your option) any later version.
c 
c This program is distributed in the hope that it will be useful, but
c WITHOUT ANY WARRANTY; without even the implied warranty of
c MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
c General Public License for more details.
c 
c You should have received a copy of the GNU General Public License
c along with this program; if not, write to the Free Software
c Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
c USA.

c **************************************************************
c
c This file contains the subroutines:  metropolis
c
c Copyright 2003-2005  Frank Eisenmenger, U.H.E. Hansmann,
c                      Shura Hayryan, Chin-Ku Hu
c
c **************************************************************

      subroutine metropolis(eol,currtem,acepta)

c===============================================================
c== SUBROUTINE FOR METROPOLIS UPDATE OF CONFIGURATIONS 	  ==
c==									 	  ==
c== CALLS: energy,addang,(rand),					  ==
c== dummy (function provided as argument)				  ==
c===============================================================


      include 'INCL.H'

	     real*8 e
        dimension aux_vlvr(mxvr)

        
c================================
c== Get Proposal configuration ==
c================================
	call evolutionary()     

         enw = energy()
         delta = enw - eol 

c================================
c== check acceptance criteria ===
c================================
        if (delta.LE.0.0d0) then
          eol = enw
          acepta = acepta + 1.0d0
        elseif (exp(-delta/currtem).GT.rand()) then
          eol = enw
          acepta = acepta + 1.0d0
        else
          vlvr = aux_vlvr
        endif
      return
      end


