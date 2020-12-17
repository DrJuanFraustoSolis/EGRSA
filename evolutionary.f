c *************axvr*************************************************
c
c This file contains the subroutines:  evolutionary
c
c Copyright 2018-2020  Juan Paulo Sánchez Hernández, Fanny G.
c                      Maldonado Nava, and Juan Frausto Solis
c                      
c This genetic algorithm was implemented by the perturbation process
c of golden ratio simulated annieling. We generate a set of structures and we
c apply crossover and mutation to obtain the best structure.
c **************************************************************
      subroutine evolutionary()

       include 'INCL.H'

       dimension pob_vlvr(10,mxvr)
       dimension a_pob_vlvr(10,mxvr)
       dimension a_pob_energy(10)
       dimension aux_vlvr(mxvr)
       dimension pob_energy(10)
       integer i,j,k,p,Tpob,ind1,ind2,Pc,Pm,e_min
       integer Ngen,minval
       real*8 e

        Tpob = 10
        Ngen = 1

c       Generate initial population 
       do i=1, Tpob		 	
            aux_vlvr = vlvr   
            jv = idvr(1+int(nvr*rand()))
            e = rand()
            dv = axvr(jv)*e 
            vlvr(jv) = addang(vrol,dv)
            pob_energy(i) = energy()
            pob_vlvr(i,:) = vlvr
            vlvr = aux_vlvr
       enddo
 
c      Genetic algorithm

      do i=1, Ngen

c     ----Tournament selection---
         do j=1,Tpob

           ind1 = 1+int(10*rand())
           ind2 = 1+int(10*rand())

           if (pob_energy(ind1).le.pob_energy(ind2)) then
             a_pob_vlvr(j,:) = pob_vlvr(ind1,:)
             a_pob_energy(j)= pob_energy(ind1)
           else
             a_pob_vlvr(j,:) = pob_vlvr(ind2,:)
             a_pob_energy(j)= pob_energy(ind2)           
           endif

         enddo


           pob_energy = a_pob_energy
           pob_vlvr = a_pob_vlvr
           p=1 



c     ---Crossover of popolation---
        do j=1,Tpob/2

         do while (ind1.eq.ind2)
            ind1 = 1+int(10*rand())
            ind2 = 1+int(10*rand())
         end do
     
         Pc = 1+int(nvr*rand())

         a_pob_vlvr(p,1:Pc) = pob_vlvr(ind1,1:Pc)
         a_pob_vlvr(p+1,1:Pc) = pob_vlvr(ind2,1:Pc)

         a_pob_vlvr(p,Pc+1:mxvr) = pob_vlvr(ind2,Pc+1:mxvr)
         a_pob_vlvr(p+1,Pc+1:mxvr) = pob_vlvr(ind1,Pc+1:mxvr)

         p=p+2
        enddo


      pob_vlvr = a_pob_vlvr

c     ---Mutation of population---

      do j=1,Tpob
        aux_vlvr = pob_vlvr(j,:)   
        jv = idvr(1+int(nvr*rand()))
        e = rand()
        dv = axvr(jv)*e 
        pob_vlvr(j,jv) = addang(aux_vlvr,dv)
        vlvr = pob_vlvr(j,:)
        pob_energy(j) = energy()
      enddo


      enddo
      
c      Select the best structure
        e_min= minloc(pob_energy,DIM=1)
        vlvr = pob_vlvr(e_min,:)

      end
