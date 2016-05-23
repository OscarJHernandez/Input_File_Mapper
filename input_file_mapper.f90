! This program will read in the .ho input file and will generate an output file that is easily readable by 
! my DD_Code program


program input_file_mapper
implicit none
integer::NMax,jjMax,jjfixed,ssfixed
integer::lParity_restrict,tMax,ssMax,info
integer:: ssMax_0,jjMax_0,tMax_0
integer,allocatable:: Q_table(:,:),Q_table2(:,:)
integer::st,matrixE_oper,matrixE
integer,external::dimensionOfMatrix
real*8:: hw
real*8,allocatable:: H_large(:,:),H(:,:)
character*50::input,word,word2,fname,description
character*30:: potential

    ! Now we open the input from the file: H_input.txt
    fname ="H_input.txt"
    open(unit =1, file= fname,iostat=st)
    
    ! Reads from an input file
    read(1,*) word, Nmax
    read(1,*) word,jjMax
    read(1,*) word,input
    read(1,*) word,word2
    read(1,*) word,word2

    read(1,*) word,word2
    read(1,*) word,word2
    read(1,*) word,word2
    read(1,*) word,word2
    read(1,*) word,word2
    read(1,*) word,word2
    read(1,*) word,info

    close(1)    
     
    ! The spin will either be 0 or 1, T=0,1
    ssMax =2
    tMax =1 
    
    !===============================================================================================
    ! We generate a matrix for the Ground state of the system
    
    ! These parameters are fixed for the deuteron ground state
    ssMax_0=2  ! The Spin can be either 0->1, but will be restricted to just s=1 by ssfixed
    jjMax_0=2 ! jj=2, which means j=1 
    lParity_restrict=1  ! The parity of the angular component of the wave function is even
    ssfixed=2   ! The spin of the deuteron is fixed to s=1
    jjfixed=2   ! The total angular momentum j of the deuteron is fixed to j=1
    tMax_0=0
    
    ! Calculate the total number of elements in the Quantum Number table
    matrixE = dimensionOfMatrix(matrixE,Nmax,jjMax_0,ssMax_0,tMax_0,jjfixed,ssfixed,lParity_restrict) 
    
    ! Allocate memory for the table of Quantum numbers 
    allocate(Q_table(matrixE,6))
    
    
    print *, 'Nmax: ',Nmax
    print *, 'jjMax: ', jjMax
    print *, 'hw: ', Int(hw)
    print *, 'ssmax: ', ssMax
    
    ! Fill the Q_table with quantum numbers with Ground state restrictions
    call generateTable(Q_table,matrixE,Nmax,jjMax_0,ssMax_0,tMax_0,jjfixed,ssfixed,lParity_restrict)

    
    ! Allocate memory for all the matrices
    allocate(H(matrixE,matrixE))
    
    CALL extractMatrixElements(H,input,Q_table,matrixE,hw,info,potential)
    
    description='GS'
    call print_matrix_elements(H,matrixE,hw,Nmax,description,potential)
    
    !===============================================================================================
    ! We generate a matrix for the the Excited states of the system
    !
    ! Here we tell the program to vary jj and ss
    jjfixed=-1
    ssfixed= -1
    lParity_restrict= 0 ! No restiction on the parity of states other than they have to be fermionic
     
    ! Now we calculate the dimension of the Matrix 
    matrixE_oper = dimensionOfMatrix(matrixE_oper,Nmax,jjMax,ssMax,tMax,jjfixed,ssfixed,lParity_restrict) 
    
    print *, 'Nmax', Nmax
    print *, 'Input File: ', input
    
    
    print *,'Size of Matrix: ', matrixE_oper
    
    ! Allocate memory for the quantum number table 
    allocate(Q_table2(matrixE_oper,6))
        
    ! We Generate the table of new Quantum Numbers (Q_table)
    call generateTable(Q_table2,matrixE_oper,Nmax,jjMax,ssMax,tMax,jjfixed,ssfixed,lParity_restrict)
    

     
    allocate(H_large(matrixE_oper,matrixE_oper))
        

    ! Initiliaze the Matrices
    H_large(:,:)=0.0
    

    ! The program now extracts matrix Elements from a file
    call extractMatrixElements(H_large,input,Q_table2,matrixE_oper,hw,info,potential)
    
    description='ES'
    ! Now we write all matrix elements to a file
    call print_matrix_elements(H_large,matrixE_oper,hw,Nmax,description,potential)
    
    print *, 'End of Program'

end program input_file_mapper


! H is a symmetric Matrix
subroutine print_matrix_elements(H,n,hw,Nmax,description,potential)
implicit none
integer::n,i,j,Nmax
real*8::hw
real*8::H(n,n)
character*50:: outputfile,string,string2,description
character*30::potential


write(string,"(I3.3)") Int(hw)
write(string2,"(I3.3)") Nmax
outputfile = "H_matrix"//"_"//trim(potential)//"_"//trim(description)//"_hw"//trim(string)//"_Nmax"//trim(string2)//".txt"
open(11,file=outputfile,status='unknown')

write(11,*) Nmax,hw

do i=1,n
  do j=i,n
write(11,*) i,j,H(i,j)
  end do
end do

close(11)

end subroutine print_matrix_elements


! This function calculates the number of total matrix elements
integer function dimensionOfMatrix(matrixE,nMax,jjMax,ssMax,tMax,jjfixed,ssfixed,lParity_restrict) 
implicit none
integer:: matrixE ! Total number of matrix Elements
integer:: nMax ! Maximum order of Laguerre Polynomial
integer:: jjMax ! Maximum total angular momentum value
integer::ssMax ! maximum Spin Value
integer:: tMax ! Maximum Isospin Value
integer:: jjFixed ! Fixed value of jj
integer:: ssFixed ! Fixed value of ss
integer:: lParity_restrict ! Restricted parity of l
integer::t,ss,jj,n,l,ll ,mt
integer::d
integer::sparity,lparity,tparity,parity

  ! Initialize this value
   matrixE=0
   d=0
       
       ! SS and jj are not fixed, assume no parity restrictions
       if((ssfixed.eq.-1).and.(jjfixed.eq.-1)) then
                
                t_loop: do t=0,tMax
                  mt_loop: do mt=-t,t
                    j_loop: do jj=0,jjMax,2   
                       ss_loop: do ss=0,ssMax,2
                           ll_loop: do ll = abs(jj-ss),abs(jj+ss),2
                                
                                ! We require 2*n + l <= Nmax
                                l =ll/2 
                                n = 0                       
                                
                                ! n_loop: do n=0,nHighest 
                                 n_loop: do while(n <= (Nmax-l)/2 )
                                   
                                    lparity = (1-2*mod(ll/2,2)) ! (-1)**l
                                    tparity = 2*t-1
                                    sparity = ss-1
                                    parity = lparity*tparity*sparity
                                    
                                    ! lParity_restrict = 0 means no parity restriction
                                    if(lParity_restrict.eq.0) then
                                        
                                        ! Choose only fermionic states
                                        if((parity.eq.-1)) then  
                                        matrixE=matrixE+1
                                        end if
                                    ! If there is a parity restriction
                                    else
                                        
                                        ! Choose fermionic states that obey l_parity restriction
                                        if((parity.eq.-1).and.(lparity.eq.lParity_restrict)) then  
                                        matrixE=matrixE+1
                                        end if
                                    
                                    end if
                                    
                                    ! Increment n
                                    n=n+1
                                    
                              end do n_loop
                           end do ll_loop
                        end do ss_loop
                    end do j_loop
                  end do mt_loop
                end do t_loop
      
       ! If ss is fixed and jj is fixed assume parity restrictions
      else if((ssfixed.ne.-1).and.(jjfixed.ne.(-1))) then
       

                        ! Loop over all possible quantum numbers
                        ! Assumes that mj is also fixed
                   t_loop2: do t=0,tMax
                     mt_loop2: do mt=-t,t
                       ll_loop2: do ll=abs(ssfixed-jjfixed),abs(ssfixed+jjfixed),2
                               
                                l = ll/2
                                n=0
                                
                               n_loop2: do while(n <= (Nmax-l)/2 )   
                                    
                                    lparity = (1-2*mod(ll/2,2))
                                    tparity = 2*t-1
                                    sparity = ssfixed-1
                                    parity = lparity*tparity*sparity
                                    
                                    if((parity.eq.-1).and.(lparity.eq.lParity_restrict)) then
                                    matrixE=matrixE+1
                                    end if 
                                 
                                 ! Increment n
                                 n = n+1
                               end do n_loop2
                               
                               
                        end do ll_loop2
                        end do mt_loop2
                     end do t_loop2  
           
    end if ! ends all ifs

 dimensionOfMatrix =matrixE

end function dimensionOfMatrix

! This subroutine fills the Q_table with all quantum numbers
subroutine generateTable(Q_table,matrixE,Nmax,jjMax,ssMax,tMax,jjfixed,ssfixed,lParity_restrict)
implicit none
integer:: matrixE ! Total number of matrix Elements
integer:: Nmax ! Maximum order of Laguerre Polynomial
integer:: jjMax ! Maximum angular momentum value
integer::ssMax ! maximum Spin Value
integer:: tMax ! Maximum Isospin Value
integer:: jjFixed ! Fixed value of jj
integer:: ssFixed ! Fixed value of ss
integer:: lParity_restrict ! Restricted value of l
integer::Q_table(matrixE,6) ! Stores index i, and 6 other quantum numbers
integer::n,l,i,t,jj,ss,ll,mt
integer::sparity,lparity,tparity,parity
              
      ! Initialize the counting 
           i=1
                                 
      ! SS and jj are not fixed
      ! Assume no restrictions on l
       if((ssfixed.eq.-1).and.(jjfixed.eq.-1)) then
       
       

                t_loop: do t=0,tMax
                mt_loop: do mt=-t,t
                    jj_loop: do jj=0,jjMax,2
                        ss_loop: do ss=0,ssMax,2
                          ll_loop: do ll=abs(jj-ss),abs(jj+ss),2       
                            
                                l =ll/2 
                                n = 0  
                                 
                            ! We require 2*n + l <= Nmax
                             n_loop: do while(n <= (Nmax-l)/2 )
                                  
                                    
                             ! n_loop: do n=0,nMax-d
                                    lparity = (1-2*mod(l,2))
                                             tparity = (2*t-1)
                                             sparity = (ss-1)
                                             parity = lparity*tparity*sparity
                                             
                                             
                                    ! lParity_restrict = 0 means no parity restriction
                                    if(lParity_restrict.eq.0) then
                                        
                                        ! If we have a fermionic state
                                        if(parity.eq.-1) then
                                              Q_table(i,1)=n
                                              Q_table(i,2)=l
                                              Q_table(i,3)=ss
                                              Q_table(i,4)=jj
                                              Q_table(i,5)=t
                                              Q_table(i,6)=mt

                                               i=i+1
                                        end if
                                    ! else we have a parity restriction 
                                    else
                                        
                                        ! Obeys parity restriction on l, and its a fermionic state
                                        if((parity.eq.-1).and.(lparity.eq.lParity_restrict)) then  
                                              Q_table(i,1)=n
                                              Q_table(i,2)=l
                                              Q_table(i,3)=ss
                                              Q_table(i,4)=jj
                                              Q_table(i,5)=t
                                              Q_table(i,6)=0

                                               i=i+1
                                        end if
                                    
                                    end if
                                    
                                    ! Increment n
                                    n =n+1
                                    
                                       end do n_loop         
                                    end do ll_loop
                                  end do ss_loop
                               end do jj_loop
                      end do mt_loop
                    end do t_loop
       
       ! If ss is fixed and jj is fixed
      else if((ssfixed.ne.-1).and.(jjfixed.ne.(-1))) then
       
                        ! Loop over all possible quantum numbers
                        ! Assumes mj is not fixed
                     t_loop2: do t=0,tMax  
                     mt_loop2: do mt=-t,t 
                        ll_loop2: do ll=abs(ssfixed-jjfixed),abs(ssfixed+jjfixed),2
                               
                                l =ll/2 
                                n = 0  
                               
                               n_loop2: do while(n <= (Nmax-l)/2 )
                                        
                                         lparity = (1-2*mod(ll/2,2))
                                         tparity = 2*t-1
                                         sparity = (ssfixed-1)
                                         parity = lparity*tparity*sparity
                                        
                                        
                                        if((parity.eq.-1).and.(lparity.eq.lParity_restrict)) then
                                           
                                              Q_table(i,1) =n
                                              Q_table(i,2)=(ll/2)
                                              Q_table(i,3)=ssfixed
                                              Q_table(i,4)=jjfixed
                                              Q_table(i,5)=t
                                              Q_table(i,6)=0
        
                                               i=i+1   
                                        end if 
                                        
                                        ! Increment n
                                        n = n+1
                                        
                                end do n_loop2
                        end do ll_loop2
                        end do mt_loop2
                    end do t_loop2
           
    end if ! ends all ifs
      
      
       

end subroutine generateTable

! This function looks up the index corresponding to a set of quantum numbers
! Returns -1 if the set of quantum numbers is not located in the table
integer function retrieveIndex(n1,l1,ss1,jj1,t1,mt1,Q_table,matrixE)
implicit none
integer::n1,l1,ss1,jj1,t1
integer::n,l,ss,jj,t
integer:: i,mt,mt1
integer::matrixE ! The total number of matrix Elements
integer::Q_table(matrixE,6)

retrieveIndex=-1

! Search the entire Quantum number table for 
    do i=1,matrixE
        
                           n=Q_table(i,1) 
                           l=Q_table(i,2)
                           ss=Q_table(i,3)
                           jj=Q_table(i,4)
                           t=Q_table(i,5)
                           mt=Q_table(i,6)
                           
        
        
        ! Go throught the table until a match is found
        if((n1.eq.n).and.(l1.eq.l).and.(ss.eq.ss1).and.(jj.eq.jj1).and.(t.eq.t1).and.(mt.eq.mt1)) then
            retrieveIndex = i
            return
        end if
        
        
    end do

return 
end function retrieveIndex


! This program extracts all of the required matrix elements from the code and generates a 2D matrix
! The program goes through the total number of matrix Elements
subroutine extractMatrixElements(Vmatrix,fname,Q_table,matrixE,hw,info,potential)
implicit none
real*8:: hw
integer::n_matrixElements,order
integer,intent(in)::matrixE,info
character*30,intent(in)::fname
integer,intent(in)::Q_table(matrixE,6)
integer::mt
integer:: n1,n2,l1,l2,st,jMax,isospin_channels,N_x,N_y,N
integer::ll1,ll2,jj,ss,i,x,y,t
real*8::m1,m2,m3,m4,KE
integer,external:: retrieveIndex
real*8,intent(out)::Vmatrix(matrixE,matrixE)
character*30::potential


! Initialize the Matrix and some Parameters:
Vmatrix(:,:)=0.d0
N_x=0
N_y=0


!Open the File again
open(unit =1, file= fname,iostat=st)
read(1,*) potential
read(1,*) hw, order, jMax
read(1,*) n_matrixElements, isospin_channels

 print *, potential
 print *,'hw = ', hw
      

! Initialize Counting:
N=0
        do i=1,n_matrixElements
            
            if(info.eq.1) then
            ! Give Some information about progress:
            write(*,FMT="(A1,A,t21,F6.2,A)",ADVANCE="NO") achar(13), &
            & "Constructing Matrix:  ", (real(i)/real(n_matrixElements))*100.0 , "%"
            end if
            
            
            ! Line 1: 2*j,2*s
            ! Line 2: n1 n2 l1 l2
            ! Line 3: m1 m2 m3 m4
            read(1,*) jj,ss
            read(1,*) n1,n2,ll1,ll2
            read(1,*) m1,m2,m3,m4
            
            l1 =ll1/2
            l2 =ll2/2
            
            do t=0,1
            do mt=-t,t
                                  
                    ! Retrieve the index from the table of quantum numbers
                    x=retrieveIndex(n1,l1,ss,jj,t,mt,Q_table,matrixE)
                    y=retrieveIndex(n2,l2,ss,jj,t,mt,Q_table,matrixE)
            
            
                        ! If the set of quantum numbers is in the Table
                        ! Note that Nir uses the HO basis states with an extra phase of (-1)^(n1+n2), so we multiply by the same factor to
                        ! eliminate the extra phase factor
                        if((x.ne.-1).and.(y.ne.-1)) then
                             
                            KE=0.d0
                            
                            if((l1.eq.l2)) then
                                  if (n1.eq.n2-1) KE=0.5d0*dsqrt(DBLE(n2)*(DBLE(n2)+DBLE(l1)+0.5d0))
                                  if (n1.eq.n2  ) KE=0.5d0*(2.d0*DBLE(n1)+DBLE(l1)+1.5d0)
                                  if (n2.eq.n1-1) KE=0.5d0*dsqrt(DBLE(n1)*(DBLE(n1)+DBLE(l1)+0.5d0))          
                            end if
                            
                            
                            ! We take T=0, mt=0 Channel
                            if((t.eq.0).and.(mt.eq.0)) then 
                                
                                Vmatrix(y,x)=KE*hw+m1*(1-2*mod(n1+n2,2))
                                

                                N=N+1
                             ! We take T=1, mt=-1 channel
                            else if((t.eq.1).and.(mt.eq.-1)) then
                                
                                Vmatrix(y,x)=KE*hw +m2*(1-2*mod(n1+n2,2))  !*(1-2*mod(l1+l2+1,2))

                                N=N+1

                            ! We take T=1, mt=0 channel
                            else if((t.eq.1).and.(mt.eq.0)) then
                                
                                Vmatrix(y,x)=KE*hw +m3*(1-2*mod(n1+n2,2))  !*(1-2*mod(l1+l2+1,2))

                                N=N+1
                           ! write(2,*),y,x, m3*(1-2*mod(n1+n2,2))
                            else if((t.eq.1).and.(mt.eq.1)) then
                                
                                Vmatrix(y,x)=KE*hw +m4*(1-2*mod(n1+n2,2))  !*(1-2*mod(l1+l2+1,2))

                                N=N+1
                           end if
                        
                        end if
                    
                    
                end do
                end do
           
   
        end do
        
print *, ' '        
print *,'-----------------------------------------------'
print *, 'Data from Matrix Element Extraction: '
print *, 'input File: ' , fname
print *, 'Matrix Elements in Input File: ', n_matrixElements
print *, 'Vector Length: ',matrixE
print *, 'Number of Extracted Matrix Elements: ', N
print *, 'Number of Matrix Elements Required for Hamiltonian', matrixE**(2)
print *,'-----------------------------------------------'
print *,'-----------------------------------------------'
! Close the input file
close(1)

        

end subroutine extractMatrixElements
