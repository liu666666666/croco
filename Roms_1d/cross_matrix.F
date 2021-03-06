! $Id: cross_matrix.F 985 2012-12-17 15:07:30Z gcambon $
!======================================================================
! ROMS_AGRIF is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al) 
! and Rutgers University (Arango et al) are under MIT/X style license.
! ROMS_AGRIF specific routines (nesting) are under CeCILL-C license.
! 
! ROMS_AGRIF website : http://www.romsagrif.org
!======================================================================
! This is "cross_matrix": analyzer of dependencies to build
!                         table of dependencies "Make.depend"
!------------------------------------------------------------------
! Usage:
!--------       cross_matrix list_of_files
! for example
!               cross_matrix *.F
! or
!               cross_matrix *.F *.h
! or
!               cross_matrix $(SRCS)            inside a makefile
!
! Algorithm:
!-----------
! cross_matrix reads all files from the list arguments it is given,
! and searches for the CPP command #include. If #include commands 
! are found, in interprets names within the "..." as new files to be
! open and searched for #include inside it and so on.
!
! cross_matrix investigates chained include statements, for example,
! if file1 contains statement #include "file2", and file2 contains
! statement #include "file3", therefore file3 will be included into
! list of dependencies of "file1" 
!
! Optionally, files which are not in the list of arguments AND
! which are not in the present working directory can be excluded
! from the list of dependencies (CPP switch EXCLUDE). A warning
! message will be issued about each file excluded from analysis. 
!
#define EXCLUDE

      implicit none
      integer length
      parameter (length=72)
      character*1 symbol(length), quote, double_quote
      integer last_arg,iargc, iin,iout,iocheck, i,j,k,l,n 

      parameter(iin=11, iout=12)

      integer nmax,n_disc,n_excl, max_name_lenght,max_names
      parameter (max_name_lenght=32, max_names=256)
      integer nlength(max_names), l_count
      character*32 name_string(max_names)
      character*1 name(max_name_lenght,max_names)
      equivalence(name,name_string)
      logical matrix(max_names,max_names),NEW_NAME, END_OF_FILE
!-
      last_arg=iargc()
      if (last_arg.eq.0) then
        write(*,*) 'ERROR IN cross_matrix: NO FILE NAMES ARE GIVEN'
        stop
      endif
!-
      quote=char(39)
      double_quote=char(34)
      n_disc=0
      n_excl=0
      l_count=0
      do j=1,max_names
        do i=1,max_names
          matrix(i,j)=.false.
        enddo
      enddo
!-
      do n=1,last_arg
        do i=1,max_name_lenght
          name(i,n)=' '
        enddo
        call getarg(n,name_string(n))
        i=max_name_lenght+1
 101     i=i-1
         if (name(i,n).eq.' ' .and. i.gt.1) goto 101
        nlength(n)=i
      enddo
      nmax=last_arg
!-
      open(unit=iout, file='Make.depend', form='FORMATTED',
     &                                    status='UNKNOWN')
      write(iout,'(A22,A40)') '# Make.depend: list of',
     &              ' dependencies generated by cross_matrix.' 
      write(iout,'(A26,A34)') '# WARNING: THIS IS MACHINE',
     &                      ' GENERATED FILE: DO NOT EDIT !!!!!'
      write(iout,'(A31,A14)') '# To create or update this file',
     &                                         ' use commands:'
      write(iout,'(A1 )') '#'
      write(iout,'(A25)') '#        cross_matrix *.F'
      write(iout,'(A4 )') '# or'
      write(iout,'(A29)') '#        cross_matrix *.F *.h'
      write(iout,'(A4 )') '# or'
      write(iout,'(A20)') '#        make depend' 
      write(iout,'(A1 )') '#'
!-
      n=0
  1    n=n+1
       END_OF_FILE=.false.
       open(unit=iin, file=name_string(n), form='FORMATTED',
     &                                status='OLD', err=102)
        goto 2
 102    continue
#ifdef EXCLUDE
        n_excl=n_excl+1
        write(iout,'(A1)') '#'
        write(iout,'(A22,A9,32A1)') '# WARNING: File is not',
     &                  ' found:  ',(name(i,n),i=1,nlength(n))
        write(iout,'(A23,A26)')     '# This file is excluded',
     &                         ' from the dependency list.'
        do i=1,nlength(n)
          name(i,n)=' '
        enddo
        nlength(n)=0
        goto 1
#endif
  2     symbol(1)= ' '
         read(iin,'(80A1)',iostat=iocheck,end=3)
     &                     (symbol(i), i=1,length)
         l_count=l_count+1
        if (iocheck.eq.0 .and. symbol(1).ne.'#') goto 2
        goto 4
  3     END_OF_FILE=.true.
  4     i=0
  5      i=i+1
         if (symbol(i).ne.'i') then
           if (i.lt.length) goto 5
         else
           if (symbol(i+1).eq.'n' .and. symbol(i+2).eq.'c' .and.
     &         symbol(i+3).eq.'l' .and. symbol(i+4).eq.'u'  .and.
     &         symbol(i+5).eq.'d' .and. symbol(i+6).eq.'e')  then
             j=i+6
  6           j=j+1
              if (symbol(j).ne.double_quote) then
                if (j.lt.length) goto 6
              else
                k=j
  7              k=k+1
                 if (symbol(k).ne.double_quote) then
                   if (k.lt.length) goto 7
                 else
                   nmax=nmax+1
                   nlength(nmax)=k-j-1
                   do l=1,nlength(nmax)
                     name(l,nmax)=symbol(l+j)
                   enddo
                   do l=nlength(nmax)+1,max_name_lenght
                     name(l,nmax)=' '
                   enddo
                   NEW_NAME=.true.
                   do k=1,nmax-1
                     if (nlength(k).eq.nlength(nmax)) then
                       do l=1,nlength(nmax) 
                         if (name(l,k).ne.name(l,nmax)) goto 8
                       enddo
                      NEW_NAME=.false.
                      matrix(k,n)=.true.
                     endif
   8                 continue
                   enddo
                   if (NEW_NAME) then
                     n_disc=n_disc+1
                     matrix(nmax,n)=.true.
                   else
                     nmax=nmax-1
                   endif
                 endif
              endif
           endif
         endif                
        if (.not.END_OF_FILE) goto 2
       close (unit=iin)
      if (n.lt.nmax) goto 1
!-  
      l=0                            ! Investigate the secondary  
  10  i=l                            ! dependencies. This is
      do n=1,nmax                    ! logical multiplication:
        do k=1,nmax                  ! if 
          if (matrix(k,n)) then      !      file1 depends on file2
            do j=1,nmax              ! and
              if (matrix(j,k)) then  !      file2 depends on file3
                matrix(j,n)=.true.   ! then
              endif                  !      file1 depends on file3
            enddo                    !
          endif                      ! This is an iterative
        enddo                        ! procedure, since staged 
      enddo                          ! include statements are 
      l=0                            ! possible. The number of .true.
      do n=1,nmax                    ! elements in the matrix grows
        do k=1,nmax                  ! when new dependencies are
          if (matrix(k,n)) l=l+1     ! discovered. The procedure
        enddo                        ! repeats itself until the next
      enddo                          ! iteration does not discover
      if (l.ne.i) go to 10           ! any new dependencies.
!-
      write(iout,'(A1)') '#' 
      write(iout,'(A27,A21,5x,I4)') '# Number of files given for',
     &                          ' dependency analysis: ', last_arg
      write(iout,'(A1)') '#'
      write(iout,'(A28,A8,17x,I4)') '# Number of newly discovered',
     &                                           ' files: ',n_disc
      write(iout,'(A1)') '#'
      write(iout,'(A26,A27,   I4)') '# Number of files excluded',
     &                         ' from dependency analysis: ',n_excl
      write(iout,'(A1)') '#'
      write(iout,'(A23,A28,2x,I4)') '# Total number of files',
     &                         ' analyzed for dependencies: ', nmax
      write(iout,'(A1)') '#'
      write(iout,'(A28,A15,8x,I6)') '# Total number of code lines',
     &                                   ' in all files: ', l_count
      write(iout,'(A1)') '#'

!-
      do n=1,nmax
        if (nlength(n).gt.0) then
         write(iout,'(A1)') '#'
         k=0
  11      continue
           do i=1,nlength(n)
             symbol(i)=name(i,n)
           enddo
           j=nlength(n)
           if (symbol(j).eq.'F') symbol(j)='o'
           j=j+1
           symbol(j)=':'
  12       j=j+1
            symbol(j)=' '
!-
             if (k.eq.0) then
              k=1
               if (name(nlength(n),n).eq.'F') then
                 do i=j+1,j+nlength(n)
                   symbol(i)=name(i-j,n)
                 enddo
                 j=j+nlength(n)+1
                 symbol(j)=' '
               endif
             endif
!-
  13         continue
             if (matrix(k,n)) then
               if (j+nlength(k).lt.length) then
                 do i=j+1,j+nlength(k)
                   symbol(i)=name(i-j,k)
                 enddo
                 j=j+nlength(k)
                 k=k+1
                 goto 12
               else
                 write(iout,'(80A1)') (symbol(i),i=1,j)
                 k=k+1
                 goto 11
               endif
             elseif (k.lt.nmax) then
               k=k+1
               goto 13
             else
               write(iout,'(80A1)') (symbol(i),i=1,j)
               if (j.gt.nlength(n)+2) then
                 write(iout,'(80A1)') (symbol(i),i=1,nlength(n)+1)
               endif
             endif
        endif
      enddo 
      close (iout)
      stop
      end


