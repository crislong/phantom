!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module sort_particles
!
! sorts the particles so neighbours are also close in memory
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: io, linklist, part, sortutils
!
 implicit none
 public :: sort_part_radius, sort_part_id, sort_part

 private

contains

! do not even compile with this routine if not sorting
!----------------------------------------------------------------
!+
!  this version sorts the particles using the neighbour lists
!+
!----------------------------------------------------------------
subroutine sort_part
 use io,       only:iprint,fatal
 use part,     only:reorder_particles,npart,ll,xyzh,vxyzu,isdead
 use linklist, only:set_linklist,ncells,ifirstincell
 integer         :: i,ipart,iprev,ifirst
 integer(kind=8) :: icell
 real            :: t0,t1,t2

 call cpu_time(t0)
 call set_linklist(npart,npart,xyzh,vxyzu)  ! don't include MPI ghosts
 call cpu_time(t1)
 write(iprint,*) '> sorting particles...',t1-t0,'s for linklist'

 ipart = 0
 iprev = 0
 ifirst = ifirstincell(1)

 do icell=1,ncells
!    call get_neighbour_list(icell,listneigh,nneigh,xyzh,xyzcache,0)

    i = ifirstincell(icell)
    if (ifirst==0) ifirst = i
    !
    !--link end of last cell to start of this cell
    !
    if (iprev /= 0) ll(iprev) = i
    !
    !--loop to the end of this cell
    !
    do while (i /= 0 .and. ipart < npart)
       ipart = ipart + 1
       !print*,'ipart = ',ipart,'iprev = ',iprev
       if (i <= 0) call fatal('sort','internal error: i<=0',i)
       iprev = i
       i = ll(i)
    enddo
 enddo

 print*,'finished',ipart,npart
 print*,'setting',iprev,ifirst
!
!--make sure we haven't missed any particles
!
 if (ipart /= npart) then
    write(iprint,*) ' shifting ',npart-ipart,' particles not in neighbour list to end of arrays'
    do i=1,npart
       if (isdead(i)) then
          ipart = ipart + 1
          ll(iprev) = i
          iprev = i
       endif
    enddo
 endif
 if (ipart /= npart) call fatal('sort','ipart /= npart')
!
!--complete the loop by linking the last particle of the last cell
!  to the first of the first cell
!
 ll(iprev) = ifirst

 if (any(ll(1:npart) <= 0)) then
    ipart = 0
    do i=1,npart
       if (ll(i)==0) ipart = ipart + 1
    enddo
    print*,'number of errors = ',ipart
    call fatal('sort','errors in index array')
 endif

! print*,' ncells = ',ncells,' ipart = ',ipart,' npart = ',npart
 write(iprint,*) ' copying arrays...'
!
!--copy arrays into correct order
!
 call reorder_particles(ll,npart)

 call cpu_time(t2)
 write(iprint,*) '> sort completed in ',t2-t1,'s'

end subroutine sort_part

!----------------------------------------------------------------
!+
!  this version sorts the particles by radius (e.g. for a disc)
!+
!----------------------------------------------------------------
subroutine sort_part_radius(np)
 use io,        only:iprint,error
 use part,      only:xyzh,reorder_particles,npart,ll
 use sortutils, only:indexxfunc,r2func
 integer, intent(in) :: np
 real :: t1,t2

 call cpu_time(t1)
 write(iprint,*) '> sorting particles in radius...'
 if (np /= npart) call error('sort','np /= npart')

 call indexxfunc(npart,r2func,xyzh,ll)

 write(iprint,*) ' copying arrays...'
!
!--copy arrays into correct order
!
 call reorder_particles(ll,npart)

 call cpu_time(t2)
 write(iprint,*) '> sort completed in ',t2-t1,'s'

end subroutine sort_part_radius

!----------------------------------------------------------------
!+
!  this version sorts the particles by ID
!+
!----------------------------------------------------------------
subroutine sort_part_id
 use io,        only:iprint,error
 use part,      only:reorder_particles,npart,ll,iorig
 use sortutils, only:indexx
 real :: t1,t2

 call cpu_time(t1)
 write(iprint,*) '> sorting particles by ID...'

 call indexx(npart,iorig,ll)

 write(iprint,*) ' copying arrays...'
 !
 !--copy arrays into correct order
 !
 call reorder_particles(ll,npart)

 call cpu_time(t2)
 write(iprint,*) '> sort completed in ',t2-t1,'s'

end subroutine sort_part_id

end module sort_particles
