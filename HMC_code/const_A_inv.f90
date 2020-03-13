!=======================================================================
subroutine inbreeding_coefficient (n,ped,f)
!	THE Mewissen and Luo Z, 1992.
!	3/16/2011 developed by A. Arakawa
!	5/24/2011 defraged by A. Arakawa
!	This routine is Meuwissen and Ruo (1992) algorithm
!-----------------------------------------------------------------------
	implicit none

	integer,intent(in)::n
	integer::ped(2,n)
	real(kind=8)::f(0:n)

	integer::i,j,k
	integer,allocatable::point(:)
	integer::isire,idam,ksire,kdam
	real(kind=8),allocatable::l(:),d(:)
	real(kind=8)::fi,r
!-----------------------------------------------------------------------
	allocate(point(n),l(n),d(n))
	f=0.d0; l=0.d0; d=0.d0; point=0; 
	f(0)=-1.d0
	do i=1,n
		isire=ped(1,i); idam=ped(2,i)
		ped(1,i)=max(isire,idam); ped(2,i)=min(isire,idam)
		d(i)=0.5d0-0.25d0*(f(isire)+f(idam))
		if((isire==0).or.(idam==0))then
			f(i)=0.d0
		else
			fi=-1.d0; l(i)=1.d0; j=i
			do while(j/=0)
				k=j
				r=0.5d0*l(k)
				ksire=ped(1,k); kdam=ped(2,k)
				if(ksire>0)then
					do while(point(k) > ksire)
						k=point(k)
					enddo
					l(ksire)=l(ksire)+r
					if(ksire/=point(k))then
						point(ksire)=point(k)
						point(k)=ksire
					endif
					if(kdam>0)then
						do while(point(k)>kdam)
							k=point(k)
						enddo
						l(kdam)=l(kdam)+r
						if(kdam/=point(k)) then
							point(kdam)=point(k)
							point(k)=kdam
						endif
					endif
				endif
				fi=fi+l(j)*l(j)*d(j)
				l(j)=0.d0
				k=j
				j=point(j)
				point(k)=0
			enddo
			f(i)=fi
		endif
	enddo
	deallocate(point,l,d)

return
end subroutine inbreeding_coefficient
!=======================================================================



!=======================================================================
subroutine add_g_add(n,ped,f,Ainv_ija,diagAinv)
!	3/16/2011 developed by A. Arakawa
!	This subroutine is based on Dr Misztal's code
	use sparsem
	
!-----------------------------------------------------------------------
	implicit none
	integer,intent(in)::n
	real(kind=8),intent(in)::f(0:n)
	integer,intent(in)::ped(2,n)
	real(kind=8)::diagAinv(n)

	integer :: i,j,k,p(3)
	real(r8) ::t(3),mendel
	type(sparse_hashm):: Ainv_hash
	type(sparse_ija):: Ainv_ija
!-----------------------------------------------------------------------

	call init(Ainv_hash)
	call zerom(Ainv_hash,n)
	call init(Ainv_ija)
	call zerom(Ainv_ija,n)
  
	t=(/1., -.5, -.5/)

	do i=1,n
		p(1)=i
		p(2)=ped(1,i)
		p(3)=ped(2,i)
		mendel=1
		if(p(2)/=0) mendel=mendel-0.25-0.25*f(p(2))
		if(p(3)/=0) mendel=mendel-0.25-0.25*f(p(3))
     
		do j=1,3
			do k=1,3
				if (p(j)/= 0 .and. p(k)/=0) then
					if(p(j)==p(k)) diagAinv(p(j))=diagAinv(p(j))+t(j)*t(k)/mendel
					call addm(t(j)*t(k)/mendel,p(j),p(k),Ainv_hash)
				endif
			enddo
		enddo
	enddo
	Ainv_ija=ainv_hash
	call convert_hash_ija_general(Ainv_ija,Ainv_hash,conv_upper_full)  !convert to full-stored Ainv matrix
	call reset(Ainv_hash)
	
	return
	end subroutine
!=======================================================================
