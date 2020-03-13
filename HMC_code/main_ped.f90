program main
!	#30/07/2019	#Starting developing in my vacation!!
!	#01/08/2019	#finished HMC & GS
!	#02/08/2019	#finished RMHMC
!	#21/08/2019	#Pedigree information for EAAP 2019

	use kinds
	use sparsem

	implicit none

	integer::i,j,k,l,m
	integer::io
	integer::iwk

	integer::n_file
	integer::n_SNP
	integer::n_efc
	integer::n_f_efc
	integer::n_r_efc
	integer::n_tmp

	logical::found
	character::cwk

	integer::pos_trait
	integer,allocatable::pos_efc(:)
	character(len=10),allocatable::fac_efc(:),dis_efc(:),met_efc(:)
	character(len=20)::data_file,genome_file,pedigree_file


	real(kind=8),allocatable::dat(:,:)

	integer::n_g_file
	integer,allocatable::g_ID(:)
	real(kind=8),allocatable::g_dat(:,:),g_tmp(:),g_dat_adj(:,:),d_dat_adj(:,:)

	integer,allocatable::g_pos(:),p_pos(:)

	real(kind=8),allocatable::freq(:)
	real(kind=8)::var_freq_add,var_freq_dom,var_NOIA
	real(kind=8),allocatable::g_mat(:,:),d_mat(:,:),aa_mat(:,:),ad_mat(:,:),dd_mat(:,:)
	real(kind=8),allocatable::g_inv(:,:),d_inv(:,:),aa_inv(:,:),ad_inv(:,:),dd_inv(:,:)
	logical::add_flag,dom_flag,aa_flag,ad_flag,dd_flag,ped_flag

!	pedigree
	integer::n_p_file
	integer,allocatable::ped_dat(:,:)
	real(kind=8),allocatable::f(:),diagAinv(:)
	type(sparse_hashm):: Ainv_hash
	type(sparse_ija):: Ainv_ija


	integer,allocatable::size_levels(:)
	integer::t_n_fac,n_d,n_rnd

	character(len=10),allocatable::method_e(:),method_g(:)

	integer::LL
	integer::pos1,pos2
	real(kind=8)::mm,vv
	real(kind=8),allocatable::x(:,:),z(:),y(:),e(:),efc(:),xx(:)
	real(kind=8)::rhs,lhs,lambda

	real(kind=8),allocatable::var(:)
	real(kind=8)::var_r,quadratic,val

	integer::n_iter
	integer,allocatable::acc_efc(:),acc_var(:)

	real(kind=8),allocatable::tau(:)
	real(kind=8),allocatable::v_r(:),S2(:)
	real(kind=8),allocatable::v(:)

	real(kind=8)::random_gamma

!	Lapack inverse
	integer :: info
	integer :: LWORK
	integer, allocatable :: IPIV(:) ! dimension N
	real(kind=8), allocatable :: WORK(:) ! dimension LWORK

!	Time routine
	real(kind=8)::time_begin,time_end

!	Reading parameter file
	open(111,file="parameter_file")
	read(111,*,iostat=io)data_file
	if(io/=0)then
		print*,"Error data file name"
		stop
	endif
	read(111,*,iostat=io)pedigree_file
	if(io/=0)then
		print*,"Error pedigree file name"
		stop
	endif
	read(111,*,iostat=io)genome_file
	if(io/=0)then
		print*,"Error genome file name"
		stop
	endif
	read(111,*,iostat=io)n_SNP
	if(io/=0)then
		print*,"Error number of SNP"
		stop
	endif
	read(111,*,iostat=io)pos_trait
	read(111,*,iostat=io)n_efc
	allocate(pos_efc(n_efc),fac_efc(n_efc),dis_efc(n_efc))
	allocate(method_e(n_efc),method_g(n_efc+1),tau(n_efc),v_r(n_efc+1),S2(n_efc+1))
	n_rnd=0
	add_flag=.false.; dom_flag=.false.
	aa_flag=.false.; ad_flag=.false.; dd_flag=.false.
	ped_flag=.false.
	do i=1,n_efc
		read(111,*,iostat=io)pos_efc(i),fac_efc(i),dis_efc(i),method_e(i)
		if(io/=0)then
			print*,"Error ",i," effect"
			stop
		endif
		if(fac_efc(i)=="main")then
			read(111,*,iostat=io)tau(i)
			if(io/=0)then
				print*,"Error ",i," effect"
				stop
			endif
		elseif(fac_efc(i)=="cov")then
			dis_efc(i)="fixed"
			read(111,*,iostat=io)tau(i)
			if(io/=0)then
				print*,"Error ",i," effect"
				stop
			endif
		elseif(fac_efc(i)=="add")then
			dis_efc(i)="random"
			add_flag=.true.
			read(111,*,iostat=io)v_r(i),S2(i),method_g(i)
			if(io/=0)then
				print*,"Error ",i," effect"
				stop
			endif
		elseif(fac_efc(i)=="dom")then
			dis_efc(i)="random"
			dom_flag=.true.
			read(111,*,iostat=io)v_r(i),S2(i),method_g(i)
			if(io/=0)then
				print*,"Error ",i," effect"
				stop
			endif
		elseif(fac_efc(i)=="add_add")then
			dis_efc(i)="random"
			aa_flag=.true.
			read(111,*,iostat=io)v_r(i),S2(i),method_g(i)
			if(io/=0)then
				print*,"Error ",i," effect"
				stop
			endif
		elseif(fac_efc(i)=="add_dom")then
			dis_efc(i)="random"
			ad_flag=.true.
			read(111,*,iostat=io)v_r(i),S2(i),method_g(i)
			if(io/=0)then
				print*,"Error ",i," effect"
				stop
			endif
		elseif(fac_efc(i)=="dom_dom")then
			dis_efc(i)="random"
			dd_flag=.true.
			read(111,*,iostat=io)v_r(i),S2(i),method_g(i)
			if(io/=0)then
				print*,"Error ",i," effect"
				stop
			endif
		elseif(fac_efc(i)=="ped")then
			dis_efc(i)="random"
			ped_flag=.true.
			read(111,*,iostat=io)v_r(i),S2(i),method_g(i)
			if(io/=0)then
				print*,"Error ",i," effect"
				stop
			endif
		else
			print*,"No effect type ",i
			stop
		endif

		if(dis_efc(i)=="random")n_rnd=n_rnd+1
	enddo

	read(111,*,iostat=io)v_r(n_efc+1),S2(n_efc+1),method_g(n_efc+1)
	if(io/=0)then
		print*,"Error residual variance"
		stop
	endif
	read(111,*,iostat=io)n_iter
	if(io/=0)then
		print*,"Error no iteration"
		stop
	endif

	close(111)

	n_tmp = max(pos_trait,maxval(pos_efc))

!	Reading data file
	open(222,file=data_file)
	n_file=0
	do
		read(222,*,iostat=io)cwk
		if(io/=0)exit
		n_file=n_file+1
	enddo

	rewind(222)
	allocate(dat(n_file,n_tmp))
	do i=1,n_file
		read(222,*)dat(i,1:n_tmp)
	enddo
	close(222)

!	Reading pedigree file
	if(ped_flag)then
		open(444,file=pedigree_file)
		n_p_file=0
		do
			read(444,*,iostat=io)cwk
			if(io/=0)exit
			n_p_file=n_p_file+1
		enddo
		rewind(444)

		allocate(ped_dat(3,n_p_file),f(0:n_p_file),diagAinv(n_p_file))
		diagAinv=0.d0

		do i=1,n_p_file
			read(444,*)ped_dat(1:3,i)
		enddo
		close(444)

		call inbreeding_coefficient(n_p_file,ped_dat(2:3,:),f)
		call add_g_add(n_p_file,ped_dat(2:3,:),f,Ainv_ija,diagAinv)

		deallocate(ped_dat,f)
		n_d=n_file
	endif

!	Reading genome file
	if(add_flag .or. dom_flag)then
		open(333,file=genome_file)
		n_g_file=0
		do
			read(333,*,iostat=io)cwk
			if(io/=0)exit
			n_g_file=n_g_file+1
		enddo
		rewind(333)

		allocate(g_ID(n_g_file),g_dat(n_g_file,n_SNP))
		do i=1,n_g_file
			read(333,*)g_ID(i),g_dat(i,1:n_SNP)
		enddo
		close(333)

!	match g_ID with ID
		do i=1,n_efc
			if(fac_efc(i)=="add" .or. fac_efc(i)=="dom" .or. fac_efc(i)=="add_add" .or. &
				fac_efc(i)=="add_dom" .or. fac_efc(i)=="dom_dom")then
				iwk=pos_efc(i)
				exit
			endif
		enddo

		allocate(g_pos(n_file),p_pos(n_file))
		g_pos=0; n_d=0
		do i=1,n_file
			found=.false.
			do j=1,n_g_file
				if(dat(i,iwk)==g_ID(j))then
					n_d=n_d+1
					g_pos(n_d)=j
					p_pos(n_d)=i
					found=.true.
					exit
				endif
			enddo
		enddo
		n_d=n_file
		allocate(freq(n_SNP))
		freq=0.d0; var_freq_add=0.d0; var_freq_dom=0.d0
		do i=1,n_SNP
			do j=1,n_d
				freq(i)=freq(i)+g_dat(g_pos(j),i)
			enddo
			freq(i)=freq(i)/(2*n_d)
			var_freq_add=var_freq_add+2*freq(i)*(1-freq(i))
			var_freq_dom=var_freq_dom+4*freq(i)**2*(1-freq(i))**2
		enddo
!	Construct genomic relationship matrix
		if(add_flag)then
			allocate(g_mat(n_d,n_d),g_dat_adj(n_d,n_SNP))
			g_mat=0.d0
			do i=1,n_SNP
				do j=1,n_d
					g_dat_adj(j,i)=int(g_dat(g_pos(j),i)+0.5)-2*freq(i)
				enddo
			enddo
			call dgemm("N","T",n_d,n_d,n_SNP,1.d0,g_dat_adj,n_d,g_dat_adj,n_d,0.d0,g_mat,n_d)
			g_mat=g_mat/var_freq_add
			do i=1,n_d
				g_mat(i,i)=g_mat(i,i)+0.01
			enddo

			LWORK=64*n_d
			allocate(g_inv(n_d,n_d),IPIV(n_d),WORK(LWORK))
			g_inv=g_mat
		    call DGETRF(n_d, n_d, g_inv, n_d, IPIV(1:n_d), info)  ! factorize
		    call DGETRI(n_d, g_inv, n_d, IPIV, WORK, LWORK, info)  ! inverse
			deallocate(g_dat_adj,IPIV,WORK)
		endif

!	Construct dominace relationship matrix
		if(dom_flag)then
			allocate(d_mat(n_d,n_d),d_dat_adj(n_d,n_SNP))
			d_mat=0.d0
			do i=1,n_SNP
				do j=1,n_d
					if(int(g_dat(g_pos(j),i)+0.5)==0)then
						d_dat_adj(j,i)=-2*freq(i)**2
					elseif(int(g_dat(g_pos(j),i)+0.5)==1)then
						d_dat_adj(j,i)=2*freq(i)*(1-freq(i))
					elseif(int(g_dat(g_pos(j),i)+0.5)==2)then
						d_dat_adj(j,i)=-2*(1-freq(i))**2
					else
						print*,i,j,int(g_dat(g_pos(j),i)+0.5)
						stop
					endif
				enddo
			enddo
			call dgemm("N","T",n_d,n_d,n_SNP,1.d0,d_dat_adj,n_d,d_dat_adj,n_d,0.d0,d_mat,n_d)
			d_mat=d_mat/var_freq_dom

			do i=1,n_d
				d_mat(i,i)=d_mat(i,i)+0.001
			enddo

			LWORK=64*n_d
			allocate(d_inv(n_d,n_d),IPIV(n_d),WORK(LWORK))
			d_inv=d_mat
		    call DGETRF(n_d, n_d, d_inv, n_d, IPIV(1:n_d), info)  ! factorize
		    call DGETRI(n_d, d_inv, n_d, IPIV, WORK, LWORK, info)  ! inverse
			deallocate(d_dat_adj,IPIV,WORK)

		endif

!	Construct add_add relationship matrix
		if(aa_flag)then
			allocate(aa_mat(n_d,n_d))

			aa_mat=g_mat*g_mat
			var_NOIA=0.d0
			do i=1,n_d
				var_NOIA=var_NOIA+aa_mat(i,i)
			enddo
			var_NOIA=var_NOIA/n_d
			aa_mat=aa_mat/var_NOIA

			do i=1,n_d
				aa_mat(i,i)=aa_mat(i,i)+0.001
			enddo

			LWORK=64*n_d
			allocate(aa_inv(n_d,n_d),IPIV(n_d),WORK(LWORK))
			aa_inv=aa_mat
		    call DGETRF(n_d, n_d, aa_inv, n_d, IPIV(1:n_d), info)  ! factorize
		    call DGETRI(n_d, aa_inv, n_d, IPIV, WORK, LWORK, info)  ! inverse
			deallocate(aa_mat,IPIV,WORK)
		endif

!	Construct add_dom relationship matrix

		if(ad_flag)then
			allocate(ad_mat(n_d,n_d))
			ad_mat=g_mat*d_mat
			var_NOIA=0.d0
			do i=1,n_d
				var_NOIA=var_NOIA+ad_mat(i,i)
			enddo
			var_NOIA=var_NOIA/n_d
			ad_mat=ad_mat/var_NOIA

			do i=1,n_d
				ad_mat(i,i)=ad_mat(i,i)+0.001
			enddo

			LWORK=64*n_d
			allocate(ad_inv(n_d,n_d),IPIV(n_d),WORK(LWORK))
			ad_inv=ad_mat
		    call DGETRF(n_d, n_d, ad_inv, n_d, IPIV(1:n_d), info)  ! factorize
		    call DGETRI(n_d, ad_inv, n_d, IPIV, WORK, LWORK, info)  ! inverse
			deallocate(IPIV,WORK)

		endif

!	Construct add_dom relationship matrix

		if(dd_flag)then
			allocate(dd_mat(n_d,n_d))
			dd_mat=d_mat*d_mat
			var_NOIA=0.d0
			do i=1,n_d
				var_NOIA=var_NOIA+dd_mat(i,i)
			enddo
			var_NOIA=var_NOIA/n_d
			dd_mat=dd_mat/var_NOIA

			do i=1,n_d
				dd_mat(i,i)=dd_mat(i,i)+0.001
			enddo

			LWORK=64*n_d
			allocate(dd_inv(n_d,n_d),IPIV(n_d),WORK(LWORK))
			dd_inv=dd_mat
		    call DGETRF(n_d, n_d,dd_inv, n_d, IPIV(1:n_d), info)  ! factorize
		    call DGETRI(n_d, dd_inv, n_d, IPIV, WORK, LWORK, info)  ! inverse
			deallocate(IPIV,WORK)
		endif

		if(allocated(g_mat))deallocate(g_mat)
		if(allocated(d_mat))deallocate(d_mat)
		if(allocated(aa_mat))deallocate(aa_mat)
		if(allocated(ad_mat))deallocate(ad_mat)
		if(allocated(dd_mat))deallocate(dd_mat)
	endif

!	Preparing McMC
	allocate(size_levels(n_efc),x(n_d,n_efc),y(n_d),e(n_d))
	e=0.d0

	do i=1,n_efc
		if(fac_efc(i)=="cov")then
			mm=0.d0
			do j=1,n_d
				x(j,i)=dat(j,pos_efc(i))
				mm=mm+dat(j,pos_efc(i))
			enddo
			mm=mm/n_d
			x(:,i)=x(:,i)-mm
			size_levels(i)=1
		elseif(fac_efc(i)=="main")then
			iwk=0
			do j=1,n_d
				x(j,i)=dat(j,pos_efc(i))
				if(x(j,i)>iwk)iwk=x(j,i)
			enddo
			size_levels(i)=iwk
		elseif(fac_efc(i)=="add" .or. fac_efc(i)=="dom" .or. fac_efc(i)=="add_add" .or. &
				fac_efc(i)=="add_dom" .or. fac_efc(i)=="dom_dom")then
			do j=1,n_d
				x(j,i)=j
			enddo
			size_levels(i)=n_d

		elseif(fac_efc(i)=="ped")then
			allocate(z(n_p_file))
			z=0
			do j=1,n_d
				z(int(dat(j,pos_efc(i))))=j
				x(j,i)=dat(j,pos_efc(i))
			enddo
			size_levels(i)=n_p_file
		endif
	enddo
	mm=0.d0; vv=0.d0
	do i=1,n_d
		y(i)=dat(i,pos_trait)
		e(i)=dat(i,pos_trait)
		mm=mm+e(i)
		vv=vv+e(i)**2
	enddo
	vv=vv/n_d-(mm/n_d)**2
	deallocate(dat)

	t_n_fac=sum(size_levels)
	allocate(efc(t_n_fac),xx(t_n_fac),var(n_efc),acc_efc(t_n_fac),acc_var(n_efc+1))
	efc=0.d0; xx=0.d0; var=0.d0; acc_efc=0; acc_var=0
	iwk=0; var_r=vv/2
	do i=1,n_efc
		if(fac_efc(i)=="cov")then
			do j=1,n_d
				xx(iwk+1)=xx(iwk+1)+x(j,i)**2
			enddo
		elseif(fac_efc(i)=="main")then
			do j=1,size_levels(i)
				do k=1,n_d
					if(x(k,i)==j) xx(iwk+j)=xx(iwk+j)+1.d0
				enddo
			enddo
		elseif(fac_efc(i)=="add" .or. fac_efc(i)=="dom" .or. fac_efc(i)=="add_add" .or. &
				fac_efc(i)=="add_dom" .or. fac_efc(i)=="dom_dom")then
			do j=1,n_d
				xx(iwk+j)=1.d0
			enddo

		elseif(fac_efc(i)=="ped" )then
			do j=1,n_d
				xx(iwk+int(x(j,i)))=1.d0
			enddo
		endif
		iwk=iwk+size_levels(i)
		if(dis_efc(i)=="random")then
			var(i)=(vv/2)/n_rnd
		else
			var(i)=tau(i)
		endif
	enddo

	open(7777,file="sample_effects")
	open(8888,file="sample_variances")
	call cpu_time(time_begin)
!	Sampling Start
	do LL=1,n_iter
		iwk=0
		do i=1,n_efc
			do j=1,size_levels(i)
				rhs=0.d0; lhs=0.d0
				if(fac_efc(i)=="cov")then
					rhs=rhs+dot_product(x(:,i),e)
				elseif(fac_efc(i)=="main")then
					do k=1,n_d
						if(x(k,i)==j)rhs=rhs+e(k)
					enddo
				elseif(fac_efc(i)=="add" .or. fac_efc(i)=="dom" .or. fac_efc(i)=="add_add" .or. &
					fac_efc(i)=="add_dom" .or. fac_efc(i)=="dom_dom")then
					rhs=rhs+e(int(x(j,i)))
				elseif(fac_efc(i)=="ped")then
					if(z(j)/=0) rhs=rhs+e(int(z(j)))
				endif
				rhs=rhs+xx(iwk+j)*efc(iwk+j)
				if(dis_efc(i)=="random")then
					if(fac_efc(i)=="add")then
						pos1=iwk+1; pos2=iwk+size_levels(i)
						lambda=var_r/var(i)
						rhs=rhs+g_inv(j,j)*efc(iwk+j)*lambda
						rhs=rhs-dot_product(g_inv(:,j),efc(pos1:pos2))*lambda
						lhs=xx(iwk+j)+g_inv(j,j)*lambda
					elseif(fac_efc(i)=="dom")then
						pos1=iwk+1; pos2=iwk+size_levels(i)
						lambda=var_r/var(i)
						rhs=rhs+d_inv(j,j)*efc(iwk+j)*lambda
						rhs=rhs-dot_product(d_inv(:,j),efc(pos1:pos2))*lambda
						lhs=xx(iwk+j)+d_inv(j,j)*lambda
					elseif(fac_efc(i)=="add_add")then
						pos1=iwk+1; pos2=iwk+size_levels(i)
						lambda=var_r/var(i)
						rhs=rhs+aa_inv(j,j)*efc(iwk+j)*lambda
						rhs=rhs-dot_product(aa_inv(:,j),efc(pos1:pos2))*lambda
						lhs=xx(iwk+j)+aa_inv(j,j)*lambda
					elseif(fac_efc(i)=="add_dom")then
						pos1=iwk+1; pos2=iwk+size_levels(i)
						lambda=var_r/var(i)
						rhs=rhs+ad_inv(j,j)*efc(iwk+j)*lambda
						rhs=rhs-dot_product(ad_inv(:,j),efc(pos1:pos2))*lambda
						lhs=xx(iwk+j)+ad_inv(j,j)*lambda
					elseif(fac_efc(i)=="dom_dom")then
						pos1=iwk+1; pos2=iwk+size_levels(i)
						lambda=var_r/var(i)
						rhs=rhs+dd_inv(j,j)*efc(iwk+j)*lambda
						rhs=rhs-dot_product(dd_inv(:,j),efc(pos1:pos2))*lambda
						lhs=xx(iwk+j)+dd_inv(j,j)*lambda
					elseif(fac_efc(i)=="ped")then
						pos1=iwk
						lambda=var_r/var(i)
						rhs=rhs+diagAinv(j)*efc(pos1+j)*lambda
						do k=Ainv_ija%ia(j),Ainv_ija%ia(j+1)-1
							rhs=rhs-Ainv_ija%a(k)*efc(pos1+Ainv_ija%ja(k))*lambda
						enddo
						lhs=xx(iwk+j)+diagAinv(j)*lambda
					else
						lambda=var_r/var(i)
						lhs=xx(iwk+j)+lambda
					endif
				else
					lambda=var_r/var(i)
					lhs=xx(iwk+j)+lambda
				endif

				call sample_normal(rhs/var_r,lhs/var_r,efc(iwk+j),val,acc_efc(iwk+j),method_e(i))
				if(fac_efc(i)=="cov")then
					e=e-x(:,i)*(val-efc(iwk+j))
				elseif(fac_efc(i)=="main")then
					do k=1,n_d
						if(x(k,i)==j)e(k)=e(k)-(val-efc(iwk+j))
					enddo
				elseif(fac_efc(i)=="add" .or. fac_efc(i)=="dom" .or. fac_efc(i)=="add_add" .or. &
						fac_efc(i)=="add_dom" .or. fac_efc(i)=="dom_dom")then
					e(int(x(j,i)))=e(int(x(j,i)))-(val-efc(iwk+j))
				elseif(fac_efc(i)=="ped")then
					if(z(j)/=0) e(int(z(j)))=e(int(z(j)))-(val-efc(iwk+j))
				endif
				efc(iwk+j)=val
			enddo

			if(dis_efc(i)=="random")then
				quadratic=v_r(i)*S2(i)
				if(fac_efc(i)=="main")then
					pos1=iwk+1; pos2=iwk+size_levels(i)
					quadratic=dot_product(efc(pos1:pos2),efc(pos1:pos2))
				elseif(fac_efc(i)=="add")then
					pos1=iwk+1; pos2=iwk+size_levels(i)
					quadratic=dot_product(matmul(g_inv,efc(pos1:pos2)),efc(pos1:pos2))
				elseif(fac_efc(i)=="dom")then
					pos1=iwk+1; pos2=iwk+size_levels(i)
					quadratic=dot_product(matmul(d_inv,efc(pos1:pos2)),efc(pos1:pos2))
				elseif(fac_efc(i)=="add_add")then
					pos1=iwk+1; pos2=iwk+size_levels(i)
					quadratic=dot_product(matmul(aa_inv,efc(pos1:pos2)),efc(pos1:pos2))
				elseif(fac_efc(i)=="add_dom")then
					pos1=iwk+1; pos2=iwk+size_levels(i)
					quadratic=dot_product(matmul(ad_inv,efc(pos1:pos2)),efc(pos1:pos2))
				elseif(fac_efc(i)=="dom_dom")then
					pos1=iwk+1; pos2=iwk+size_levels(i)
					quadratic=dot_product(matmul(dd_inv,efc(pos1:pos2)),efc(pos1:pos2))
				elseif(fac_efc(i)=="ped")then
					pos1=iwk
					do j=1,size_levels(i)
						do k=Ainv_ija%ia(j),Ainv_ija%ia(j+1)-1
							quadratic=quadratic+efc(pos1+j)*Ainv_ija%a(k)*efc(pos1+Ainv_ija%ja(k))
						enddo
					enddo
				endif

!	Sampling from inverse gamma
				call sample_i_gamma(real(size_levels(i)+v_r(i),kind=8),quadratic,var(i),val,acc_var(i),method_g(i))
				var(i)=val
			endif
			iwk=iwk+size_levels(i)
		enddo

!	Sampling from inverse gamma
		quadratic=dot_product(e,e)+v_r(n_efc+1)*S2(n_efc+1)
		call sample_i_gamma(real(n_d+v_r(n_efc+1),kind=8),quadratic,var_r,val,acc_var(n_efc+1),method_g(n_efc+1))
		var_r=val

		print*,LL,var,var_r
		write(7777,'(10000F15.10)')efc(1:t_n_fac)
		write(8888,'(100F15.8)')var,var_r
	enddo
	call cpu_time(time_end)
	print *, "cpu time:", time_end-time_begin, "seconds."
	open(1111,file="computing_time")
	write(1111,*)"cpu time:", time_end-time_begin, "seconds."
	close(1111)
	open(2222,file="acceptance_ratio_effects")
	do i=1,t_n_fac
		write(2222,*)acc_efc(i)
	enddo
	close(2222)
	open(3333,file="acceptance_ratio_variances")
	do i=1,n_efc+1
		write(3333,*)acc_var(i)
	enddo
	close(3333)

	end


subroutine sample_normal(rhs,lhs,val,val_new,acc,method)

	implicit none

	real(kind=8),intent(in)::rhs,lhs
	real(kind=8)::val,val_new
	real(kind=8)::r_normal

	integer,intent(inout)::acc

	integer::LL,leap
	real(kind=8)::epsilon,p
	real(kind=8)::U_before,U_after
	real(kind=8)::K_before,K_after
	real(kind=8)::H_before,H_after

!	RMHMC
	integer::fix_loop
	real(kind=8)::G

	real(kind=8)::unif

	character(len=10),intent(in)::method

	if(method=="gibbs")then
		acc=acc+1
		val_new = rhs/lhs+sqrt(1.d0/lhs)*dble(r_normal())

	elseif(method=="HMC")then
		val_new=val
		LL=7
		epsilon=sqrt(1.d0/lhs)/(0.1589825*20)

		p=r_normal()
		U_before=-(rhs*val_new-0.5*lhs*val_new**2)
		K_before=p**2/2
		H_before=(U_before+K_before)

		do leap=1,LL
			p=p+0.5*epsilon*(rhs-lhs*val_new)
			val_new=val_new+epsilon*p
			p=p+0.5*epsilon*(rhs-lhs*val_new)
		enddo

		U_after=-(rhs*val_new-0.5*lhs*val_new**2)
		K_after=p**2/2
		H_after=(U_after+K_after)

		call random_number(unif)

		if(unif<exp(H_before-H_after))then
			acc=acc+1
		else
			val_new=val
		endif


	elseif(method=="RMHMC")then
		val_new=val
		LL=20
		epsilon=0.1

		G=lhs
		p=sqrt(G)*r_normal()
		U_before=-(rhs*val_new-0.5*lhs*val_new**2)
		K_before=p**2/(2*G)
		H_before=(U_before+K_before)

		do leap=1,LL
			p=p+0.5*epsilon*(rhs-lhs*val_new)
			val_new=val_new+epsilon*p/G
			p=p+0.5*epsilon*(rhs-lhs*val_new)
		enddo

		U_after=-(rhs*val_new-0.5*lhs*val_new**2)
		K_after=p**2/(2*G)
		H_after=(U_after+K_after)
		call random_number(unif)

		if(unif<exp(H_before-H_after))then
			acc=acc+1
		else
			val_new=val
		endif
	endif

	return
end subroutine


subroutine sample_i_gamma(nb,xx,val,val_new,acc,method)

	implicit none

	real(kind=8),intent(in)::nb,xx
	real(kind=8)::val,val_new

	real(kind=8)::r_normal
	real(kind=8)::random_gamma

	integer,intent(inout)::acc

	integer::LL,leap
	real(kind=8)::epsilon,p
	real(kind=8)::U_before,U_after
	real(kind=8)::K_before,K_after
	real(kind=8)::H_before,H_after

!	RMHMC
	integer::fix_loop
	real(kind=8)::G,val_new1,p1
	real(kind=8)::dLdx,dGdX,dHdp,dHdx
	real(kind=8)::dHdp1

	real(kind=8)::unif

	character(len=10),intent(in)::method

	if(method=="gibbs")then
		acc=acc+1
		val_new = 1.d0/random_gamma(nb,1/xx,.true.)

	elseif(method=="HMC")then
		val_new=val
		LL=7
		epsilon=sqrt(xx**2/((nb-1)**2*(nb-2)))/(0.112485939*20)

		p=r_normal()
		U_before=-(-(2+nb)*log(val_new)-xx/val_new)/2.d0
		K_before=p**2/2
		H_before=(U_before+K_before)

		do leap=1,LL
			p=p+0.5*epsilon*(-(2.d0+nb)/(2.d0*val_new)+xx/(2.d0*val_new**2))
			val_new=val_new+epsilon*p
			p=p+0.5*epsilon*(-(2.d0+nb)/(2.d0*val_new)+xx/(2.d0*val_new**2))
		enddo

		U_after=-(-(2+nb)*log(val_new)-xx/val_new)/2.d0
		K_after=p**2/2
		H_after=(U_after+K_after)

		call random_number(unif)

		if(unif<exp(H_before-H_after))then
			acc=acc+1
		else
			val_new=val
		endif

	elseif(method=="RMHMC")then
		val_new=val
		LL=20
		epsilon=0.1

		G=-0.5*(-nb+2)/(val_new**2)
		p=sqrt(G)*r_normal()
		U_before=-(-(2+nb)*log(val_new)-xx/val_new)/2.d0
		K_before=p**2/(2*G)
		H_before=(U_before+K_before)+log(G)/2

		dLdx=-(2.d0+nb)/(2.d0*val_new)+xx/(2.d0*val_new**2)
		dGdx=(-nb+2)/(val_new**3)
		dHdp=p/G
		dHdx=0.5*(dGdx/G-dHdp*dGdx*dHdp)-dLdx

		do leap=1,LL
			do fix_loop=1,10
				p1=p-0.5*epsilon*dHdx
				dHdp=p1/G
				dHdx=0.5*(dGdx/G-dHdp*dGdx*dHdp)-dLdx
			enddo
			p=p1

			dHdp=p/G
			dHdp1=dHdp

			do fix_loop=1,10
				val_new1=val_new+0.5*epsilon*(dHdp+dHdp1)
				G=-0.5*(-nb+2)/(val_new**2)
				dHdp1=p/G
			enddo

			val_new=val_new1
			dLdx=-(2.d0+nb)/(2.d0*val_new)+xx/(2.d0*val_new**2)
			dGdx=(-nb+2)/(val_new**3)
			dHdp=p/G
			dHdx=0.5*(dGdx/G-dHdp*dGdx*dHdp)-dLdx
			p=p-0.5*epsilon*dHdx

		enddo

		U_after=-(-(2+nb)*log(val_new)-xx/val_new)/2.d0
		K_after=p**2/(2*G)
		H_after=(U_after+K_after)+log(G)/2
		call random_number(unif)

		if(unif<exp(H_before-H_after))then
			acc=acc+1
		else
			val_new=val
		endif
	endif

	return
end subroutine
