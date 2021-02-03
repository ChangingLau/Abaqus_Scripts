! This is a fortran program for Newton Method
		PROGRAM newtonMethod
		implicit none 
		real,parameter :: e=1e-4
		integer :: i=1
		real :: xn,xnp1,f,fd,root,realE
		! Require input from user
		write(*,*)"Input the initiate x:"
		read(*,*)xn
		if (f(xn)==0) then
			root = xn
		else
			do 
				xnp1=xn-(f(xn)/fd(xn))
				realE=abs((xnp1-xn)/xnp1)
				write(*,*)"[",i,"], x_n=",xn,", x_n+1=",xnp1,", relaE=",realE
				if (realE<e) exit
				xn=xnp1
				i=i+1
			end do
			write(*,*)"root=",xn, ", iterations=",i
		end if
	END PROGRAM newtonMethod

! functions below should be modified for various expressions
	! initial function
	real function f(x)
	real :: x
	f=x**3-x+1
	end

	! function's derivative
	real function fd(x)
	real :: x
	fd=3*x**2-1
	end
