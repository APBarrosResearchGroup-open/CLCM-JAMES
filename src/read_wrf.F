module read_wrf
implicit none
public :: read2D, read2Dinteger, read3D, getbdytime
CONTAINS
subroutine read2D(value,varname,file_prefix,myi,myj,si,ei,sj,ej)
    integer, intent(in) :: si,ei,sj,ej
    real, dimension(si:ei,sj:ej),intent(out) :: value
    ! double precision, dimension(si:ei,sj:ej) :: value_temp
    character(10), intent(in) :: varname 
    character(30), intent(in) :: file_prefix
    character(100) :: filename
    character(20) :: myistr, myjstr
    integer :: unit, ios
    character(100) :: line
    integer :: myi, myj
    integer :: idim1, idim2
    integer :: ierr
    integer :: iunit = 240414
    integer :: ind_value = 0
    write(myistr,'(I4.4)') myi
    write(myjstr,'(I4.4)') myj
    write(filename,'(A)') trim(file_prefix)//trim(myistr)//"_"//trim(myjstr)
    ! if (myi.eq.1 .and. myj.eq.1) then
        ! print *, "reading ", varname, " from ",trim(filename)
    ! end if
    open(unit=iunit, file=trim(filename),status='old', action='read',iostat=ierr)
    if (ierr /= 0) then
        print *, "Error opening file ", trim(filename)
        call stopcm1
    end if
    do
        read(iunit, '(A)', iostat=ios) line
        if (ios /= 0) then
            print *, "ERROR: varname not found ", varname," in ",filename
            call stopcm1
        end if
        if (trim(line) .eq. trim(varname)) then

            do idim1 = si, ei
                do idim2 = sj, ej
                    ! read(iunit, '(A)', iostat=ios) line
                    ! if (trim(varname) .eq. "piBXS") then
                    !     print *, line
                    ind_value = ind_value + 1
                    read (iunit,*,iostat=ios) value(idim1,idim2)
                        ! print *, value(idim1,idim2)
                    ! end if
                    if (ios /= 0) then
                        print *, "ERROR (","ios=",ios, "): reading ", varname, " at ", idim1, idim2," from ", trim(filename), " ind_value=",ind_value
                        call stopcm1
                    end if
                end do

            end do
            exit
        end if
    end do
    ! if (myi.eq.1 .and. myj.eq.1) then
    !     print *, "Finish reading ", varname
    ! end if

    close(iunit)
end subroutine read2D

subroutine read3D(value,varname,file_prefix,myi,myj,si,ei,sj,ej,sk,ek)
    integer, intent(in) :: si,ei,sj,ej,sk,ek
    real, dimension(si:ei,sj:ej,sk:ek),intent(out) :: value
    character(10), intent(in) :: varname
    character(30), intent(in) :: file_prefix
    character(100) :: filename
    character(20) :: name, myistr, myjstr
    integer :: unit, ios
    character(100) :: line
    integer :: myi, myj
    integer :: idim1, idim2, idim3
    integer :: ierr
    integer :: iunit = 240414
    if (myi.ge.0 .and. myj.ge.0 ) then
        write(myistr,'(I4.4)') myi
        write(myjstr,'(I4.4)') myj
        write(filename,'(A)') trim(file_prefix)//trim(myistr)//"_"//trim(myjstr)
    else
        write(filename,'(A)') trim(file_prefix)
    end if

    if (myi.eq.1 .and. myj.eq.1) then
        print *, "reading ", varname, " from ",trim(filename)
    end if

    open(unit=iunit, file=trim(filename),status='old', action='read',iostat=ierr)
    if (ierr /= 0) then
        print *, "Error opening file ", trim(filename)
        call stopcm1
    end if
    ! print *, "reading "//trim(varname)//" from "//trim(filename)
    do
        read(iunit, '(A)', iostat=ios) line
        if (ios /= 0) then
            print *, "ERROR: varname not found ", varname
            call stopcm1
        end if
        if (trim(line) .eq. trim(varname)) then

            do idim1 = si, ei
                do idim2 = sj, ej
                    do idim3 = sk, ek
                        read (iunit,*,iostat=ios) value(idim1,idim2,idim3)
                        if (ios /= 0) then
                            print *, "ERROR: reading ", varname, " at ", idim1, idim2,idim3," from ", trim(filename)
                            call stopcm1
                        end if
                    end do
                end do

            end do
            exit
        end if
    end do

    close(iunit)
    if (myi.eq.1 .and. myj.eq.1) then
        print *, "Finish reading ", varname, " from ",trim(filename)
    end if
end subroutine read3D

subroutine read2Dinteger(value,varname,file_prefix,myi,myj,si,ei,sj,ej)
    integer, intent(in) :: si,ei,sj,ej
    integer, dimension(si:ei,sj:ej),intent(out) :: value
    character(10), intent(in) :: varname
    character(30), intent(in) :: file_prefix
    character(100) :: filename
    character(20) :: name, myistr, myjstr
    integer :: unit, ios
    character(100) :: line
    integer :: myi, myj
    integer :: idim1, idim2
    integer :: ierr
    integer :: iunit = 240414
    write(myistr,'(I4.4)') myi
    write(myjstr,'(I4.4)') myj
    write(filename,'(A)') trim(file_prefix)//trim(myistr)//"_"//trim(myjstr)
    if (myi.eq.1 .and. myj.eq.1) then
        print *, "reading ", varname, " from ",trim(filename)
    end if
    open(unit=iunit, file=trim(filename),status='old', action='read',iostat=ierr)
    if (ierr /= 0) then
        print *, "Error opening file ", trim(filename)
        call stopcm1
    end if
    ! print *, "reading "//trim(varname)//" from "//trim(filename)
    do
        read(iunit, '(A)', iostat=ios) line
        if (ios /= 0) then
            print *, "ERROR: varname not found ", varname
            call stopcm1
        end if
        if (trim(line) .eq. trim(varname)) then

            do idim1 = si, ei
                do idim2 = sj, ej
                    read (iunit,*,iostat=ios) value(idim1,idim2)
                    if (ios /= 0) then
                        print *, "ERROR: reading ", varname, " at ", idim1, idim2," from ", trim(filename)
                        call stopcm1
                    end if
                end do

            end do
            exit
        end if
    end do

    close(iunit)
end subroutine read2Dinteger

subroutine getbdytime(t_bdy_curr,t_bdy_next,t_cm1)
    integer, intent(out) :: t_bdy_curr,t_bdy_next
    double precision, intent(in) :: t_cm1
    integer :: time
    character(100) :: filename
    integer :: unit, ios, ierr
    integer :: iunit = 240422
    filename = "DATA_BDY/time"
    open(unit=iunit, file=trim(filename),status='old', action='read',iostat=ierr)
    if (ierr /= 0) then
        print *, "Error opening file ", trim(filename)
        call stopcm1
    end if
    ! print *, "reading "//trim(varname)//" from "//trim(filename)
    do
        read(iunit, *, iostat=ios) time
        if (ios /= 0) then

            print *, "ERROR: CM1 time is not in range of boundary time "
            print *, "CM1 time:", t_cm1
            call stopcm1

        end if
        if (time+0d0 <= t_cm1) then
            t_bdy_curr = time
        end if
        if (time+0d0 > t_cm1) then
            t_bdy_next = time
            exit
        end if
    end do

    close(iunit)
end subroutine getbdytime

end module read_wrf