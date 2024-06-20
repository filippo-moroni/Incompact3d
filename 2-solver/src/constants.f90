module constants

    use decomp_2d, only : mytype
    use param,     only : onehundredeighty

    ! Mathematical constants
    real(mytype), parameter :: pi = 3.14159265358979323846_mytype
    real(mytype), parameter :: conrad = pi / onehundredeighty
    real(mytype), parameter :: condeg = onehundredeighty / pi

end module constants
