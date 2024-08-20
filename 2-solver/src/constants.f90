!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module constants

    use decomp_2d, only : mytype
    use param,     only : onehundredeighty

    ! Mathematical constants
    real(mytype), parameter :: pi = 3.14159265358979323846_mytype
    real(mytype), parameter :: conrad = pi / onehundredeighty
    real(mytype), parameter :: condeg = onehundredeighty / pi

end module constants
