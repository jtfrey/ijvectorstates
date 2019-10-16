Program FTest

    Use IJVectorStates

    Type(IJVectorStatesType)        :: ij_states
    Integer, Parameter              :: N = 33, M = 12
    Integer                         :: IJ(N, M), ii, ji, IS
    Logical                         :: S

    IJ = Reshape( &
            (/  40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, &
                40, 41, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, &
                40, 41, 42, 43, 44, 45, 46, 47, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, &
                40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 75, &
                40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 74, &
                40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 80, 81, 82, 83, 84, &
                43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 82, 83, 84, &
                43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 67, 68, 69, 70, 71, 72, 82, 83, 84, 85, &
                40, 41, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, &
                43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 82, 83, 87, &
                43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 67, 68, 69, 70, 71, 72, 82, 83, 84, 88, &
                40, 41, 42, 43, 44, 45, 46, 47, 48, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 64, 65, 66, 67, 68, 69, 70, 71, 72, 81, 83 /), &
                Shape(IJ) &
            )

    ! Initialize the IJ state array:
    Call IJVectorStatesInit(ij_states, M)

    ! What's the state array look like initially?
    Write(*,*) ''
    Write(*,*) 'Initial state array'
    Write(*,*) ''
    Do ii = 1, M
        Do ji = 1, M
            IS = IJVectorStatesCheckQuick(ij_states, ii, ji)
            Write(*,1000,advance="no") IS
        End Do
        Write(*,*) ''
    End Do

    !
    ! Generate a visual array of which IJ-combinations have <= 2 perturbations.
    ! Note that the calls to IJVectorStatesCheck() will initially have to actually compare
    ! the two vectors so they will take longer -- but the result will be cached so
    ! that subsequent calls against the same index-pair (i,j) or (j,i) will use the
    ! cached state.
    !
    ! IJVectorStatesCheck() returns .True. if the two vectors differ by [0,2] elements,
    ! .False. otherwise.
    !
    Write(*,*) ''
    Write(*,*) 'First pass with IJVectorStatesCheck()'
    Write(*,*) ''
    Do ii = 1, M
        Do ji = 1, M
            S = IJVectorStatesCheck(ij_states, N, IJ(:, ii), ii, IJ(:,ji), ji)
            if (S) then
                Write(*,1000,advance="no") 1
            else
                Write(*,1000,advance="no") 0
            end if
        End Do
        Write(*,*) ''
    End Do

1000 Format(I1,' ')

    ! Now this time we'll be hitting the cached values entirely:
    Write(*,*) ''
    Write(*,*) 'All cached values, hopefully (no 3s)'
    Write(*,*) ''
    Do ii = 1, M
        Do ji = 1, M
            IS = IJVectorStatesCheckQuick(ij_states, ii, ji)
            Write(*,1000,advance="no") IS
        End Do
        Write(*,*) ''
    End Do

End Program
