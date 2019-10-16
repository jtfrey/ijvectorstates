!
!
!
Module IJVectorStates
Implicit None

    Integer, Parameter          :: IJVectorStateUnchecked = 3
    Integer, Parameter          :: IJVectorStateTrue = 1
    Integer, Parameter          :: IJVectorStateFalse = 0

    Type IJVectorStatesType
        Integer                 :: n_vectors;
        Integer                 :: n_per_integer;
        Integer, Allocatable    :: bit_array(:);
    End Type IJVectorStatesType


    Public                      :: IJVectorStatesInit, IJVectorStatesCheck, IJVectorStatesCheckQuick
    Private                     :: IJVectorStatesCompareIJ

Contains

    Subroutine IJVectorStatesInit(tsa, n_vectors)
        Implicit None

        Type(IJVectorStatesType), Intent(InOut) :: tsa
        Integer, Intent(In)                     :: n_vectors
        Integer                                 :: bits_per_integer
        Integer                                 :: array_dim, integer_dim, init_value, i

        ! Already allocated?
        if ( Allocated(tsa%bit_array) ) Deallocate(tsa%bit_array)

        ! Set the vector count on the object:
        tsa%n_vectors = n_vectors

        ! We don't need the diagonal, we implicitly know a vector is
        ! equal to itself:
        array_dim = (n_vectors * (n_vectors - 1)) / 2

        ! How many bits fit in an integer?
        bits_per_integer = bit_size(bits_per_integer)

        ! How many 2-bit element tags per integer?
        tsa%n_per_integer = bits_per_integer / 2

        ! Initialize the initial value:
        init_value = IJVectorStateUnchecked
        Do i = 2, tsa%n_per_integer
            init_value = IOr(IShft(init_value, 2), IJVectorStateUnchecked)
        End Do

        ! How many integers do we need?
        array_dim = (array_dim + tsa%n_per_integer) / tsa%n_per_integer

        ! Allocate and initialize the array:
        Allocate(tsa%bit_array(array_dim))
        tsa%bit_array(:) = init_value
    End Subroutine

    Function IJVectorStatesCheckQuick(tsa, IIndex, JIndex)
        Implicit None

        Integer                                 :: IJVectorStatesCheckQuick
        Type(IJVectorStatesType), Intent(InOut) :: tsa
        Integer, Intent(In)                     :: IIndex, JIndex
        Integer                                 :: ij_index, integer_index, bit_index

        if ( IIndex == JIndex ) then
            ! Index is the same, then it's implicitly enabled:
            IJVectorStatesCheckQuick = IJVectorStateTrue
        else
            ! Calculate the index of vectors IIndex, JIndex:
            if ( IIndex < JIndex ) then
                ij_index = (tsa%n_vectors*(tsa%n_vectors-1)/2) - ((tsa%n_vectors-(IIndex - 1))*((tsa%n_vectors-(IIndex - 1))-1))/2 + (JIndex - 1) - IIndex
            else
                ij_index = (tsa%n_vectors*(tsa%n_vectors-1)/2) - ((tsa%n_vectors-(JIndex - 1))*((tsa%n_vectors-(JIndex - 1))-1))/2 + (IIndex - 1) - JIndex
            end if

            ! Translate the I-J index into an integer and bit offset in the
            ! bit array:
            integer_index = 1 + (ij_index / tsa%n_per_integer)
            bit_index = 2 * Modulo(ij_index, tsa%n_per_integer)

            ! Isolate the state value:
            IJVectorStatesCheckQuick = IAnd(IShft(tsa%bit_array(integer_index), -bit_index), 3)
        end if
    End Function

    Function IJVectorStatesCheck(tsa, N, I, IIndex, J, JIndex)
        Implicit None

        Logical                                 :: IJVectorStatesCheck
        Type(IJVectorStatesType), Intent(InOut) :: tsa
        Integer, Intent(In)                     :: N, I(:), IIndex, J(:), JIndex
        Integer                                 :: ij_index, integer_index, bit_index
        Integer                                 :: state_value

        if ( IIndex == JIndex ) then
            ! Index is the same, then it's implicitly enabled:
            IJVectorStatesCheck = .True.
        else
            ! Calculate the index of vectors IIndex, JIndex:
            if ( IIndex < JIndex ) then
                ij_index = (tsa%n_vectors*(tsa%n_vectors-1)/2) - ((tsa%n_vectors-(IIndex - 1))*((tsa%n_vectors-(IIndex - 1))-1))/2 + (JIndex - 1) - IIndex
            else
                ij_index = (tsa%n_vectors*(tsa%n_vectors-1)/2) - ((tsa%n_vectors-(JIndex - 1))*((tsa%n_vectors-(JIndex - 1))-1))/2 + (IIndex - 1) - JIndex
            end if

            ! Translate the I-J index into an integer and bit offset in the
            ! bit array:
            integer_index = 1 + (ij_index / tsa%n_per_integer)
            bit_index = 2 * Modulo(ij_index, tsa%n_per_integer)

            ! Isolate the state value:
            state_value = IAnd(IShft(tsa%bit_array(integer_index), -bit_index), 3)

            ! Already set or not?
            if (state_value == IJVectorStateUnchecked) then
                ! We haven't checked these two vectors yet:
                if (IJVectorStatesCompareIJ(N, I, J) <= 2) then
                    state_value = IJVectorStateTrue
                else
                    state_value = IJVectorStateFalse
                end if

                ! Cache that value back into the bit array:
                tsa%bit_array(integer_index) = IOr(IAnd(tsa%bit_array(integer_index), Not(IShft(3, bit_index))), IShft(state_value, bit_index))
            end if

            ! Set return value according to state_value:
            if (state_value == IJVectorStateFalse) then
                IJVectorStatesCheck = .False.
            else
                IJVectorStatesCheck = .True.
            end if
        end if
    End Function

    Function IJVectorStatesCompareIJ(N, I, J)
        Implicit None

        Integer                 :: IJVectorStatesCompareIJ
        Integer, Intent(In)     :: N, I(:), J(:)
        Integer                 :: Ii, Ji, diffs

        ! Initialize indices and diff count:
        Ii = 1
        Ji = 1
        diffs = 0

        ! Loop over the two vectors simultaneously:
        do while ((diffs < 3) .and. (Ii <= N) .and. (Ji < N))
            if (I(Ii) == J(Ji)) then
                Ii = Ii + 1
                Ji = Ji + 1
            else if (I(Ii) < J(Ji)) then
                ! this means I(Ii) is NOT in j, so we move ahead on i
                ! but NOT j and note the difference:
                diffs = diffs + 1
                Ii = Ii + 1
            else
                ! this means J(Ji) is NOT in i, so we move ahead on j
                ! but NOT i and note the difference:
                diffs = diffs + 1
                Ji = Ji + 1
            end if
        end do

        ! If we haven't yet hit 3 differences yet, then we can absorb however
        ! many indices were not traversed on each vector:
        if ((diffs < 3) .and. (Ii <= N)) then
            ! Additional diffs on i are (N - Ii + 1)
            diffs = Min(3, diffs + (N - Ii + 1))
        end if
        if ((diffs < 3) .and. (Ji <= N)) then
            ! Additional diffs on j are (N - Ji + 1)
            diffs = Min(3, diffs + (N - Ji + 1))
        end if
        IJVectorStatesCompareIJ = diffs
    End Function

End Module
