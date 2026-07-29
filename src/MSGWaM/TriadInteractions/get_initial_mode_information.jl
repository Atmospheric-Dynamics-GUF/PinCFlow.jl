function get_initial_mode_information end

function get_initial_mode_information(
    peak_action::AbstractMatrix{<:AbstractFloat};
    support_tol::AbstractFloat = 0.0,
    diagonal_connectivity::Bool = false,
)

    nkp, nm =
        size(peak_action)

    visited =
        falses(nkp, nm)

    stack =
        CartesianIndex{2}[]

    ModeInformation =
        NamedTuple{
            (
                :cell_count,
                :peak_index,
                :kpi_min,
                :kpi_max,
                :mi_min,
                :mi_max,
                :peak_action,
            ),
            Tuple{
                Int,
                CartesianIndex{2},
                Int,
                Int,
                Int,
                Int,
                Float64,
            },
        }

    mode_information =
        ModeInformation[]

    #----------------------------------------------------------
    # Spectral-cell connectivity
    #----------------------------------------------------------

    if diagonal_connectivity
        neighbour_offsets = (
            CartesianIndex(-1, -1),
            CartesianIndex(-1, 0),
            CartesianIndex(-1, 1),
            CartesianIndex(0, -1),
            CartesianIndex(0, 1),
            CartesianIndex(1, -1),
            CartesianIndex(1, 0),
            CartesianIndex(1, 1),
        )
    else
        neighbour_offsets = (
            CartesianIndex(-1, 0),
            CartesianIndex(1, 0),
            CartesianIndex(0, -1),
            CartesianIndex(0, 1),
        )
    end

    #----------------------------------------------------------
    # Identify all connected active components
    #----------------------------------------------------------

    @inbounds for mi in 1:nm,
        kpi in 1:nkp

        if visited[kpi, mi] ||
           peak_action[kpi, mi] <= support_tol

            continue
        end

        empty!(stack)

        initial_index =
            CartesianIndex(kpi, mi)

        push!(
            stack,
            initial_index,
        )

        visited[initial_index] =
            true

        mode_cell_count = 0

        mode_peak_action = 0.0

        mode_peak_index =
            initial_index

        kpi_min = kpi
        kpi_max = kpi

        mi_min = mi
        mi_max = mi

        #------------------------------------------------------
        # Depth-first search of one connected component
        #------------------------------------------------------

        while !isempty(stack)

            current_index =
                pop!(stack)

            current_kpi =
                current_index[1]

            current_mi =
                current_index[2]

            current_action =
                peak_action[current_index]

            mode_cell_count += 1

            if current_action > mode_peak_action

                mode_peak_action =
                    current_action

                mode_peak_index =
                    current_index
            end

            kpi_min =
                min(kpi_min, current_kpi)

            kpi_max =
                max(kpi_max, current_kpi)

            mi_min =
                min(mi_min, current_mi)

            mi_max =
                max(mi_max, current_mi)

            #--------------------------------------------------
            # Add unvisited active neighbours
            #--------------------------------------------------

            for neighbour_offset in neighbour_offsets

                neighbour_index =
                    current_index +
                    neighbour_offset

                neighbour_kpi =
                    neighbour_index[1]

                neighbour_mi =
                    neighbour_index[2]

                if neighbour_kpi < 1 ||
                   neighbour_kpi > nkp ||
                   neighbour_mi < 1 ||
                   neighbour_mi > nm

                    continue
                end

                if visited[neighbour_index] ||
                   peak_action[neighbour_index] <= support_tol

                    continue
                end

                visited[neighbour_index] =
                    true

                push!(
                    stack,
                    neighbour_index,
                )
            end
        end

        push!(
            mode_information,
            (
                cell_count =
                    mode_cell_count,

                peak_index =
                    mode_peak_index,

                kpi_min =
                    kpi_min,

                kpi_max =
                    kpi_max,

                mi_min =
                    mi_min,

                mi_max =
                    mi_max,

                peak_action =
                    Float64(mode_peak_action),
            ),
        )
    end

    return mode_information
end