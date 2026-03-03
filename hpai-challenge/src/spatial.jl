# ── 1. Spatial utilities ──────────────────────────────────────────────────────
# Point-in-polygon, HRZ membership, neighbour precomputation.

using JSON

# Ray-casting algorithm for point-in-polygon test
function point_in_polygon(px::Float64, py::Float64, ring::Vector)::Bool
    n = length(ring)
    inside = false
    j = n
    @inbounds for i in 1:n
        xi, yi = ring[i][1], ring[i][2]
        xj, yj = ring[j][1], ring[j][2]
        if ((yi > py) != (yj > py)) &&
           (px < (xj - xi) * (py - yi) / (yj - yi) + xi)
            inside = !inside
        end
        j = i
    end
    return inside
end

function point_in_multipolygon(px::Float64, py::Float64, mp_coords::Vector)::Bool
    for polygon in mp_coords
        # polygon[1] is the outer ring; polygon[2:end] are holes
        outer = polygon[1]
        if point_in_polygon(px, py, outer)
            in_hole = false
            for h in 2:length(polygon)
                if point_in_polygon(px, py, polygon[h])
                    in_hole = true
                    break
                end
            end
            if !in_hole
                return true
            end
        end
    end
    return false
end

"""
    identify_hrz_farms(x, y) -> BitVector

Returns a BitVector indicating which farms fall inside the HRZ boundary.
`x` and `y` are coordinate vectors (EPSG:32626).
"""
function identify_hrz_farms(x::Vector{Float64}, y::Vector{Float64})
    geojson = JSON.parsefile(HRZ_GEOJSON)
    features = geojson["features"]

    hrz = falses(length(x))
    for feat in features
        geom = feat["geometry"]
        coords = geom["coordinates"]
        gtype = geom["type"]
        for i in eachindex(x)
            hrz[i] && continue
            if gtype == "MultiPolygon"
                hrz[i] = point_in_multipolygon(x[i], y[i], coords)
            elseif gtype == "Polygon"
                hrz[i] = point_in_polygon(x[i], y[i], coords[1])
            end
        end
    end
    return hrz
end

# ── 2. Neighbour precomputation ──────────────────────────────────────────────

"""
    compute_neighbours(x, y; max_dist) -> (nbr_idx, nbr_dist)

Brute-force O(N²) pairwise distance computation within `max_dist` metres.
Returns adjacency lists of neighbour indices and distances.
"""
function compute_neighbours(x::Vector{Float64}, y::Vector{Float64};
                            max_dist::Float64 = 50_000.0)
    N = length(x)
    nbr_idx = [Int[] for _ in 1:N]
    nbr_dist = [Float64[] for _ in 1:N]
    max_dist_sq = max_dist^2

    @inbounds for i in 1:N
        for j in (i+1):N
            dx = x[i] - x[j]
            dy = y[i] - y[j]
            dsq = dx * dx + dy * dy
            if dsq <= max_dist_sq
                d = sqrt(dsq)
                push!(nbr_idx[i], j)
                push!(nbr_dist[i], d)
                push!(nbr_idx[j], i)
                push!(nbr_dist[j], d)
            end
        end
    end
    return nbr_idx, nbr_dist
end
