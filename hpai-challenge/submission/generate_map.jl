using CairoMakie
using CSV, DataFrames, JSON

const DATA_DIR = joinpath(@__DIR__, "..", "data")
const OUTPUT_DIR = @__DIR__

# Load population (all farms)
pop = CSV.read(joinpath(DATA_DIR, "canonical", "population.csv"), DataFrame)

# Load cases
cases = CSV.read(joinpath(@__DIR__, "..", "output", "spatial_cases.csv"), DataFrame)

# Load county boundaries
counties = JSON.parsefile(joinpath(DATA_DIR, "canonical", "counties_32626.geojson"))

# Load HRZ boundary
hrz = JSON.parsefile(joinpath(DATA_DIR, "canonical", "hrz_32626.geojson"))

# Helper: extract polygon rings from a GeoJSON geometry (returns Vector of Nx2 matrices)
function extract_rings(geom)
    rings = Vector{Matrix{Float64}}()
    coords = geom["coordinates"]
    polygons = geom["type"] == "MultiPolygon" ? coords : [coords]
    for polygon in polygons
        for ring in polygon
            mat = reduce(vcat, [Float64[p[1] p[2]] for p in ring])
            push!(rings, mat)
        end
    end
    return rings
end

# Set up figure
fig = Figure(size = (700, 900))
ax = Axis(fig[1, 1];
    xlabel = "Easting (m)",
    ylabel = "Northing (m)",
    title = "HPAI confirmed cases — Jolly Island (Phase 1)\n22 Dec 2025 – 13 Jan 2026, N = 103",
    aspect = DataAspect(),
)

# Draw HRZ shading
for feature in hrz["features"]
    for ring in extract_rings(feature["geometry"])
        poly!(ax, ring[:, 1], ring[:, 2]; color = (:red, 0.08), strokecolor = :red,
            strokewidth = 1.0, linestyle = :dash)
    end
end

# Draw county boundaries
for feature in counties["features"]
    for ring in extract_rings(feature["geometry"])
        lines!(ax, ring[:, 1], ring[:, 2]; color = :gray60, linewidth = 0.5)
    end
end

# Plot all farms as small grey dots
scatter!(ax, pop.x, pop.y; markersize = 1, color = (:gray80, 0.3), rasterize = true)

# Plot case farms by species
chicken = filter(:species => ==("chicken"), cases)
duck = filter(:species => ==("duck"), cases)

sc = scatter!(ax, chicken.x, chicken.y;
    markersize = 8, color = :red, marker = :circle,
    strokecolor = :black, strokewidth = 0.3,
    label = "Chicken ($(nrow(chicken)))")

sd = scatter!(ax, duck.x, duck.y;
    markersize = 8, color = :steelblue, marker = :rect,
    strokecolor = :black, strokewidth = 0.3,
    label = "Duck ($(nrow(duck)))")

# Legend — use axislegend which pulls from labelled plot elements
axislegend(ax; position = :lt, labelsize = 12)

save(joinpath(OUTPUT_DIR, "q1_spatial_map.pdf"), fig)
println("Saved q1_spatial_map.pdf")
