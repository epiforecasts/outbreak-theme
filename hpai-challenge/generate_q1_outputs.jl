## Q1 descriptive output generation for HPAI modelling challenge
## Generates summary tables and epicurves from phase-1 case data

using CSV, DataFrames, Dates, Statistics

basedir = @__DIR__
datadir = joinpath(basedir, "data")
outdir = joinpath(basedir, "output")
mkpath(outdir)

# ── Load data ──────────────────────────────────────────────────────────────────

cases = CSV.read(
    joinpath(datadir, "phase-1", "cases.csv"), DataFrame;
    dateformat="yyyy-mm-dd"
)
population = CSV.read(
    joinpath(datadir, "canonical", "population.csv"), DataFrame
)
activity = CSV.read(
    joinpath(datadir, "phase-1", "activity.csv"), DataFrame;
    dateformat="yyyy-mm-dd"
)

# Merge case data with population metadata
merged = leftjoin(cases, population; on=:farm_id)

# Identify active farms: those with at least one record where date_end is missing
active_farms = unique(activity[ismissing.(activity.date_end), :farm_id])

# ── 1. table_cases_by_type.csv ─────────────────────────────────────────────────
# Distribution by species and production type

println("Generating table_cases_by_type.csv …")

table_rows = DataFrame(
    species=String[], production=String[],
    n_cases=Int[], n_farms=Int[], n_active_farms=Int[],
    attack_rate_pct=Float64[], attack_rate_active_pct=Float64[]
)

for g in groupby(population, [:species, :production])
    sp = first(g.species)
    pr = first(g.production)
    n_farms_total = nrow(g)
    farm_ids = g.farm_id
    n_active = count(fid -> fid in active_farms, farm_ids)
    n_cases_here = count(fid -> fid in Set(merged.farm_id), farm_ids)

    ar_pct = round(100.0 * n_cases_here / n_farms_total; digits=2)
    ar_active_pct = n_active > 0 ? round(100.0 * n_cases_here / n_active; digits=2) : 0.0

    push!(table_rows, (sp, pr, n_cases_here, n_farms_total, n_active, ar_pct, ar_active_pct))
end

sort!(table_rows, [:species, :production])
CSV.write(joinpath(outdir, "table_cases_by_type.csv"), table_rows)

# ── 2. epicurve_confirmed.csv ──────────────────────────────────────────────────
# Daily incidence by date_confirmed, split by chicken/duck

println("Generating epicurve_confirmed.csv …")

date_range = minimum(merged.date_confirmed):Day(1):maximum(merged.date_confirmed)

epi_confirmed = DataFrame(
    date=Date[], n_total=Int[], n_chicken=Int[], n_duck=Int[],
    cum_total=Int[], cum_chicken=Int[], cum_duck=Int[]
)

let cum_t = 0, cum_c = 0, cum_d = 0
    for d in date_range
        day_cases = filter(r -> r.date_confirmed == d, merged)
        nt = nrow(day_cases)
        nc = count(r -> r.species == "chicken", eachrow(day_cases))
        nd = count(r -> r.species == "duck", eachrow(day_cases))
        cum_t += nt
        cum_c += nc
        cum_d += nd
        push!(epi_confirmed, (d, nt, nc, nd, cum_t, cum_c, cum_d))
    end
end

CSV.write(joinpath(outdir, "epicurve_confirmed.csv"), epi_confirmed)

# ── 3. epicurve_suspicious.csv ─────────────────────────────────────────────────
# Daily incidence by date_suspicious for passive surveillance only

println("Generating epicurve_suspicious.csv …")

passive = filter(r -> r.detection_method == "passive", merged)

date_range_susp = minimum(passive.date_suspicious):Day(1):maximum(passive.date_suspicious)

epi_suspicious = DataFrame(
    date=Date[], n_total=Int[], n_chicken=Int[], n_duck=Int[]
)

for d in date_range_susp
    day_cases = filter(r -> r.date_suspicious == d, passive)
    nt = nrow(day_cases)
    nc = count(r -> r.species == "chicken", eachrow(day_cases))
    nd = count(r -> r.species == "duck", eachrow(day_cases))
    push!(epi_suspicious, (d, nt, nc, nd))
end

CSV.write(joinpath(outdir, "epicurve_suspicious.csv"), epi_suspicious)

# ── 4. spatial_cases.csv ───────────────────────────────────────────────────────
# All case farms with coordinates and metadata

println("Generating spatial_cases.csv …")

spatial = select(merged,
    :farm_id, :x, :y, :county, :district, :species, :production, :capacity,
    :date_suspicious, :date_confirmed, :detection_method, :cull_status
)
sort!(spatial, :date_confirmed)

CSV.write(joinpath(outdir, "spatial_cases.csv"), spatial)

# ── 5. cases_by_county.csv ────────────────────────────────────────────────────
# County-level summary

println("Generating cases_by_county.csv …")

county_rows = DataFrame(
    county=String[], n_cases=Int[], n_farms=Int[], n_active_farms=Int[],
    attack_rate_pct=Float64[], n_chicken_cases=Int[], n_duck_cases=Int[]
)

for g in groupby(population, :county)
    cty = first(g.county)
    n_farms_total = nrow(g)
    farm_ids = g.farm_id
    n_active = count(fid -> fid in active_farms, farm_ids)

    # Cases in this county
    county_cases = filter(r -> r.county == cty, merged)
    nc = nrow(county_cases)
    n_chick = count(r -> r.species == "chicken", eachrow(county_cases))
    n_duck = count(r -> r.species == "duck", eachrow(county_cases))

    ar_pct = round(100.0 * nc / n_farms_total; digits=2)

    push!(county_rows, (cty, nc, n_farms_total, n_active, ar_pct, n_chick, n_duck))
end

sort!(county_rows, :county)
CSV.write(joinpath(outdir, "cases_by_county.csv"), county_rows)

# ── 6. epicurve_by_production.csv ─────────────────────────────────────────────
# Daily incidence by date_confirmed broken down by production type

println("Generating epicurve_by_production.csv …")

prod_types = ["broiler_1", "broiler_2", "layer", "conventional", "organic"]
date_range_prod = minimum(merged.date_confirmed):Day(1):maximum(merged.date_confirmed)

epi_prod = DataFrame(date=Date[])
for pt in prod_types
    epi_prod[!, Symbol(pt)] = Int[]
end
epi_prod[!, :total] = Int[]

for d in date_range_prod
    day_cases = filter(r -> r.date_confirmed == d, merged)
    row = Dict{Symbol,Any}(:date => d)
    let total = 0
        for pt in prod_types
            n = count(r -> r.production == pt, eachrow(day_cases))
            row[Symbol(pt)] = n
            total += n
        end
        row[:total] = total
    end
    push!(epi_prod, row)
end

CSV.write(joinpath(outdir, "epicurve_by_production.csv"), epi_prod)

println("All Q1 descriptive output files generated in: $outdir")
