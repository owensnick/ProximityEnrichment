module ProximityEnrichment

using GenomicFeatures, Statistics, StatsBase, HypothesisTests




"""
    intervaldist(intA, intB)

    intA        Genomic Interval A
    intB        Genomic Interval B

    returns distance between two intervals, 0 if intersecting.
"""
function intervaldist(intA, intB)
    (seqname(intA) != seqname(intB)) && return typemax(Int)

    ΔAB = leftposition(intB) - rightposition(intA)
    ΔBA = leftposition(intA) - rightposition(intB)

    (ΔBA ≤ 0) && (ΔAB ≤ 0) && return 0

    min(abs(ΔAB), abs(ΔBA))
end

"""
    closestinterval(intervals, i)

    intervals       Set of genomic locations
    i               Single genomic location

    returns distance between i and member of intervals closest to i
"""
function closestinterval(intervals, i)

    ind = searchsorted(intervals, i)

    Δ = typemax(Int)
    (1 ≤ ind.start ≤ length(intervals)) && (Δ = min(Δ, intervaldist(intervals[ind.start], i)))
    (1 ≤ ind.stop  ≤ length(intervals)) && (Δ = min(Δ, intervaldist(intervals[ind.stop],  i)))

    ifelse(Δ == typemax(Int), -1, Δ)
end


"""
    closest_tss_peak(gene_intervals, peak_intervals)

    gene_intervals   Intervals describing TSS
    peak_intervals   Intervals describing peak positions

    returns vector of distances
"""
closest_tss_peak(gene_intervals, peak_intervals) = [closestinterval(peak_intervals, iv) for iv in gene_intervals]


"""
    getx(d, xp)

    Aligns distances between peaks and genes that exist with positions to calculate fisher

    returns new_x and indexing of xp
"""
function getx(d, xp)
    ud = sort!(unique(d[d .> xp[1]]))
    dud = diff(ud)
    δ = step(xp)

    ind = findfirst(d -> d < δ, dud)

    if isnothing(ind)
        return xp, ind
    else
        xpi = Int(floor((ud[ind] - xp[1])/δ)) + 2
        return [ud[1:ind] ; xp[xpi:end]], ind
    end
end

"""
    fishertest(indA, indB)

    indA   logical vector selecting items in list A
    indB   logical vector selecting items in list B

    returns pvalues and oddsratio
"""
function fishertest(indA, indB)
    a = sum( indA  .&  indB )
    b = sum( indA  .& .!indB )
    c = sum(.!indA .&  indB )
    d = sum(.!indA .& .!indB )

    if iszero(a) || iszero(b) || iszero(c) || iszero(d)
       return 1.0, 1.0
    end
    or = (a/c)/(b/d)

    pvalue(FisherExactTest(a, b, c, d), tail=:right), or
end


"""
    fisherdist(xp, genes, peaks, geneind=trues(length(genes)), peakind=trues(length(peaks)); trans = x -> log10(x + 1), label="")

    xp          points to calculate enrichment at
    genes       Vector of Intervals describing TSS coordinates
    peaks       Vector of Intervals describing peaks
    geneind     logical vector selecting the foreground, i.e. genes of interest
    peakind     logical vector selecting the peaks to


    returns named tuple of relevant stats

"""
function fisherdist(xp, genes, peaks, geneind=trues(length(genes)), peakind=trues(length(peaks)); trans = x -> log10(x + 1), label="")
    ### Calculate distance closest peak to each tss
    Δb = trans.(closest_tss_peak(genes, peaks[peakind]))     ## background
    Δg = Δb[geneind]                                         ## foreground

    ### Calculate subset of xp to avoid calculating pvalues at distances below/above closest/furtherest gene
    xg, xgi = getx(Δg, xp)

    ### enrichment over background, not necessary for fisher calculation, but useful to plot
    fc = ecdf(Δg)(xg)./ecdf(Δb)(xg)

    ### Fisher tests for number of genes within δ against gene list
    count = [sum(Δb .≤ δ) for δ ∈ xg]
    ft    = [fishertest(geneind, Δb .≤ δ) for δ ∈ xg]

    pvalue = first.(ft)
    or     = last.(ft)
    (xg=xg, xgi=xgi, fc=fc, pvalue=pvalue, or=or, n=sum(geneind), np=size(peaks[peakind]), count=count, label=label)
end

end # module
