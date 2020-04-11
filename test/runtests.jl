using ProximityEnrichment
using Test
using GenomicFeatures

intA = Interval("chr1", 10, 100, '.')
intB = Interval("chr1", 80, 200, '.')
intC = Interval("chr1", 110, 500, '.')
intD = Interval("chr2", 100, 200, '.')
intE = Interval("chr2", 400, 600, '.')
intF = Interval("chrX", 10, 600, '-')

ivs = [intA, intB, intC, intD, intE, intF]





@testset "ProximityEnrichment.jl" begin
    # Write your own tests here
    @test ProximityEnrichment.intervaldist(intA, intB) == 0
    @test ProximityEnrichment.intervaldist(intB, intA) == 0
    @test ProximityEnrichment.intervaldist(intA, intC) == 10
    @test ProximityEnrichment.intervaldist(intC, intA) == 10
    @test ProximityEnrichment.intervaldist(intB, intC) == 0

    @test ProximityEnrichment.intervaldist(Interval("chr1", 1:20, '.'), Interval("chr2", 1:20, '.')) == typemax(Int)
    @test ProximityEnrichment.intervaldist(Interval("chr1", 1:20,), Interval("chr1", 100:120)) == 80

    @test ProximityEnrichment.closestinterval(ivs, Interval("chr1", 1, 8, '.')) == 2
    @test ProximityEnrichment.closestinterval(ivs, Interval("chr1", 1, 9, '.')) == 1
    @test ProximityEnrichment.closestinterval(ivs, Interval("chr1", 5, 15, '.')) == 0
    @test ProximityEnrichment.closestinterval(ivs, Interval("chr1", 50, 150, '.')) == 0

    @test ProximityEnrichment.closestinterval(ivs, Interval("chr1", 600, 750, '.')) == 100
    @test ProximityEnrichment.closestinterval(ivs, Interval("chr2", 10, 15, '.')) == 85
    @test ProximityEnrichment.closestinterval(ivs, Interval("chr2", 300, 300, '.')) == 100
    @test ProximityEnrichment.closestinterval(ivs, Interval("chr2", 300, 300, '.')) == 100
    @test ProximityEnrichment.closestinterval(ivs, Interval("chr2", 1000, 2000, '+')) == 400
    @test ProximityEnrichment.closestinterval(ivs, Interval("chr3", 1000, 2000, '+')) == -1
    @test ProximityEnrichment.closestinterval(ivs, Interval("chrY", 1000, 2000, '+')) == -1

end


intA = Interval("chr1", 10, 100, '.')
intB = Interval("chr1", 110, 500, '.')
intC = Interval("chr1", 80, 200, '.')

@show ProximityEnrichment.intervaldist(intA, intB)
