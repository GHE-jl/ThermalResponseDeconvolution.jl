using Test
using DelimitedFiles
using ThermalResponseDeconvolution

@testset "ThermalResponseDeconvolution.jl" begin

    @testset "convolution" begin
        f = [1.0, 2.0, 3.0]
        g = [1.0, 0.0, 0.0]

        # Convolving with a unit impulse returns f unchanged
        @test convolution(f, g) ≈ f

        # Exact discrete convolution formula, truncated to length(f)
        f2, g2 = [1.0, 2.0, 3.0], [4.0, 5.0, 6.0]
        expected = [
            f2[1] * g2[1],
            f2[1] * g2[2] + f2[2] * g2[1],
            f2[1] * g2[3] + f2[2] * g2[2] + f2[3] * g2[1],
        ]
        @test convolution(f2, g2) ≈ expected

        # Result length always matches length(f)
        @test length(convolution(f2, g2)) == length(f2)

        # Linearity: conv(f, a*g1 + b*g2) == a*conv(f,g1) + b*conv(f,g2)
        g3 = [0.5, 1.5, -1.0]
        a, b = 2.0, 3.0
        @test convolution(f2, a .* g2 .+ b .* g3) ≈ a .* convolution(f2, g2) .+ b .* convolution(f2, g3)

        # Mismatched lengths error
        @test_throws ArgumentError convolution([1.0, 2.0], [1.0, 2.0, 3.0])
    end

    @testset "rms" begin
        @test rms(zeros(5)) == 0.0
        @test rms(fill(3.0, 10)) ≈ 3.0
        @test rms([3.0, -4.0]) ≈ sqrt((9.0 + 16.0) / 2)
    end

    @testset "set_nodes" begin
        nData, n0 = 10080, 50
        id = set_nodes(nData, n0)

        @test length(id) == n0
        @test length(unique(id)) == n0
        @test issorted(id)
        @test id[1] == 1
        @test id[end] == nData
        @test all(1 .<= id .<= nData)

        # Requesting more nodes than data points is not meaningful and must not hang
        @test length(set_nodes(10, 5)) == 5
    end

    @testset "deconvolution on sample data" begin
        data_dir = joinpath(@__DIR__, "..", "data")
        data, _ = readdlm(joinpath(data_dir, "data.csv"), ',', header=true)
        gref, _ = readdlm(joinpath(data_dir, "gRef.csv"), ',', header=true)

        t = Float64.(data[:, 1])
        Tin = Float64.(data[:, 2])
        Tout = Float64.(data[:, 3])
        gRef = Float64.(gref[:, 2])

        f = diff([0.0; Tin .- Tout])
        Texp = Tout .- Tout[1]
        n = 50

        ĝ, gOpt = deconvolution(t, f, Texp; n=n, c=2)

        # Same length as the input data, node values reduced to n
        @test length(ĝ) == length(Texp)
        @test length(gOpt) == n

        # Recovers the reference thermal response function closely
        @test rms(ĝ .- gRef) < 0.05

        # Reconstructing the temperature from ĝ reproduces the measured data closely
        T̂ = collect(convolution(f, ĝ))
        @test rms(T̂ .- Texp) < 0.1

        # The interpolation passes exactly through the optimized node values
        id = set_nodes(length(Texp), n)
        @test ĝ[id] ≈ gOpt atol=1e-6

        # Constraint choices run without error and produce a monotonic response
        for c in (0, 1, 2)
            ĝc, _ = deconvolution(t, f, Texp; n=n, c=c)
            @test length(ĝc) == length(Texp)
        end
        ĝ1, _ = deconvolution(t, f, Texp; n=n, c=1)
        @test all(diff(ĝ1) .>= -1e-8)   # positivity constraint: (near) non-decreasing
    end

end
