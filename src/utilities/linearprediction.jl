abstract type AbstractPredictionScheme end

# see PHYSICAL REVIEW B 79, 245101 (2009), BarthelWhite2009
# also Appendix in PHYSICAL REVIEW B 90, 115124 (2014), WolfSchollwock2014b
struct LinearPrediction{T<:Number} <: AbstractPredictionScheme
	stepsize::Float64
	nob::Int
	nfit::Int
	p::Int
	obs::Vector{T}
	ws::Vector{Float64}
	a::Vector{T}
end

function LinearPrediction(obs::Vector{T}, ws::Vector{<:Real}; stepsize::Real, nfit=length(obs), p::Int=div(nfit, 2)) where {T<:Number}
	nobs = length(obs)
	(nfit <= nobs) || throw(ArgumentError("nfit must be less than nobs")) 
	(p <= nfit) || throw(ArgumentError("p must be less than nfit"))
	(nobs == length(ws)) || throw(DimensionMismatch())

	Rmat = zeros(T, p, p)
	rvec = zeros(T, p)
	n_start = nobs - nfit + 1

	for n in n_start:nobs
		for i in 1:p
			for j in 1:p
				if (i < n) && (j < n)
					Rmat[j, i] += (conj(obs[n-j]) * obs[n-i]) / ws[n]
				end
			end
		end
	end
	# a diagonal shift to avoid numeric instability of inverting
	δ = 1.0e-6
	for i in 1:p
		Rmat[i, i] += δ
	end
	for n in n_start:nobs
		for j in 1:p
			if j < n
				rvec[j] += (conj(obs[n-j]) * obs[n]) / ws[n]
			end
		end
	end

	# println("Rmat is $Rmat")
	# println("rvec is $rvec")

	a = - Rmat \ rvec

	# println("a is $a")

	return LinearPrediction{T}(convert(Float64, stepsize), length(obs), nfit, p, copy(obs), ws, a)
end
LinearPrediction(obs::Vector{<:Number}; kwargs...) = LinearPrediction(obs, ones(length(obs)); kwargs...)

current_size(x::LinearPrediction) = length(x.obs)

function compute_next!(x::LinearPrediction{T}) where T
	obs_new = zero(T)
	n = current_size(x)
	for i in 1:x.p
		obs_new += x.a[i] * x.obs[n-i+1]
	end
	push!(x.obs, -obs_new)
end

function Base.getindex(x::LinearPrediction, j::Int)
	if j <= current_size(x)
		return x.obs[j]
	else
		compute_next!(x)
		return x[j]
	end	
end
Base.getindex(x::LinearPrediction, j::AbstractRange{Int64}) = [x[i] for i in j]

# make a simple linear extrapolation
function (x::LinearPrediction)(t::Real)
	(t >= zero(t)) || throw(ArgumentError("t must be nonnegative"))
	ns = t / x.stepsize
	n = floor(Int, ns)
	# (n ≈ ns) && return x[n]
	dif = ns - n
	(dif ≈ zero(t)) && return x[n+1]
	return x[n+1] * (1-dif) + x[n+2] * dif
end

function linear_predict(obs::Vector{<:Number}, stepsize::Real; nfit=length(obs), p::Int=div(nfit, 2), 
						δt::Real=stepsize, maxiter::Int=10000, tol::Real=1.0e-6, verbosity::Int=1)
	x = LinearPrediction(obs, stepsize=stepsize, nfit=nfit, p=p)
	iter = 0
	delta = abs(x.obs[end])
	while (delta > tol) && (iter < maxiter)
		compute_next!(x)
		iter += 1
		delta = abs(x.obs[end])
	end
    if (verbosity >= 2) && (iter < maxiter)
        println("early converge in $iter-th iterations with error $(delta)")
    end
    if (verbosity > 0) && (delta >= tol)
        println("fail to converge, required precision: $(tol), actual precision $delta in $iter iterations")
    end

    t_all = (length(x.obs)-1) * stepsize
    n = round(Int, t_all / δt)
    return [x((i-1) * δt) for i in 1:n+1]
end

