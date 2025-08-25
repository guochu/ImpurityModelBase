# bogliubov transformation of BCS bath

_dispersion(ϵ::Real, Δ::Number) = sqrt(ϵ^2 + abs2(Δ))

function _u(ϵ::Real, Δ::Number)
	(Δ == zero(Δ)) && return (ϵ >= 0) ? one(ϵ) : zero(ϵ)
	return sqrt((1 + ϵ/_dispersion(ϵ, Δ))/2)
end 
function _v(ϵ::Real, Δ::Real)
	(Δ == zero(Δ)) && return (ϵ >= 0) ? zero(ϵ) : one(ϵ)
	return sqrt((1 - ϵ/_dispersion(ϵ, Δ))/2)
end 

function _v(ϵ::Real, Δ::Complex)
	(Δ == zero(Δ)) && return (ϵ >= 0) ? zero(ϵ) : one(ϵ)
	phi = angle(Δ)
	return sqrt((1 - ϵ/_dispersion(ϵ, Δ))/2) * exp(im*phi)
end 


# function _u2(ϵ, Δ) 
# 	x = _u(ϵ, Δ)
# 	return conj(x) * x
# end
# function _v2(ϵ, Δ)
# 	x = _v(ϵ, Δ)
# 	return conj(x) * x
# end 
# function _uv(ϵ, Δ)
# 	return conj(_u(ϵ, Δ)) * _v(ϵ, Δ)
# end 
# function _vu(ϵ, Δ)
# 	return conj(_v(ϵ, Δ)) * _u(ϵ, Δ)
# end



function bogoliubov_freehamiltonian(h::BCSToulouse)
	T = eltype(h)
	L = div(num_sites(h), 2)
	ham = GenericQuadraticHamiltonian(T, num_sites(h))
	bath = h.bath
	ws = frequencies(bath)
	Δ = bath.Δ
	ϵ_d = h.ϵ_d

	push!(ham, adaga(1,1, coeff=ϵ_d))
	push!(ham, adaga(2,2, coeff=ϵ_d))

	for i in 1:L-1
		wj = _dispersion(ws[i], Δ)
		# wj = ws[i]
		push!(ham, adaga(2i+1,2i+1, coeff=wj))
		push!(ham, adaga(2i+2, 2i+2, coeff=wj))

	end
	return ham
end

function bogoliubov_hamiltonian(h::BCSToulouse)
	ham = bogoliubov_freehamiltonian(h)
	bath = h.bath
	ws, fs = frequencies(bath), spectrumvalues(bath)
	Δ = bath.Δ
	L = div(num_sites(h), 2)
	for i in 1:L-1
		uk = _u(ws[i], Δ)
		vk = _v(ws[i], Δ)
		t = adaga(1, 2i+1, coeff=fs[i]*conj(uk))
		push!(ham, t)
		push!(ham, t')
		t = aa(2, 2i+1, coeff=fs[i]*conj(vk))
		push!(ham, t)
		push!(ham, t')
		t = adaga(2, 2i+2, coeff=fs[i]*conj(uk))
		push!(ham, t)
		push!(ham, t')
		t = aa(1, 2i+2, coeff=-fs[i]*conj(vk))
		push!(ham, t)
		push!(ham, t')
	end
	return ham
end
bogoliubov_cmatrix(h::BCSToulouse) = cmatrix(bogoliubov_hamiltonian(h))


function bogoliubov_hamiltonian(bath::AbstractDiscreteBCSBath)
	ws = frequencies(bath)
	Δ = bath.Δ
	L = div(num_sites(bath), 2)

	T = eltype(bath)
	data = QuadraticTerm{T}[]
	for i in 1:L-1
		v = _dispersion(ws[i], Δ)
		push!(ham, adaga(2i-1, 2i-1, coeff=v))
		push!(ham, adaga(2i, 2i, coeff=v))
	end
	return ham
end
bogoliubov_cmatrix(h::AbstractDiscreteBCSBath) = cmatrix(bogoliubov_hamiltonian(h))


function bogoliubov_thermocdm(m::BCSToulouse)
	h = bogoliubov_cmatrix(m)
	cache = eigencache(h)
	return fermionicthermocdm(cache, β=m.bath.β)
end

bogoliubov_separablecdm(m::BCSToulouse; nsys::Real=fermidirac(m.bath.β, 0, m.ϵ_d)) = bogoliubov_separablecdm(m, nsys)
function bogoliubov_separablecdm(m::BCSToulouse, nsys::Real)
	L = div(num_sites(m), 2)
	ρ₀ = thermocdm(m.bath)

	ρ = zeros(eltype(ρ₀), 4L, 4L)
	ρ[1, 1] = nsys
	ρ[2, 2] = nsys
	ρ[2L+1, 2L+1] = 1-nsys
	ρ[2L+2, 2L+2] = 1-nsys	

	bath = m.bath
	ws, fs = frequencies(bath), spectrumvalues(bath)
	ws2 = [_dispersion(w, f) for (w, f) in zip(ws, fs)]
	ns = [fermidirac(bath.β, bath.μ, w) for w in ws2]
	N = div(num_sites(bath), 2)
	for i in 1:N
		ρ[2i+1, 2i+1] = ns[i]
		ρ[2i+2, 2i+2] = ns[i]

		ρ[2L+2i+1, 2L+2i+1] = 1-ns[i]
		ρ[2L+2i+2, 2L+2i+2] = 1-ns[i]		
	end
	return ρ
end
