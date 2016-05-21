function NaivePDE.periodic_makes_boundary!{T1 <: Number,n}(
  x:: Array{T1,n}, p :: MPI_areaConfig)
  return NaivePDE.periodic_makes_boundary!(x, p.OL)
end

function periodic_expand_boundary{T1 <: Number,n}(
  xcore:: Array{T1,n}, p :: MPI_areaConfig)
  return NaivePDE.periodic_expand_boundary(xcore, p.OL)
end

function NaivePDE.get_core{T1 <: Number,n}(
  x:: Array{T1,n}, p :: MPI_areaConfig)
  return NaivePDE.get_core(x, p.OL)
end

function NaivePDE.NaN_expand{T1 <: Number,n}(
  xcore:: Array{T1,n}, p :: MPI_areaConfig)
  return NaivePDE.NaN_expand(xcore, p.OL)
end
