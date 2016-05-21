immutable MPI_2D_Config <: MPI_areaConfig
  ##  Voodoo numbers controlling data layout.
  ##  sNx :: No. X points in sub-grid.
  ##  sNy :: No. Y points in sub-grid.
  ##  OLx :: Overlap extent in X.
  ##  OLy :: Overlat extent in Y.
  ##  nPx :: No. of processes to use in X.
  ##  nPy :: No. of processes to use in Y.
  sN :: Tuple{Int64,Int64}
  OL :: Tuple{Int64,Int64}
  nP :: Tuple{Int64,Int64}
  rank :: Int64
  size :: Int64
  function MPI_2D_Config(n :: Tuple{Int64,Int64}, l :: Tuple{Int64,Int64}, p :: Tuple{Int64,Int64})
    @assert n[1] > l[1]
    @assert n[2] > l[2]
    comm = MPI.COMM_WORLD
    new(n,l,p, MPI.Comm_rank(comm), MPI.Comm_size(comm))
  end
end

function MPI2D_InitialArray!(p :: MPI_2D_Config )
  return fill(0.0, )
end

function MPI2D_GatherAll!(p :: MPI_2D_Config, T :: Array{Float64,2}, Twhole :: Array{Float64,2})
  comm = MPI.COMM_WORLD
  (OLx, OLy) = p.OL
  (sNx, sNy) = p.sN
  (nPx, nPy) = p.nP
  rank = p.rank

  if (rank!=0)
    ##  Sending to Coupler
    sreq = MPI.Isend(T,0, rank+32,comm)
    MPI.Waitall!([sreq;])
  else
    Ti = fill(NaN, size(T))
    Twhole[1:sNx, 1:sNy] = T[OLx+(1:sNx),OLy+(1:sNy)]
    for ip = 1:(nPx*nPy-1)
      ix = mod(ip,nPx)
      iy = round(Int64, (ip-ix)/nPx)
      rreq = MPI.Irecv!(Ti, ip, ip+32, comm)
      MPI.Waitall!([rreq;])
      Twhole[ix*sNx+ (1:sNx),iy*sNy+(1:sNy)] = Ti[OLx+(1:sNx),OLy + (1:sNy)]
    end
  end
  return Twhole
end

function MPI2D_DistributeCore!(p :: MPI_2D_Config, T :: Array{Float64,2},  Twhole :: Array{Float64,2})
  comm = MPI.COMM_WORLD
  (OLx, OLy) = p.OL
  (sNx, sNy) = p.sN
  (nPx, nPy) = p.nP
  rank = p.rank

  if (rank != 0)
    ##  Receiving From Coupler
    Ti = fill(NaN, p.sN)
    csreq = MPI.Irecv!(Ti,0, rank+32,comm)
    MPI.Waitall!([csreq;])
    T[OLx+(1:sNx), OLy+(1:sNy)] = Ti
  else
    Ti = fill(NaN, size(T))
    T[(1+OLx):(OLx+sNx),(1+OLy):(OLy+sNy)] =  Twhole[1:sNx, 1:sNy]
    ## Send to Workers
    for ip = 1:(nPx*nPy-1)
      ix = mod(ip,nPx)
      iy = round(Int64, (ip-ix)/nPx)
      Ti = fill(NaN, p.sN)
      Ti[:,:] = Twhole[ix*sNx+ (1:sNx),iy*sNy+(1:sNy)]
      sreq = MPI.Isend(Ti, ip, ip+32, comm)
      MPI.Waitall!([sreq;])
    end
  end
  return Twhole
end

function MPI2D_DistributeBoundary!(p :: MPI_2D_Config, T :: Array{Float64,2},  Twhole :: Array{Float64,2})
  comm = MPI.COMM_WORLD
  #= MPI.Barrier(comm) =#
  (OLx, OLy) = p.OL
  (sNx, sNy) = p.sN
  (nPx, nPy) = p.nP
  rank = p.rank

  if (rank!=0)

    crreq = Array(MPI.Request,4)
    Twest = fill(NaN, (OLx, sNy))
    Tnorth = fill(NaN,(sNx, OLy))
    Teast = fill(NaN, (OLx, sNy))
    Tsouth = fill(NaN, (sNx, OLy))

    ##  Receiving from Coupler
    crreq[1] = MPI.Irecv!(Twest,0, 4*rank+1, comm)
    crreq[2] = MPI.Irecv!(Tnorth,0, 4*rank+2, comm)
    crreq[3] = MPI.Irecv!(Teast,0, 4*rank+3, comm)
    crreq[4] = MPI.Irecv!(Tsouth,0, 4*rank+4, comm)
    MPI.Waitall!(crreq)
    T[1:OLx,OLy+(1:sNy)] = Twest
    T[OLx+(1:sNx),sNy+OLy+(1:OLy)] = Tnorth
    T[sNx+OLx+(1:OLx),OLy+(1:sNy)] = Teast
    T[OLx+(1:sNx),1:OLy] = Tsouth
  else

    Twest = Twhole[sNx*nPx-OLx+(1:OLx),1:sNy]
    Tnorth = Twhole[1:sNx, mod(sNy+(1:OLy)-1, sNy*nPy)+1]
    Teast  = Twhole[mod(sNx+(1:OLx)-1, sNx*nPx )+1,1:sNy]
    Tsouth = Twhole[1:sNx, mod(sNy*nPy-OLy+(1:OLy)-1,sNy*nPy)+1] 

    T[1:OLx,OLy+(1:sNy)] = Twest
    T[OLx+(1:sNx),sNy+OLy+(1:OLy)] = Tnorth
    T[sNx+OLx+(1:OLx),OLy+(1:sNy)] = Teast
    T[OLx+(1:sNx),1:OLy] = Tsouth

    sreq = Array(MPI.Request,4)
    for ip = 1:(nPx*nPy-1)
      ix = mod(ip,nPx)
      iy = round(Int64, (ip-ix)/nPx)

      # Sending West Boundary
      iwx = mod((ix-1), nPx)
      iwy = iy
      pendingx = (iwx+1)*sNx-OLx
      pendingy = (iwy)*sNy
      Twest = Twhole[pendingx+(1:OLx), pendingy+(1:sNy)]
      sreq[1] = MPI.Isend(Twest, ip, ip*4+1, comm)

      # Sending North Boundary
      iwx = ix
      iwy = mod(iy+1,nPy)
      pendingx = iwx*sNx
      pendingy = iwy*sNy 
      Tnorth = Twhole[pendingx+(1:sNx), pendingy+(1:OLy)]
      sreq[2] = MPI.Isend(Tnorth, ip, ip*4+2, comm)
      
      # Sending East Boundary
      iwx = mod(ix+1,nPx)
      iwy = iy
      pendingx = iwx*sNx
      pendingy = iwy*sNy 
      Twest = Twhole[pendingx+(1:OLx), pendingy+(1:sNy)]
      sreq[3] = MPI.Isend(Teast, ip, ip*4+3, comm)
      
      # Sending South Boundary
      iwx = ix
      iwy = mod(iy-1,nPy)
      pendingx = iwx*sNx
      pendingy = (iwy+1)*sNy-OLy
      Twest = Twhole[pendingx+(1:sNx), pendingy+(1:OLy)]
      sreq[4] = MPI.Isend(Tsouth, ip, ip*4+4, comm)
      
      MPI.Waitall!(sreq)
     
    end
  end
  return T
end

function MPI2D_GatherAll!(p :: MPI_2D_Config, T :: Array{Float64,3}, Twhole :: Array{Float64,3})
  for ii = 1:size(T,3)
    T1 = T[:,:,ii]
    Twhole1 = Twhole[:,:,ii]
    MPI2D_GatherAll!(p,T1, Twhole1)
  end
  return Twhole1
end

function MPI2D_DistributeCore!(p :: MPI_2D_Config, T :: Array{Float64,3}, Twhole :: Array{Float64,3})
  for ii = 1:size(T,3)
    T1 = T[:,:,ii]
    Twhole1 = Twhole[:,:,ii]
    MPI2D_DistributeCore!(p,T1,  Twhole1)
  end
  return Twhole1
end

function MPI2D_DistributeBoundary!(p :: MPI_2D_Config, T :: Array{Float64,3},  Twhole :: Array{Float64,3})
  for ii = 1:size(T,3)
    T1 = T[:,:,ii]
    Twhole1 = Twhole[:,:,ii]
    MPI2D_DistributeBoundary!(p,T1, Twhole1)
  end
  return Twhole1
end

export MPI_2D_Config
export MPI2D_GatherAll!
export MPI2D_DistributeCore!
export MPI2D_DistributeBoundary!
