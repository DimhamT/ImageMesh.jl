module ImageMesh

using Images
using CairoMakie
using Gridap.Helpers
using Gridap.Geometry
using Gridap.ReferenceFEs
using Gridap.Adaptivity
using Gridap.Arrays


export img2mesh

"""
    img2mesh(imgfile::String,[outfile::String]; kwargs...)

Converts an input image into a colorful simplicial mesh and saves the result to a file.

# Arguments
- `imgfile::String`: Path to the input image file.
- `outfile::String`: Path to save the output mesh image. Defaults to a derived path based on `imgfile`.
- `resolution=nothing`: Target resolution of the output image mesh. If not provided, it is calculated based on the aspect ratio of the input image.
- `base::Int`: Base scaling factor for the mesh grid. Default is 1.
- `levels::AbstractVector{Int}`: A vector of refinement levels for mesh generation. Each level determines the triangle size (smaller values produce finer meshes).
- `counts::AbstractVector{Int}`: A vector specifying how many times to apply refinement at each level. Must align with `levels`.

# Workflow
1. Loads the input image and optionally resizes it to the specified resolution.
2. Converts the image to grayscale for mesh generation.
3. Creates an initial Cartesian grid model and refines it into a simplicial mesh using the **Longest-edge bisection** algorithm.
4. Iteratively applies mesh refinement based on the `levels` and `counts` parameters.
5. Computes the color of each mesh edge based on the original image.
6. Uses CairoMakie to plot the mesh and saves the result to the output file.

# Example Usage
```julia
img2mesh("input.jpg", "output.png"; resolution=(800, 600), levels=[256, 128, 64], counts=[1, 2, 1])
```
"""
function img2mesh(
  imgfile::String,
  outfile::String=_default_output_path(imgfile); 
  resolution=nothing, 
  base::Int=1,
  levels::AbstractVector{Int}=[256, 221, 181, 141, 101, 66, 36, 16],
  counts::AbstractVector{Int}=[1,2,2,2,1,1,1,1])
  @check all( map(l->0<l<=256, levels) )
  @check sort(levels, rev=true) == levels
  @check length(levels) == length(counts)
  original = rotr90(load(imgfile))
  if isnothing(resolution)
    sz = size(original)
    r = round(Int, max(sz...)/min(sz...) * 10)
    if sz[1] > sz[2]
      resolution = (100r, 1000)
    else
      resolution = (1000, 100r)
    end
  end
  partition = broadcast((f,r)->f(r), (numerator, denominator), reduce(//, resolution))

  original = imresize(original, resolution)
  gray = convert(Matrix{Float64}, convert(Matrix{Gray}, original))
  cart = CartesianDiscreteModel((0, partition[1], 0, partition[2]), partition .* base)
  model = refine(simplexify(cart), refinement_method="nvb")

  for i = eachindex(levels)
		level = levels[i]
		for _ = 1:counts[i]
			model = img2mesh(gray, model, level/256, partition)
		end
	end

  edges = Table(get_face_coordinates(get_grid_topology(model), 1))
  seges = map(x -> Point2f(Tuple(x)), edges.data)
  bary = mapreduce(vcat, edges.data) do point
    transpose(collect(Tuple(point)))
  end
  mn = collect(size(gray) .- 1)
  mns = floor.(Int, broadcast(÷, broadcast(*, bary, mn'), collect(partition)')) .+ 1
  colors = map(eachrow(mns)) do x
    c = original[CartesianIndex(x...)]
    if c isa RGBA
      c = iszero(c.alpha) ? RGBA(0, 0, 0) : c
    elseif c isa RGB
      c = all(map(>(240/255),Tuple(c))) ? RGBA(0, 0, 0) : c 
    else
      @notimplemented
    end 
  end 

  fig = Figure(size=resolution,figure_padding=0)
  ax = CairoMakie.Axis(
    fig[1,1], 
    aspect=Rational(partition...), 
    yautolimitmargin = (0, 0), 
    xautolimitmargin = (0, 0)
  )
  hidespines!(ax)
  hidedecorations!(ax)
  linesegments!(ax, seges; color=colors)

  save(outfile, fig)
end

function img2mesh(gray::Matrix{Float64}, model, level, partition)
	MN = size(gray) .- 1
	cell_coords = get_cell_coordinates(model)
	cache = array_cache(cell_coords)
	I = Vector{Bool}(undef, length(cell_coords))
	@inbounds for i ∈ eachindex(cell_coords)
    coords = Arrays.getindex!(cache, cell_coords, i)
		bary = sum(coords) / length(coords)
		mn = @. floor(Int, bary.data * MN / partition) .+ 1
		if gray[CartesianIndex(mn)] < level
			I[i] = true
		else
			I[i] = false
		end
	end

	return refine(model, refinement_method="nvb", cells_to_refine=findall(I))
end

function _default_output_path(infile)
  path, file = splitdir(infile)
  name, _ = splitext(file)
  newfile = name * "_mesh.png"
  outfile = joinpath(path,newfile)
  return isfile(outfile) ? _default_output_path(outfile) : outfile
end

end
