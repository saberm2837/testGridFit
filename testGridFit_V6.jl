#using Pkg
#Pkg.add(url="https://github.com/WiredBrains-Lab/GridFit.jl")
using GridFit 
using StatsBase
using JSON3 
using DataFrames 
using CSV
using XLSX

function read_data(xlsx_data_fname::String,grid_dim::Dict)
    # All data is stored in an MS Excel file. This section of the code reads data from an Excel file. 
    orig_coords = Dict()
    for sub_n in [1,2,3,5,6] # These are sub#
        # Find the starting point (row, col) and endpoint (row, col) of data in each Excel sheet
        xlsx_data_start_row = 3 # default start row
        xlsx_data_start_col = 'B' # default start col
        xlsx_data_end_row = xlsx_data_start_row+(grid_dim[sub_n][1]-1)
        xlsx_data_end_col = xlsx_data_start_col+(grid_dim[sub_n][2]-1)
        # Create a string containing data location in each Excel sheet
        xlsx_data_loc = string(xlsx_data_start_col)*string(xlsx_data_start_row)*":"*string(xlsx_data_end_col)*string(xlsx_data_end_row)
        # Fetch Excel data from the i-th Excel sheet (using the XLSX package)
        grid_data = XLSX.readdata(xlsx_data_fname,string(sub_n),xlsx_data_loc)
        # Split each and every coordinate (x,y,z) and convert them into Float64 data type
        for row = 1:grid_dim[sub_n][1]
            for col = 1:grid_dim[sub_n][2] 
                grid_data[row,col] = parse.(Float64,split(grid_data[row,col], ", "))
            end
        end
        # Create a grid object using the RectGrid function (using the GridFit package)
        grid = RectGrid(grid_dim[sub_n][2],grid_dim[sub_n][1],spacing=10.0)
        # Create a dictionary from the i-th subject data (original coordinates)
        orig_coords[sub_n]=Dict(Pair.(grid.layout, grid_data))
    end
    
    # Create an internal storage to store the grid layout and coordinates in JSON file format (using the JSON3 package)
    open("grid_coords.json", "w") do f
        JSON3.pretty(f, JSON3.write(orig_coords))
    end
end

# The purpose of this function is to find all valid neighboring points (left, right, top, down)
function find_neighbors(grid::RectGrid,point::Int)
    adjacent_points = []
    nRows = size(grid.layout,1)
    nCols = size(grid.layout,2)
    
    # Check the validity of the left point 
    if(mod(point,nCols)!=0)
        push!(adjacent_points,(point+1))
    end
    
    # Check the validity of the right point 
    if(mod(point,nCols)!=1)
        push!(adjacent_points,(point-1))
    end
    
    # Check the validity of the top point 
    if(point>nCols)
        push!(adjacent_points,(point-nCols))
    end
    
    # Check the validity of the bottom point 
    if(point<=((nRows-1)*nCols))
        push!(adjacent_points,(point+nCols))
    end
    
    return adjacent_points
end

# This function helps to find random electrodes from a given set of grid points 
function select_points(grid::RectGrid,nPoints::Int=3)
    points = []
    nRows = size(grid.layout,1)
    nCols = size(grid.layout,2)
    
    # If the total number of points is less than 3 or outside the total number of electrodes, then this parameter will be set to default
    if (nPoints<3 || nPoints>nRows*nCols)
        nPoints = 3
    end
       
    # Randomly select the first point from a set of grid points
    push!(points,rand(1:nRows*nCols))

    # Seaching for the rest of the points (second, third, ...)
    adj_point_list = []
    for i=2:nPoints
        adj_points = find_neighbors(grid,points[i-1])
        for adj_p in adj_points
            if !(adj_p in points) && !(adj_p in adj_point_list)
                push!(adj_point_list,adj_p)
            end
        end
        new_point = rand(adj_point_list)
        push!(points,new_point)
        deleteat!(adj_point_list, findall(x->x==new_point,adj_point_list))
    end

    return points    
end

# This function finds the centroid from a given set of coordinates
function find_centroid(fixed_coords::Dict{Int64, Vector{Float64}})
    centroid_fixed_coords = [0.0, 0.0, 0.0]
    for point in keys(fixed_coords)
        centroid_fixed_coords += fixed_coords[point]
    end
    centroid_fixed_coords/=length(fixed_coords)
    return centroid_fixed_coords
end

# This function finds the nearest distance from a given point (coordinate) to a set of fixed coordinates  
function find_nearest_fixed_point_dist(point::Int64,orig_grid_coords::Dict{Int64,Vector{Float64}},fixed_coords::Dict{Int64, Vector{Float64}})   
    nearest_fixed_point_dist = Inf 
    if (point in keys(fixed_coords))
        nearest_fixed_point_dist = 0
    else
        for fixed_point in keys(fixed_coords)
            if (nearest_fixed_point_dist > L2dist(fixed_coords[fixed_point],orig_grid_coords[point]))
                nearest_fixed_point_dist = L2dist(fixed_coords[fixed_point],orig_grid_coords[point])
            end
        end 
    end
    return nearest_fixed_point_dist
end

# This function calculates errors to measure the performance of the optimization
# Performance analysis (Error calculation): 1. Average Error (AvgErr) 2. Maximum Error (MaxErr)
function calc_error(orig_grid_coords::Dict{Int64, Vector{Float64}},opt_grid_coords::Dict{Int64, Vector{Float64}},fixed_coords::Dict{Int64, Vector{Float64}},op::String="centroidFP")
    # Initialize variables
    # Three distance categories: short-distance (0-20 mm), mid-distance (20-40 mm), long-distance (40-60 mm), long-long-distance(60+ mm)
    err_dist_0to20 = []
    err_dist_20to40 = []
    err_dist_40to60 = []
    err_dist_60plus = []
    
    for point in keys(opt_grid_coords) 
        # Calculate the Error of i-th coordinate
        iErr=L2dist(opt_grid_coords[point],orig_grid_coords[point])
              
        # Calculate distance from either the nearest fixed coordinate or fixed coordinate centroid
        fixed_point_dist = 0.0 
        if (cmp(op,"nearestFP")==0)      
            fixed_point_dist = find_nearest_fixed_point_dist(point,orig_grid_coords,fixed_coords)
        else
            centroid_fixed_point = find_centroid(fixed_coords)
            fixed_point_dist = L2dist(centroid_fixed_point,orig_grid_coords[point]) 
        end
        
        # Sort error by distance
        if (fixed_point_dist <= 20)
            push!(err_dist_0to20,iErr)
        elseif (fixed_point_dist <= 40)
            push!(err_dist_20to40,iErr)
        elseif (fixed_point_dist <= 60)
            push!(err_dist_40to60,iErr)
        else
            push!(err_dist_60plus,iErr)
        end 
    end 

    if(length(err_dist_0to20)==0) 
        avgErr_dist_0to20 = missing
        stdErr_dist_0to20 = missing
        maxErr_dist_0to20 = missing
    else
        avgErr_dist_0to20 = round(mean(err_dist_0to20),digits=2)
        stdErr_dist_0to20 = round(std(err_dist_0to20),digits=2)
        maxErr_dist_0to20 = round(maximum(err_dist_0to20),digits=2)
    end
    if(length(err_dist_20to40)==0)
        avgErr_dist_20to40 = missing
        stdErr_dist_20to40 = missing
        maxErr_dist_20to40 = missing
    else
        avgErr_dist_20to40 = round(mean(err_dist_20to40),digits=2)
        stdErr_dist_20to40 = round(std(err_dist_20to40),digits=2)
        maxErr_dist_20to40 = round(maximum(err_dist_20to40),digits=2)
    end
    if(length(err_dist_40to60)==0)
        avgErr_dist_40to60 = missing
        stdErr_dist_40to60 = missing
        maxErr_dist_40to60 = missing
    else
        avgErr_dist_40to60 = round(mean(err_dist_40to60),digits=2)
        stdErr_dist_40to60 = round(std(err_dist_40to60),digits=2)
        maxErr_dist_40to60 = round(maximum(err_dist_40to60),digits=2)
    end
    if(length(err_dist_60plus)==0)
        avgErr_dist_60plus = missing
        stdErr_dist_60plus = missing
        maxErr_dist_60plus = missing
    else
        avgErr_dist_60plus = round(mean(err_dist_60plus),digits=2)
        stdErr_dist_60plus = round(std(err_dist_60plus),digits=2)
        maxErr_dist_60plus = round(maximum(err_dist_60plus),digits=2)
    end
    
    # Collect all errors in a dictionary object
    Errors = Dict(
       "avgErr_dist_0to20" => avgErr_dist_0to20,
        "stdErr_dist_0to20" => stdErr_dist_0to20,
        "maxErr_dist_0to20" => maxErr_dist_0to20,
        "count_dist_0to20" => length(err_dist_0to20),
        "avgErr_dist_20to40" => avgErr_dist_20to40,
        "stdErr_dist_20to40" => stdErr_dist_20to40,
        "maxErr_dist_20to40" => maxErr_dist_20to40,
        "count_dist_20to40" => length(err_dist_20to40),
        "avgErr_dist_40to60" => avgErr_dist_40to60,
        "stdErr_dist_40to60" => stdErr_dist_40to60,
        "maxErr_dist_40to60" => maxErr_dist_40to60,
        "count_dist_40to60" => length(err_dist_40to60),
        "avgErr_dist_60plus" => avgErr_dist_60plus,
        "stdErr_dist_60plus" => stdErr_dist_60plus,
        "maxErr_dist_60plus" => maxErr_dist_60plus,
        "count_dist_60plus" => length(err_dist_60plus)
    )

    return Errors
end

# Main function to run GridFit optimization
function test_grid_opt(sub_n::Int64,select_n_points::Int64,nRep::Int64)
    # The grid_matrix_dim dictionary object contains all (subject, grid matrix dimension) pairs 
    grid_matrix_dim=Dict(
        1 => (8,8),
        2 => (4,6),
        3 => (6,8),
        5 => (5,8),
        6 => (6,8)
    )
    
    # read data from the given Excel sheet
    input_fname = pwd()*"/sub1-6_grid_coords.xlsx"
    if(isfile("grid_coords.json")==false)
        read_data(input_fname,grid_matrix_dim) 
    end
    
    # Anatomy dataset filenames
    anat_dset = Dict(
        1 => "sub1_ss2.nii.gz",
        2 => "sub2_ss2.nii.gz",
        3 => "sub3_ss2.nii.gz",
        5 => "sub1_ss2.nii.gz",
        6 => "sub1_ss2.nii.gz"
    )   
         
    # Select random initial points
    grid = RectGrid(grid_matrix_dim[sub_n][2],grid_matrix_dim[sub_n][1],spacing=10.0)
    points = select_points(grid,select_n_points)
    
    # Find coordinates from selected points
    read_grid_coords = JSON3.read(pwd()*"/grid_coords.json", Dict{Int64,Dict{Int64,Vector{Float64}}})
    selected_coords = Dict{Int64,Vector{Float64}}()
    for i in points
        push!(selected_coords, i => read_grid_coords[sub_n][i])
    end
    
    # Define an empty data frame to store results    
    df_header = ["sub","points","run","avgErr_dist_0to20","stdErr_dist_0to20","maxErr_dist_0to20","count_dist_0to20","avgErr_dist_20to40","stdErr_dist_20to40","maxErr_dist_20to40","count_dist_20to40","avgErr_dist_40to60","stdErr_dist_40to60","maxErr_dist_40to60","count_dist_40to60","avgErr_dist_60plus","stdErr_dist_60plus","maxErr_dist_60plus","count_dist_60plus"] 
    Results = DataFrame([[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]], df_header)
    
    # Run the optimization
    println("Running grid opt for subject "*string(sub_n)*" with "*string(select_n_points)*" fixed points (Run # "*string(nRep)*")")
    opt = run_gridopt(grid,anat_dset[sub_n],selected_coords)

    # Calculate the error of the optimization (AvgErr and MaxErr)
    errors = calc_error(read_grid_coords[sub_n],opt.grid_coords,selected_coords)
    
    if(errors!=nothing)
        # Store results in a previously defined data frame
        push!(Results,[sub_n select_n_points nRep errors[df_header[4]] errors[df_header[5]] errors[df_header[6]] errors[df_header[7]] errors[df_header[8]] errors[df_header[9]] errors[df_header[10]] errors[df_header[11]] errors[df_header[12]] errors[df_header[13]] errors[df_header[14]] errors[df_header[15]] errors[df_header[16]] errors[df_header[17]] errors[df_header[18]] errors[df_header[19]]])  
    
        # Save fixed coordinates, optimal coordinates, and errors in three different CSV files
        CSV.write("sub"*string(sub_n)*"_pt"*string(select_n_points)*"_run"*string(nRep)*"_GridFitFixedCoords.csv",selected_coords)  
        CSV.write("sub"*string(sub_n)*"_pt"*string(select_n_points)*"_run"*string(nRep)*"_GridFitOptCoords.csv",opt.grid_coords)    
        CSV.write("sub"*string(sub_n)*"_pt"*string(select_n_points)*"_run"*string(nRep)*"_GridFitErrbyDist.csv", Results)
    end
end
