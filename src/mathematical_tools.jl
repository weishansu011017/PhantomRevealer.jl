"""
Useful mathematical toolkit.
    by Wei-Shan Su,
    June 27, 2024
"""

"""
    _cart2cylin(point::Vector)
Transfer a cartesian coordinate into cylindrical coordinate.

# Parameters
- `point :: Vector`: The point with a set of cartesian coordinate.

# Returns
- `Vector`: The point with a set of cylindrical coordinate.
"""
function _cart2cylin(point::Vector)
    s = sqrt(point[1]^2 + point[2]^2)
    theta = (atan(point[2], point[1]) + 2 * pi) % (2 * pi)
    if length(point) == 2
        return [s, theta]
    else
        return [s, theta, point[3]]
    end
end

"""
    _cylin2cart(point::Vector)
Transfer a cylindrical coordinate into cartesian coordinate.

# Parameters
- `point :: Vector`: The point with a set of cylindrical coordinate..

# Returns
- `Vector`: The point with a set of cartesian coordinate
"""
function _cylin2cart(point::Vector)
    x = point[1] * cos(point[2])
    y = point[1] * sin(point[2])
    if length(point) == 2
        return [x, y]
    else
        return [x, y, point[3]]
    end
end

"""
    _Integral_1d(x::AbstractVector ,y::AbstractVector, inteval::Vector)
Intergral a discrete function f(x) in a interval in 1D. 

# Parameters
- `x :: AbstractVector`: The x value.
- `y :: AbstractVector`: The f(x) value.
- `inteval :: Vector`: The intergral range.

# Returns
- `Float64`: The integral result.
"""
function _Integral_1d(x::AbstractVector, y::AbstractVector, inteval::Vector)
    spline = CubicSplineInterpolation(x, y, extrapolation_bc = Line())
    f_interp(x) = spline(x)
    integral, error = quadgk(f_interp, inteval[1], inteval[2])
    return integral
end

"""
    value2closestvalueindex(array::AbstractVector, target::Float64)
Find the index of value which is the closest value to a given target in a array.

# Parameters
- `array::AbstractVector`: The array.
- `target::Float64`: The target for finding.

# Returns
- `Int64`: The index of closest value.
"""
function value2closestvalueindex(array::AbstractVector, target::Float64)
    target_index = argmin(abs.(target .- array))
    return target_index
end

"""
    find_array_max_index(y::AbstractVector)
Find the index of maximum value in a array.

# Parameters
- `y :: AbstractVector`: The array.

# Returns
- `Int64`: The index of maximum value.
"""
function find_array_max_index(y::AbstractVector)
    arange = 1:length(y)
    spline(x) = CubicSplineInterpolation(arange, y, extrapolation_bc = Line())(x)
    dspline(x) = ForwardDiff.derivative(spline, x)
    ddspline(x) = ForwardDiff.derivative(dspline, x)

    closest_to_zero = Inf
    target_index = 0
    for i in arange
        dspline_value = dspline(i)
        if abs(dspline_value) < closest_to_zero
            ddspline_value = ddspline(i)
            signature = sign(ddspline_value)
            if signature <= 0
                closest_to_zero = abs(dspline_value)
                target_index = i
            end
        end
    end

    if target_index == 0
        error("SearchMaxError: The maximum value is not found.")
    end

    return target_index
end
